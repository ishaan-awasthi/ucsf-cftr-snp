from chimerax.core.commands import run
from chimerax.core.logger import StringPlainTextLog
import re
import csv

# ========================== PARAMETERS ==========================
mutations = [
    "p.Arg31Cys", "p.Gly85Glu", "p.Ile148Thr", "p.Arg170His", "p.Glu217Gly", "p.Ser256Gly", "p.Ile269Thr",
    "p.Ile285Phe", "p.Arg297Gln", "p.Ala309Gly", "p.Ala399Val", "p.Val456Ala", "p.Leu467Phe", "p.Cys491Phe",
    "p.Ser549Asn", "p.Tyr569Asp", "p.Gly622Asp", "p.Val855Ile", "p.Tyr919Cys", "p.Asp924Asn", "p.Leu997Phe",
    "p.Ile1139Val", "p.Asn1303Ile"
]

structures = [
    "5UAK", "5W81", "5UAR", "6MSM", "6O1V", "6O2P", "7SVR", "7SVD", "7SV7", "8EJ1", "8EIG", "8EIQ", "8EIO", "8FZQ"
]
output_csv_path = "./distances_data.csv"

# ========================== DATA EXTRACTION ==========================

def extract_section(log_output, start_marker, end_marker):
    section_data = []
    in_section = False
    for line in log_output.splitlines():
        if start_marker in line:
            in_section = True
            continue
        if end_marker in line:
            break
        if in_section:
            section_data.append(line.strip())
    return section_data

def extract_atoms(log_output):
    atom_data = []
    atom_pattern = re.compile(r'atom id (/A:(\d+)@(\S+)) idatm_type')

    for match in atom_pattern.finditer("\n".join(extract_section(log_output, "--start-atoms--", "--end-atoms--"))):
        full_id, residue_num, atom_name = match.groups()
        atom_data.append((int(residue_num), atom_name, full_id))

    atom_priority = ["CA", "N", "CB", "CG1", "CG2", "CD", "CE"]
    best_atoms = {}
    
    for residue, atom_name, full_id in atom_data:
        if residue not in best_atoms:
            best_atoms[residue] = (atom_name, full_id)
        else:
            current_best = best_atoms[residue][0]
            if atom_name in atom_priority and (
                current_best not in atom_priority or atom_priority.index(atom_name) < atom_priority.index(current_best)
            ):
                best_atoms[residue] = (atom_name, full_id)

    return [atom_id for _, atom_id in best_atoms.values()]

def extract_distances(log_output):
    print(f"[DEBUG] Raw log before extraction:\n{log_output}")

    distance_pattern = re.compile(r"Distance between .*?: ([\d\.]+)Å")
    matches = distance_pattern.findall("\n".join(extract_section(log_output, "--start-distance--", "--end-distance--")))

    distances = [float(match) for match in matches] if matches else []

    print(f"[DEBUG] Extracted raw distances: {distances}")
    return distances

# ========================== STRUCTURE VALIDATION ==========================

def check_ligands():
    with StringPlainTextLog(session.logger) as log:
        run(session, "select ligand")
        run(session, "info selection")
        log_output = log.getvalue()
    return "No atoms selected" not in log_output

# ========================== DISTANCE CALCULATION ==========================

def compute_distances(atom_list1, atom_list2, structure, mutation):
    count_in_range = 0

    print(f"\n[DEBUG] Computing distances between {len(atom_list1)} and {len(atom_list2)} atoms for {structure}_{mutation}...")

    if not atom_list1 or not atom_list2:
        print(f"[ERROR] One of the atom lists is empty! Skipping distance calculations for {structure}_{mutation}.")
        return 0

    # Delete existing distance measurements to avoid conflicts
    run(session, "distance delete all")

    with StringPlainTextLog(session.logger) as log:
        run(session, "log text --start-distance--")

        measured_pairs = set()  # To keep track of already measured distances

        for atom1 in atom_list1:
            for atom2 in atom_list2:
                if atom1 == atom2 or (atom1, atom2) in measured_pairs or (atom2, atom1) in measured_pairs:
                    print(f"[DEBUG] Skipping duplicate/self-distance: {atom1} to {atom2}")
                    continue

                print(f"[DEBUG] Measuring distance: {atom1} to {atom2}")
                run(session, f"distance {atom1} {atom2}")
                measured_pairs.add((atom1, atom2))  # Store the pair to prevent duplicate calculations

        run(session, "log text --end-distance--")

        log_output = log.getvalue()
        distance_values = extract_distances(log_output)

        for d in distance_values:
            print(f"[DEBUG] Checking distance {d}Å for {structure}_{mutation}")
            if d < 5.0:
                print(f"[DEBUG] ✅ Distance under 5Å: {d}Å for {structure}_{mutation}")
                count_in_range += 1

    print(f"[DEBUG] Finished computing distances for {structure}_{mutation}. Found {count_in_range} pairs under 5Å.")
    return count_in_range


# ========================== MAIN PROCESSING ==========================

def process_structure(structure, mutation):
    print(f"\n========== Processing {structure} with mutation {mutation} ==========")
    residue_number = re.search(r'\d+', mutation).group()

    print(f"[DEBUG] Opening structure {structure}...")
    run(session, f"open {structure}")

    print(f"[DEBUG] Checking for ligands in {structure}...")
    if not check_ligands():
        print(f"[DEBUG] No ligands found in {structure}. Skipping...")
        run(session, "close all")
        return 0  

    print(f"[DEBUG] Selecting atoms within 4.5Å of ligand...")
    try:
        with StringPlainTextLog(session.logger) as log:
            run(session, "log text --start-atoms--")
            run(session, "select zone ligand 4.5 protein residues true;")
            run(session, "info atoms sel")
            run(session, "log text --end-atoms--")
            log_output = log.getvalue()

        if "No atoms selected" in log_output:
            print(f"[DEBUG] No atoms found in 4.5Å ligand zone for {structure}. Skipping...")
            run(session, "close all")
            return 0  

        atom_list1 = extract_atoms(log_output)
        print(f"[DEBUG] Extracted {len(atom_list1)} atoms in 4.5Å zone of ligand.")

    except Exception as e:
        print(f"[ERROR] Ligand zone selection failed for {structure}: {str(e)}")
        run(session, "close all")
        return 0  

    print(f"[DEBUG] Deselecting all atoms...")
    run(session, "~sel")

    print(f"[DEBUG] Selecting atoms near residue {residue_number} in {structure}...")
    try:
        run(session, f"sel :{residue_number}")

        with StringPlainTextLog(session.logger) as log:
            run(session, "log text --start-atoms--")
            run(session, "select sel@<5")
            run(session, "info atoms sel")
            run(session, "log text --end-atoms--")
            log_output = log.getvalue()

        if "No atoms selected" in log_output:
            print(f"[DEBUG] No atoms found near mutation site {residue_number} for {structure}. Skipping...")
            run(session, "close all")
            return 0  

        atom_list2 = extract_atoms(log_output)
        print(f"[DEBUG] Extracted {len(atom_list2)} atoms in 5Å zone of residue {residue_number}.")

    except Exception as e:
        print(f"[ERROR] Mutation site selection failed for {structure}: {str(e)}")
        run(session, "close all")
        return 0  

    print(f"[DEBUG] Starting distance computation between {len(atom_list1)} and {len(atom_list2)} atoms...")
    result = compute_distances(atom_list1, atom_list2, structure, mutation)
    print(f"[DEBUG] Finished processing {structure} {mutation}. Found {result} pairs under 5Å.")

    run(session, "close all")

    return result

# ========================== EXECUTION & CSV OUTPUT ==========================

print("Starting analysis...")
results = []

for structure in structures:
    for mutation in mutations:
        count = process_structure(structure, mutation)
        results.append([f"{structure}_{mutation}", count])
        print(f"{structure} {mutation}: {count} atom pairs under 5Å")

with open(output_csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Structure_Mutation", "Pairs"])
    writer.writerows(results)

print(f"Results saved to {output_csv_path}")
print("Analysis complete.")
