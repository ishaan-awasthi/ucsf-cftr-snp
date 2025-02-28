from chimerax.core.commands import run
from chimerax.core.logger import StringPlainTextLog
import re
import csv

# ========================== PARAMETERS ==========================
structures = ["6O2P"]
mutations = ["p.Val470Met"]
output_csv_path = "./distances_data.csv"

# ========================== DATA EXTRACTION ==========================

def extract_section(log_output, start_marker, end_marker):
    """
    Extract content between start and end markers in ChimeraX log output.
    """
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
    """
    Extract atom IDs and select the best atom per residue based on priority.
    """
    atom_data = []
    atom_pattern = re.compile(r'atom id (/A:(\d+)@(\S+)) idatm_type')

    for match in atom_pattern.finditer("\n".join(extract_section(log_output, "--start-atoms--", "--end-atoms--"))):
        full_id, residue_num, atom_name = match.groups()
        atom_data.append((int(residue_num), atom_name, full_id))

    # Prioritize atom selection per residue
    atom_priority = ["CA", "N", "CB", "CG1", "CG2", "CD", "CE"]
    best_atoms = {}
    
    for residue, atom_name, full_id in atom_data:
        if residue not in best_atoms:
            best_atoms[residue] = (atom_name, full_id)  # First atom by default
        else:
            current_best = best_atoms[residue][0]
            if atom_name in atom_priority and (
                current_best not in atom_priority or atom_priority.index(atom_name) < atom_priority.index(current_best)
            ):
                best_atoms[residue] = (atom_name, full_id)

    return [atom_id for _, atom_id in best_atoms.values()]  # Return optimized atom list

def extract_distances(log_output):
    """
    Extract distance values from ChimeraX log output.
    """
    print(f"[DEBUG] Raw log before extraction:\n{log_output}")  

    # Updated regex to match the actual format in logs
    distance_pattern = re.compile(r"Distance between .*?: ([\d\.]+)Å")

    matches = distance_pattern.findall("\n".join(extract_section(log_output, "--start-distance--", "--end-distance--")))

    distances = [float(match) for match in matches] if matches else []

    print(f"[DEBUG] Extracted raw distances: {distances}")  
    return distances

# ========================== STRUCTURE VALIDATION ==========================

def check_ligands():
    """
    Check if ligands exist in the opened structure before selecting them.
    """
    with StringPlainTextLog(session.logger) as log:
        run(session, "select ligand")
        run(session, "info selection")
        log_output = log.getvalue()
    return "No atoms selected" not in log_output

# ========================== DISTANCE CALCULATION ==========================

def compute_distances(atom_list1, atom_list2, structure, mutation):
    """
    Compute distances between optimized atom lists and count pairs under 5Å.
    """
    count_in_range = 0

    print(f"\n[DEBUG] Computing distances between {len(atom_list1)} and {len(atom_list2)} atoms for {structure}_{mutation}...")

    with StringPlainTextLog(session.logger) as log:
        run(session, "log text --start-distance--")

        for atom1 in atom_list1:
            for atom2 in atom_list2:
                run(session, f"distance {atom1} {atom2}")

        run(session, "log text --end-distance--")

        log_output = log.getvalue()
        distance_values = extract_distances(log_output)

        # 🔥 Debug each distance before counting
        for d in distance_values:
            print(f"[DEBUG] Checking distance {d}Å for {structure}_{mutation}")
            if d < 5.0:
                print(f"[DEBUG] ✅ Distance under 5Å: {d}Å for {structure}_{mutation}")
                count_in_range += 1

    print(f"[DEBUG] Finished computing distances for {structure}_{mutation}. Found {count_in_range} pairs under 5Å.")

    return count_in_range

# ========================== MAIN PROCESSING ==========================

def process_structure(structure, mutation):
    """
    Process a single structure-mutation pair with detailed debug logs.
    """
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

# Save results to CSV
with open(output_csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Structure_Mutation", "Pairs"])
    writer.writerows(results)

print(f"Results saved to {output_csv_path}")
print("Analysis complete.")
