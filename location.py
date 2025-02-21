from chimerax.core.commands import run
from chimerax.core.logger import StringPlainTextLog
import re
import csv

# ========================== PARAMETERS ==========================
structures = ["5UAK", "5W81", "5UAR", "6MSM", "6O1V", "6O2P", "7SVR", "7SVD", "7SV7", "8EJ1", "8EIG", "8EIQ", "8EIO", "8FZQ"]
mutations = ["p.Arg31Cys", "p.Gly85Glu", "p.Ile148Thr", "p.Arg170His", "p.Glu217Gly"]  # Add full list here
output_csv_path = "distances_data.csv"

def extract_section(log_output, start_marker, end_marker):
    """
    Extract content between start and end markers in ChimeraX log output.
    """
    section_data = []
    in_section = False
    lines = log_output.splitlines()
    for line in lines:
        if start_marker in line:
            in_section = True
            continue
        if end_marker in line:
            in_section = False
            break
        if in_section:
            section_data.append(line.strip())
    return section_data

def extract_atoms(log_output):
    """
    Extract atom IDs from a marked section of ChimeraX log output.
    """
    atom_pattern = re.compile(r'atom id (/A:\d+@\S+) idatm_type')
    return [match.group(1) for match in atom_pattern.finditer("\n".join(extract_section(log_output, "--start-atoms--", "--end-atoms--")))]

def extract_distances(log_output):
    """
    Extract distance values from ChimeraX log output.
    """
    distance_pattern = re.compile(r'distance.*? ([\d\.]+)A')
    return [float(match.group(1)) for match in distance_pattern.finditer("\n".join(extract_section(log_output, "--start-distance--", "--end-distance--")))]

def check_ligands():
    """
    Check if ligands exist in the opened structure before selecting them.
    """
    with StringPlainTextLog(session.logger) as log:
        run(session, "select ligand")
        run(session, "info selection")
        log_output = log.getvalue()
    return "No atoms selected" not in log_output

def compute_distances(atom_list1, atom_list2):
    """
    Compute distances between atoms in two lists and count how many fall in the 3-4Å range.
    """
    count_in_range = 0
    print("Computing distances...")
    with StringPlainTextLog(session.logger) as log:
        run(session, "log text --start-distance--")
        for atom1 in atom_list1:
            for atom2 in atom_list2:
                run(session, f"distance {atom1} {atom2}")
        run(session, "log text --end-distance--")
        distance_values = extract_distances(log.getvalue())
    
    count_in_range = sum(1 for d in distance_values if 3.0 <= d <= 4.0)
    print(f"Finished computing distances. Found {count_in_range} pairs in range.")
    return count_in_range

def process_structure(structure, mutation):
    """
    Process a single structure-mutation pair.
    """
    print(f"Processing {structure} with mutation {mutation}...")
    residue_number = re.search(r'\d+', mutation).group()
    run(session, f"open {structure}")
    
    # Check if a ligand exists before proceeding
    if not check_ligands():
        print(f"Skipping {structure} - No ligand found.")
        return 0
    
    with StringPlainTextLog(session.logger) as log:
        run(session, "log text --start-atoms--")
        run(session, "select zone ligand 4.5 protein residues true;")
        run(session, "info atoms sel")
        run(session, "log text --end-atoms--")
        atom_list1 = extract_atoms(log.getvalue())
    
    print(f"Extracted {len(atom_list1)} atoms in 4.5A zone of ligand.")
    
    run(session, "~sel")
    run(session, f"sel :{residue_number}")
    
    with StringPlainTextLog(session.logger) as log:
        run(session, "log text --start-atoms--")
        run(session, "select sel@<5")
        run(session, "info atoms sel")
        run(session, "log text --end-atoms--")
        atom_list2 = extract_atoms(log.getvalue())
    
    print(f"Extracted {len(atom_list2)} atoms in 5A zone of residue {residue_number}.")
    
    return compute_distances(atom_list1, atom_list2)

print("Starting analysis...")
results = []

for structure in structures:
    for mutation in mutations:
        count = process_structure(structure, mutation)
        results.append([structure, mutation, count])
        print(f"{structure} {mutation}: {count} atom pairs in range")

with open(output_csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Structure", "Mutation", "Pairs in 3-4Å range"])
    writer.writerows(results)
print(f"Results saved to {output_csv_path}")
print("Analysis complete.")
