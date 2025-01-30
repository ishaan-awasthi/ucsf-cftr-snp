from chimerax.core.commands import run
from chimerax.core.session import Session
from chimerax.core.logger import StringPlainTextLog
import re
import csv

# ========================== PARAMETERS: Mutations, Structures, and Weights ==========================

mutations = [
    "p.Arg31Cys", "p.Gly85Glu", "p.Ile148Thr", "p.Arg170His", "p.Glu217Gly", "p.Ser256Gly", "p.Ile269Thr",
    "p.Ile285Phe", "p.Arg297Gln", "p.Ala309Gly", "p.Ala399Val", "p.Val456Ala", "p.Leu467Phe", "p.Cys491Phe",
    "p.Ser549Asn", "p.Tyr569Asp", "p.Gly622Asp", "p.Val855Ile", "p.Tyr919Cys", "p.Asp924Asn", "p.Leu997Phe",
    "p.Ile1139Val", "p.Asn1303Ile"
]

structureNames = [
    "5UAK", "5W81", "5UAR", "6MSM", "6O1V", "6O2P", "7SVR", "7SVD", "7SV7", "8EJ1", "8EIG", "8EIQ", "8EIO", "8FZQ"
]


# Weights associated with the structures
weightages = {
    "5UAK": 1.0,
    "5W81": 0.75,
    "5UAR": 75.0,
    "6MSM": 1.0,
    "6O1V": 0.5,
    "6O2P": 0.5,
    "7SVR": 1.0,
    "7SVD": 0.5,
    "7SV7": 0.5,
    "8EJ1": 0.5,
    "8EIG": 1.0,
    "8EIQ": 0.4,
    "8EIO": 0.4,
    "8FZQ": 1.0
}

output_directory = "/Users/ishaan/Desktop/UCSF/Final/Output"
output_csv_path = "/Users/ishaan/Desktop/UCSF/Final/metrics_data.csv"

# ========================== HELPER FUNCTIONS ==========================

# Define pattern to identify backbone atoms (CA, O, N, C)
ignore_atoms_pattern = re.compile(r'\b(CA|O|N|C)\b')

# Define the list of helical ranges
list_of_ranges = [range(10, 18), range(19, 31), range(45, 65), range(68, 78), range(78, 99), range(99, 110),
                  range(116, 167), range(168, 174), range(176, 194), range(194, 202), range(202, 216),
                  range(221, 268), range(268, 277), range(279, 324), range(324, 330), range(333, 352),
                  range(352, 377), range(463, 473), range(501, 509), range(513, 526), range(526, 534),
                  range(549, 565), range(579, 596), range(606, 614), range(629, 637), range(637, 645),
                  range(845, 855), range(857, 884), range(910, 927), range(935, 958), range(960, 968),
                  range(968, 986), range(986, 1013), range(1014, 1048), range(1048, 1061), range(1061, 1070),
                  range(1070, 1122), range(1126, 1139), range(1139, 1170), range(1249, 1260), range(1278, 1284),
                  range(1299, 1306), range(1311, 1323), range(1324, 1331), range(1339, 1344), range(1347, 1363),
                  range(1371, 1377), range(1377, 1393), range(1404, 1409), range(1427, 1437)]

def get_int(string):
    '''
    Extract the first integer in an input string.
    @string: input string
    '''
    start_index = -1
    end_index = -1
    for i in range(len(string)):
        if string[i].isdigit():
            if start_index == -1:
                start_index = i
            end_index = i
    if start_index == -1:
        return 0
    else:
        integer_str = string[start_index:end_index + 1]
        return int(integer_str)

def get_int_index(string):
    '''
    Find the index of the first integer in a string.
    @string: input string
    '''
    match = re.search(r'\d', string)
    if match:
        return match.start()
    else:
        return -1

def are_in_different_ranges(residue1, residue2):
    '''
    Check if two residue numbers belong to different helical ranges.
    If one residue is within a range and the other is not, consider it interhelical.
    @residue1: first residue number
    @residue2: second residue number
    @return: True if they are in different ranges or one is in range and the other is not, False otherwise
    '''
    in_range1 = None
    in_range2 = None
    # Find which range each residue belongs to
    for r in list_of_ranges:
        if residue1 in r:
            in_range1 = r
        if residue2 in r:
            in_range2 = r
    # Consider interhelical if the residues are in different ranges or if one is in a range and the other is not
    if (in_range1 is not None and in_range2 is None) or (in_range2 is not None and in_range1 is None) or (in_range1 != in_range2):
        return True
    return False

def extract_section(html_content, start_marker, end_marker):
    '''
    Extract the content between the start and end markers.
    '''
    section_data = []
    in_section = False
    lines = html_content.splitlines()
    for line in lines:
        line = line.strip()
        if start_marker in line:
            in_section = True
            continue
        if end_marker in line:
            in_section = False
            break
        if in_section:
            section_data.append(line)
    return section_data

def count_valid_interresidual_contacts(html_content):
    '''
    Count valid interresidual contacts from the "interresidual_contacts" section of the HTML content.
    A valid interresidual contact starts with '/' and does not contain any backbone atoms.
    '''
    interresidual_contacts_data = extract_section(html_content, "--start-interresidual_contacts--", "--end-interresidual_contacts--")
    valid_interresidual_contacts_count = 0
    for line in interresidual_contacts_data:
        line = line.strip()
        # Check if the line starts with a '/'
        if line.startswith('/'):
            # Skip the line if it contains any backbone atoms
            if ignore_atoms_pattern.search(line):
                continue  # Skip lines with backbone atoms
            # If no backbone atoms, it's a valid interresidual contact
            valid_interresidual_contacts_count += 1
    return valid_interresidual_contacts_count

def count_valid_interhelical_contacts(html_content):
    '''
    Count valid interhelical contacts from the "contacts" section of the HTML content.
    A valid interhelical contact connects two residues in different helical ranges or one is not in any range.
    '''
    interhelical_contacts_data = extract_section(html_content, "--start-interhelical_contacts--", "--end-interhelical_contacts--")
    valid_interhelical_contacts_count = 0
    for line in interhelical_contacts_data:
        line = line.strip()
        # Check if the line starts with a '/'
        if line.startswith('/'):
            # Extract the two residue numbers (assumed to be at specific positions in the line)
            try:
                res1 = int(re.search(r'/A [A-Z]{3} (\d+)', line).group(1))  # First residue number
                res2 = int(re.search(r'/A [A-Z]{3} (\d+)', line[15:]).group(1))  # Second residue number
            except:
                continue
            # Check if the residues are in different ranges
            if are_in_different_ranges(res1, res2):
                valid_interhelical_contacts_count += 1
    return valid_interhelical_contacts_count

def extract_clashes(html_content):
    '''
    Extract the number of clashes between the clash markers.
    '''
    clashes_data = extract_section(html_content, "--start-clashes--", "--end-clashes--")
    for line in clashes_data:
        if 'no clashes' in line.lower():
            return 0
        clashes_match = re.search(r'(\d+) clashes', line)
        if clashes_match:
            return int(clashes_match.group(1))
    return 0

def count_valid_hydrogen_bonds(html_content, target_residue):
    '''
    Count valid hydrogen bonds involving the side chain of the target residue.
    Only counts bonds if the atom involved for the target residue is a side-chain atom (not CA, N, O, etc.)
    @html_content: log content to parse for hydrogen bonds
    @target_residue: residue number being mutated (e.g., 1303)
    '''
    hydrogen_bonds_data = extract_section(html_content, "--start-hbonds--", "--end-hbonds--")
    valid_hbond_count = 0
    pattern = rf"/A [A-Z]{{3}} {target_residue} ([A-Z0-9]+)"

    for line in hydrogen_bonds_data:
        line = line.strip()
        if line.startswith('/'):
            # Find atom following the target residue in each bond
            match = re.search(pattern, line)
            if match:
                atom = match.group(1)
                # Only count if the atom is not a backbone atom
                if not ignore_atoms_pattern.search(atom):
                    valid_hbond_count += 1

    return valid_hbond_count

def get_all_metrics(log_file, target_residue):
    '''
    Extract all structural metrics, including the number of clashes, interresidual, interhelical contacts, and H-bonds.
    @log_file: input log file
    @target_residue: residue being mutated (e.g., 1303)
    '''
    with open(log_file, 'r') as f:
        html_content = f.read()

    # Extract Clashes
    storeClashes = extract_clashes(html_content)

    # Count valid interresidual contacts
    storeInterresidualContacts = count_valid_interresidual_contacts(html_content)

    # Count valid interhelical contacts
    storeInterhelicalContacts = count_valid_interhelical_contacts(html_content)

    # Count valid hydrogen bonds
    storeHydrogenBonds = count_valid_hydrogen_bonds(html_content, target_residue)

    return [storeClashes, storeInterresidualContacts, storeInterhelicalContacts, storeHydrogenBonds]

def read_mutations(mutations):
    '''
    Process SNP data for CFTR protein from the provided list of mutations.
    '''
    data = []
    for mutation in mutations:
        selNumber = get_int(mutation)
        numberIndex = get_int_index(mutation)
        normalAminoAcid = mutation[2:numberIndex]
        mutatedAminoAcid = mutation[numberIndex + len(str(selNumber)):len(mutation)]
        data.append([str(selNumber), normalAminoAcid, mutatedAminoAcid])
    return sorted(data, key=lambda x: int(x[0]))

def get_file_name_from_path(file_path):
    '''
    Extract the structure name and residue number from the file path using regex.
    @file_path: full path to the log file
    @return: extracted file name (structure_residue)
    '''
    match = re.search(r'([a-zA-Z0-9]+)_[A-Za-z]+_(\d+)_', file_path)
    if match:
        structure_name = match.group(1)
        residue_number = match.group(2)
        return f"{structure_name}_{residue_number}"
    else:
        raise ValueError(f"File name could not be extracted from the path: {file_path}")

def get_all_structure_info(logFiles, structures, data):
    '''
    Extract all metrics (clashes, interresidual, interhelical contacts, H-bonds) before and after for all log files.
    @logFiles: list of all log files generated from get_mutation_data()
    @data: list of mutations
    @structures: list of structure names
    '''
    final_data = {}
    # Initialize the final_data dictionary to include all structure-residue pairs with amino acid info
    for structure in structures:
        for mutation in data:
            residue_number = mutation[0]
            normal_amino_acid = mutation[1]
            mutated_amino_acid = mutation[2]
            file_name = f"{structure}_{normal_amino_acid}_{residue_number}_{mutated_amino_acid}"
            final_data[file_name] = {
                'before': {'clashes': [], 'interresidual_contacts': [], 'interhelical_contacts': [], 'hbonds': []}, 
                'after': {'clashes': [], 'interresidual_contacts': [], 'interhelical_contacts': [], 'hbonds': []}
            }

    # Process each log file
    for i in range(len(logFiles)):
        if "html" in logFiles[i]:
            structure_name, residue_number = get_file_name_from_path(logFiles[i]).split("_")
            matching_mutation = next((m for m in data if m[0] == residue_number), None)
            if matching_mutation:
                normal_amino_acid = matching_mutation[1]
                mutated_amino_acid = matching_mutation[2]
                file_name = f"{structure_name}_{normal_amino_acid}_{residue_number}_{mutated_amino_acid}"
                metrics = get_all_metrics(logFiles[i], int(residue_number))

                if "MUTATED" not in logFiles[i]:
                    final_data[file_name]['before']['clashes'].append(metrics[0])
                    final_data[file_name]['before']['interresidual_contacts'].append(metrics[1])
                    final_data[file_name]['before']['interhelical_contacts'].append(metrics[2])
                    final_data[file_name]['before']['hbonds'].append(metrics[3])
                else:
                    final_data[file_name]['after']['clashes'].append(metrics[0])
                    final_data[file_name]['after']['interresidual_contacts'].append(metrics[1])
                    final_data[file_name]['after']['interhelical_contacts'].append(metrics[2])
                    final_data[file_name]['after']['hbonds'].append(metrics[3])

    return final_data

def get_mutation_data(data, structureNames):
    '''
    Create HTML files (1 for each structure, residue location, mutation status triple) that can later be parsed for structural data.
    '''
    filePaths = []
    for i in range(len(data)):
        for j in range(len(structureNames)):
            run(session, "close session")
            run(session, "open '" + structureNames[j] + "'")
            run(session, "log clear")
            x = data[i][0]
            try:
                with StringPlainTextLog(session.logger) as log:
                    run(session, "log text --start-clashes--")
                    run(session, "sel :"+x+"; clashes sel ignoreHiddenModels true select true dashes 11 radius 0.06 reveal true log true")
                    run(session, "log text --end-clashes--")
                    run(session, "log text --start-interresidual_contacts--")
                    run(session, "sel :"+x+"; contacts sel ignoreHiddenModels true select true dashes 10 radius 0.065 reveal true log true")
                    run(session, "log text --end-interresidual_contacts--")
                    run(session, "log text --start-interhelical_contacts--")
                    run(session, "sel :"+x+"; contacts sel ignoreHiddenModels true select true dashes 10 radius 0.065 reveal true log true")
                    run(session, "log text --end-interhelical_contacts--")
                    run(session, "log text --start-hbonds--")
                    run(session, "sel :"+x+"; hbonds sel color #20ffff radius 0.05 dashes 5 twoColors true slopColor #000000 select true reveal true log true")
                    run(session, "log text --end-hbonds--")
                    output = log.getvalue()
                    normal_amino_acid = data[i][1]
                    mutated_amino_acid = data[i][2]
                    file_name = f"{structureNames[j]}_{normal_amino_acid}_{x}_{mutated_amino_acid}.html"
                    with open(f"{output_directory}/{file_name}", 'w') as f:
                        print(output, file=f)
                run(session, "log clear")
                run(session, "close session")
                run(session, "open '" + structureNames[j] + "'")
                run(session, "log clear")
                y = data[i][2]
                with StringPlainTextLog(session.logger) as log: 
                    run(session, "swapaa #!1/A:"+x+" "+y.upper()+" criteria 1 rotLib Dunbrack")
                    run(session, "log text --start-clashes--")
                    run(session, "sel :"+x+"; clashes sel ignoreHiddenModels true select true dashes 11 radius 0.06 reveal true log true")
                    run(session, "log text --end-clashes--")
                    run(session, "log text --start-interresidual_contacts--")
                    run(session, "sel :"+x+"; contacts sel ignoreHiddenModels true select true dashes 10 radius 0.065 reveal true log true")
                    run(session, "log text --end-interresidual_contacts--")
                    run(session, "log text --start-interhelical_contacts--")
                    run(session, "sel :"+x+"; contacts sel ignoreHiddenModels true select true dashes 10 radius 0.065 reveal true log true")
                    run(session, "log text --end-interhelical_contacts--")
                    run(session, "log text --start-hbonds--")
                    run(session, "sel :"+x+"; hbonds sel color #20ffff radius 0.05 dashes 5 twoColors true slopColor #000000 select true reveal true log true")
                    run(session, "log text --end-hbonds--")
                    output = log.getvalue()
                    mutated_file_name = f"{structureNames[j]}_{normal_amino_acid}_{x}_{mutated_amino_acid}_MUTATED.html"
                    with open(f"{output_directory}/{mutated_file_name}", 'w') as f:
                        print(output, file=f)
                filePaths.append(f"{output_directory}/{file_name}")
                filePaths.append(f"{output_directory}/{mutated_file_name}")
                print("filepaths made this round")
                run(session, "log clear")
            except:
                filePaths.append(f"{structureNames[j]}.cif")
                filePaths.append(f"{structureNames[j]}.cif")
        run(session, "close session")
    return filePaths

# ========================== MAIN CODE ==========================

data = read_mutations(mutations)
logFiles = get_mutation_data(data, structureNames)
final_data = get_all_structure_info(logFiles, structureNames, data)

# Save the data to CSV, including clashes, interresidual contacts, interhelical contacts, and hydrogen bonds (before and after)
with open(output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    header = ["Structure_Residue", "Before_Clashes", "After_Clashes", "Before_InterresidualContacts", "After_InterresidualContacts", "Before_InterhelicalContacts", "After_InterhelicalContacts", "Before_Hbonds", "After_Hbonds", "Delta_Clashes", "Delta_Interresidual", "Delta_Interhelical", "Delta_Hbonds"]
    writer.writerow(header)

    for key, value in final_data.items():
        before_clashes = value['before']['clashes'][0] if value['before']['clashes'] else -1
        after_clashes = value['after']['clashes'][0] if value['after']['clashes'] else -1
        before_interresidual_contacts = value['before']['interresidual_contacts'][0] if value['before']['interresidual_contacts'] else -1
        after_interresidual_contacts = value['after']['interresidual_contacts'][0] if value['after']['interresidual_contacts'] else -1
        before_interhelical_contacts = value['before']['interhelical_contacts'][0] if value['before']['interhelical_contacts'] else -1
        after_interhelical_contacts = value['after']['interhelical_contacts'][0] if value['after']['interhelical_contacts'] else -1
        before_hbonds = value['before']['hbonds'][0] if value['before']['hbonds'] else -1
        after_hbonds = value['after']['hbonds'][0] if value['after']['hbonds'] else -1
        
        # Calculate deltas (after - before)
        delta_clashes = after_clashes - before_clashes if before_clashes != -1 and after_clashes != -1 else "N/A"
        delta_interresidual = after_interresidual_contacts - before_interresidual_contacts if before_interresidual_contacts != -1 and after_interresidual_contacts != -1 else "N/A"
        delta_interhelical = after_interhelical_contacts - before_interhelical_contacts if before_interhelical_contacts != -1 and after_interhelical_contacts != -1 else "N/A"
        delta_hbonds = after_hbonds - before_hbonds if before_hbonds != -1 and after_hbonds != -1 else "N/A"
        
        # Write the row with the delta values
        row = [key, before_clashes, after_clashes, before_interresidual_contacts, after_interresidual_contacts, before_interhelical_contacts, after_interhelical_contacts, before_hbonds, after_hbonds, delta_clashes, delta_interresidual, delta_interhelical, delta_hbonds]
        writer.writerow(row)

print(f"Data saved to {output_csv_path}")