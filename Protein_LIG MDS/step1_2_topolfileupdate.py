import os

def modify_topol_file(file_path):
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Check if the lines are already present
    include_lig_itp = '#include "LIG.itp"\n'
    lig_molecule = 'LIG                 1\n'

    if include_lig_itp not in lines:
        # Find the position to insert '#include "LIG.itp"'
        for i, line in enumerate(lines):
            if line.strip() == '#include "charmm36-jul2022.ff/forcefield.itp"':
                lines.insert(i + 1, include_lig_itp)
                break

    if lig_molecule not in lines:
        # Find the position to insert 'LIG                 1'
        for i, line in enumerate(lines):
            if line.strip() == 'Protein_chain_B     1':
                lines.insert(i + 1, lig_molecule)
                break

    # Write the modified contents back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

def process_folders():
    current_directory = os.getcwd()
    for root, dirs, files in os.walk(current_directory):
        for file in files:
            if file == 'topol.top':
                file_path = os.path.join(root, file)
                print(f"Modifying file: {file_path}")
                modify_topol_file(file_path)

if __name__ == "__main__":
    process_folders()

