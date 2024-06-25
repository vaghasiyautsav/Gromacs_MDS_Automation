import os

def modify_gro_files(folder, already_processed_folders):
    lig_gro_path = os.path.join(folder, 'lig.gro')
    conf_gro_path = os.path.join(folder, 'conf.gro')

    # Read lig.gro and extract lines from line 3 to the second last line
    with open(lig_gro_path, 'r') as lig_file:
        lig_lines = lig_file.readlines()
        lig_data = lig_lines[2:-1]  # Extract lines from line 3 to the second last line
        num_lig_lines = len(lig_data)

    # Read conf.gro and check if ligand details are already present
    with open(conf_gro_path, 'r') as conf_file:
        conf_lines = conf_file.readlines()

    # Check if ligand details are already present
    if any('LIG' in line for line in conf_lines):
        already_processed_folders.append(folder)
        return

    # Insert lig_data before the last line
    conf_lines = conf_lines[:-1] + lig_data + conf_lines[-1:]

    # Calculate the new number of atoms
    num_conf_atoms = int(conf_lines[1].strip())  # Get the number from the second line
    new_num_atoms = num_conf_atoms + num_lig_lines

    # Update the number of atoms in the second line
    conf_lines[1] = f"{new_num_atoms}\n"

    # Write the modified conf.gro back to the file
    with open(conf_gro_path, 'w') as conf_file:
        conf_file.writelines(conf_lines)

def process_folders():
    current_directory = os.getcwd()
    already_processed_folders = []
    for root, dirs, files in os.walk(current_directory):
        for dir_name in dirs:
            folder_path = os.path.join(root, dir_name)
            lig_gro_path = os.path.join(folder_path, 'lig.gro')
            conf_gro_path = os.path.join(folder_path, 'conf.gro')
            if os.path.exists(lig_gro_path) and os.path.exists(conf_gro_path):
                print(f"Processing folder: {folder_path}")
                modify_gro_files(folder_path, already_processed_folders)

    if already_processed_folders:
        print("Ligand details already present in the following folders. Please run the next step: step1_4_gromacsauto.sh")
        for folder in already_processed_folders:
            print(folder)

if __name__ == "__main__":
    process_folders()

