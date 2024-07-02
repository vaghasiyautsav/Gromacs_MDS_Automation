import os
import shutil
import re

def list_folders():
    folders = [f for f in os.listdir() if os.path.isdir(f)]
    return folders

def select_folders(folders):
    print("Available folders:")
    for i, folder in enumerate(folders, 1):
        print(f"{i}) {folder}")
    print(f"{len(folders) + 1}) Done")

    selected_folders = []
    while True:
        selections = input("Select folders (enter numbers separated by spaces, press 'd' when done): ").split()
        for selection in selections:
            if selection.lower() == 'd':
                return selected_folders
            elif selection.isdigit() and 1 <= int(selection) <= len(folders):
                folder = folders[int(selection) - 1]
                if folder not in selected_folders:
                    selected_folders.append(folder)
                    print(f"Selected: {folder}")
                else:
                    print(f"Already selected: {folder}")
            else:
                print(f"Invalid option: {selection}")

def process_files(folders):
    for folder in folders:
        print(f"Processing folder: {folder}")
        os.chdir(folder)

        # Rename .itp and .pdb files
        for file in os.listdir():
            if file.endswith(".itp"):
                shutil.move(file, "LIG.itp")
            elif file.endswith(".pdb"):
                shutil.move(file, "LIG.pdb")

        # Edit LIG.itp file
        if os.path.exists("LIG.itp"):
            with open("LIG.itp", "r") as file:
                content = file.read()
            
            # Find the ligand name in the [ moleculetype ] section
            ligand_name_match = re.search(r'^\s*([^\s]+\.pdb\.mol2)\s+\d+', content, re.MULTILINE)
            if ligand_name_match:
                ligand_name = ligand_name_match.group(1)
                content = content.replace(ligand_name, "LIG")
            
            # Replace all occurrences of UNK with LIG
            content = content.replace("UNK", "LIG")
            
            with open("LIG.itp", "w") as file:
                file.write(content)

        # Edit LIG.pdb file
        if os.path.exists("LIG.pdb"):
            with open("LIG.pdb", "r") as file:
                content = file.read()
            content = content.replace("UNK", "LIG")
            with open("LIG.pdb", "w") as file:
                file.write(content)

        os.chdir("..")

def main():
    folders = list_folders()
    selected_folders = select_folders(folders)
    process_files(selected_folders)

if __name__ == "__main__":
    main()
