#!/bin/bash

# Function to display folders and get user selection
select_folders() {
    echo "Available folders:"
    options=($(ls -d */))
    for i in "${!options[@]}"; do
        echo "$((i+1))) ${options[$i]}"
    done
    echo "$(( ${#options[@]} + 1 ))) Done"

    selected_folders=()
    while true; do
        read -p "Select folders (enter numbers separated by spaces, press 'd' when done): " -a selections
        for selection in "${selections[@]}"; do
            if [[ "$selection" == "d" ]]; then
                return
            elif [[ "$selection" -ge 1 && "$selection" -le "${#options[@]}" ]]; then
                folder="${options[$((selection-1))]}"
                if [[ ! " ${selected_folders[@]} " =~ " ${folder} " ]]; then
                    selected_folders+=("$folder")
                    echo "Selected: $folder"
                else
                    echo "Already selected: $folder"
                fi
            else
                echo "Invalid option: $selection"
            fi
        done
    done
}

# Function to run the command in each selected folder
run_command_in_folders() {
    for folder in "${selected_folders[@]}"; do
        echo "Processing folder: $folder"
        cd "$folder" || continue

        # Find the .pdb file in the current folder
        pdb_file=$(find . -maxdepth 1 -name "*.pdb" | head -n 1)
        if [[ -z "$pdb_file" ]]; then
            echo "No .pdb file found in $folder"
            cd ..
            continue
        fi

        echo "Using PDB file: $pdb_file"
        echo -e "8\n1" | gmx pdb2gmx -f "$pdb_file" -ignh

        cd ..
    done
}

# Main script execution
select_folders
run_command_in_folders
