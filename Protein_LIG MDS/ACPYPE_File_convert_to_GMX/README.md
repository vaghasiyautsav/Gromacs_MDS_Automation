# Change Name to LIG Script

This Python script is designed to help users rename and modify specific files within selected folders. The script performs the following tasks:
1. Prompts the user to select folders from the current directory.
2. Renames `.itp` and `.pdb` files to `LIG.itp` and `LIG.pdb`, respectively.
3. Modifies the content of `LIG.itp` and `LIG.pdb` files:
   - Replaces the ligand name in the `[ moleculetype ]` section with `LIG`.
   - Replaces all occurrences of `UNK` with `LIG`.

## Directory Structure

Ensure your directory structure looks like this before running the script:

```
/path/to/your/directory/
│
├── folder1/
│   ├── somefile.itp
│   ├── somefile.pdb
│   └── other files...
│
├── folder2/
│   ├── anotherfile.itp
│   ├── anotherfile.pdb
│   └── other files...
│
└── change_name_LIG.py
```

## Prerequisites

- Python 3.x installed on your system.

## Usage

1. **Clone or download the script** to your local machine.

2. **Navigate to the directory** containing the script and the folders you want to process:
   ```bash
   cd /path/to/your/directory
   ```

3. **Run the script**:
   ```bash
   python3 change_name_LIG.py
   ```

4. **Follow the prompts** to select the folders you want to process. Enter the numbers corresponding to the folders, separated by spaces, and press 'd' when done.

## Example

```
Available folders:
1) folder1/
2) folder2/
3) Done
Select folders (enter numbers separated by spaces, press 'd' when done): 1 2 d
Selected: folder1/
Selected: folder2/
Processing folder: folder1/
Processing folder: folder2/
```

## Script Details

The script performs the following steps:

1. **List Folders**: Lists all directories in the current path.
2. **Select Folders**: Prompts the user to select folders by entering their corresponding numbers.
3. **Process Files**: Iterates over the selected folders, renames `.itp` and `.pdb` files to `LIG.itp` and `LIG.pdb` respectively, and edits the contents of these files:
   - Replaces the ligand name in the `[ moleculetype ]` section with `LIG`.
   - Replaces all occurrences of `UNK` with `LIG`.

## Notes

- Ensure that the `.itp` and `.pdb` files exist in the selected folders.
- The script assumes that the ligand name in the `[ moleculetype ]` section of the `.itp` file ends with `.pdb.mol2`.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Thanks to the GROMACS community for their tools and resources.

```

This README file provides clear instructions on how to use the script, what the directory structure should look like, and details about the script's functionality. It is formatted to look good on GitHub.
