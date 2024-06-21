# Automated GROMACS Molecular Dynamics Simulation

This repository contains scripts to automate the process of setting up and running GROMACS molecular dynamics simulations for different ligands and/or proteins.

## Initial Setup

### Step 1: Prepare MDP and Protein Files

1. Copy all sample MDP files from the `mdp_Files` folder to the respective ligand folder.
2. Place your protein PDB file (without the ligand) in the respective ligand folder. In this example, the protein file is `final_TLR7.pdb`.

### Step 2: Prepare Ligand Files

1. Rename your ligand `.mol2` file to `LIG.mol2`.
2. Submit `LIG.mol2` to the [SwissParam](http://www.swissparam.ch/) website to generate the parameter files for the ligand.

### Step 3: Extract Ligand Parameter Files

1. Download the `LIG.zip` file from SwissParam.
2. Extract `LIG.zip` in the same directory as the respective ligand.
3. Copy `LIG.itp` and `LIG.pdb` from the extracted files to the working directory.

## GROMACS Simulation Steps

### Step 1: Generate GRO and Topology Files

1. Update the PDB file name in `step1_1_groAndTopolFileGeneration.sh` to match your protein file (e.g., `final_TLR7.pdb`).
2. Open the current folder in a terminal.
3. Run the following command:
   ```bash
   bash step1_1_groAndTopolFileGeneration.sh
   ```

### Step 2: Update Topology Files

1. If you are using the latest GROMACS and Python in your Ubuntu system, no changes are needed in this step.
2. Run the following command:
   ```bash
   python3 step1_2_topolfileupdate.py
   ```

### Step 3: Update Configuration File

1. This script updates the GROMACS-generated `conf.gro` file of the protein by adding the GROMACS-generated `lig.gro` file of the ligand. After this step, `conf.gro` will represent the protein-ligand complex.
2. Run the following command:
   ```bash
   python3 step1_3_edit_confgrofile.py
   ```

### Step 4: Run GROMACS Simulation

1. This bash script contains the code to run GROMACS commands sequentially to perform the simulation production.
2. Edit `md.mdp` if you have any specific requirements for the production run.
3. Run the following command:
   ```bash
   bash step1_4_gromacsauto.sh
   ```

## Directory Structure

```
├── mdp_Files
│   ├── *.mdp
├── protein_folder
│   ├── final_TLR7.pdb
│   ├── LIG.itp
│   ├── LIG.pdb
│   ├── step1_1_groAndTopolFileGeneration.sh
│   ├── step1_2_topolfileupdate.py
│   ├── step1_3_edit_confgrofile.py
│   ├── step1_4_gromacsauto.sh
```

## Notes

- Ensure that the names of the files and folders match the expected names in the scripts to avoid errors.
- Modify the provided scripts as needed to fit your specific requirements.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [SwissParam](http://www.swissparam.ch/) for ligand parameter files.
- [GROMACS](http://www.gromacs.org/) for molecular dynamics simulations.

For further questions or contributions, please open an issue or submit a pull request.
