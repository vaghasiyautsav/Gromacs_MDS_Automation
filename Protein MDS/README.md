# Automated GROMACS Molecular Dynamics Simulation

This repository contains scripts to automate the process of setting up and running GROMACS molecular dynamics simulations for different ligands and/or proteins.

## Initial Setup

### Step 1: Prepare MDP and Protein Files

1. Copy all sample MDP files from the `mdp_Files` folder to the respective ligand folder.
2. Place your protein PDB file (without the ligand) in the respective ligand folder. In this example, the protein file is `final_TLR7.pdb`.

## GROMACS Simulation Steps

### Step 1: Generate GRO and Topology Files

1. Update the PDB file name in `step1_1_groAndTopolFileGeneration.sh` to match your protein file (e.g., `final_TLR7.pdb`).
2. Open the current folder in a terminal.
3. Run the following command:
   ```bash
   bash step1_1_groAndTopolFileGeneration.sh
   ```
4. Script will ask user to give the name of protein which is present in parent directory as well as in each folders.

### Step 4: Run GROMACS Simulation

1. This bash script contains the code to run GROMACS commands sequentially to perform the simulation production.
2. Edit `md.mdp` if you have any specific requirements for the production run.
3. Run the following command:
   ```bash
   bash step1_2_gromacsauto.sh
   ```

## Directory Structure

```
├── Parent Folder
│   ├── step1_1_groAndTopolFileGeneration.sh
│   ├── step1_2_gromacsauto.sh
├── Protein_folder
│   ├── final_TLR7.pdb
│   ├── *.mdp
```

## Notes

- Ensure that the names of the files and folders match the expected names in the scripts to avoid errors.
- Modify the provided scripts as needed to fit your specific requirements.

## License

This project is licensed.

## Acknowledgments

- [SwissParam](http://www.swissparam.ch/) for ligand parameter files.
- [GROMACS](http://www.gromacs.org/) for molecular dynamics simulations.

For further questions or contributions, please open an issue or submit a pull request.
