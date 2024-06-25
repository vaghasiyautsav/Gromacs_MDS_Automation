# Automated GROMACS Molecular Dynamics Simulation

This repository contains scripts to automate the process of setting up and running GROMACS molecular dynamics simulations for different ligands and/or proteins.

## Initial Setup

### Step 1: Prepare MDP and Protein Files

1. Copy all sample MDP files from the `mdp_Files` folder to the respective ligand folder.
2. Place your protein PDB file (without the ligand) in the respective ligand folder. In this example, the protein file is `final_TLR7.pdb`.

## GROMACS Simulation Steps

## Auto Gromacs 

1. Run Following command in terminal to start auto gromacs script
```
bash autoGromacs.sh
```
   1. Script will ask for Starting time and end time of simulation in Nano Seconds.
   2. Then user will be asked to select folders from the list in number reprented. Selected all folders and put "d" at the end and press Enter.
   Example 
   ```
(base) vaxine@vgpu1:/media/vaxine/Expansion/Gromacs Base/Gromacs_MDS_Automation/Gromacs_MDS_Automation/Protein MDS$ bash autoGromacs.sh 
Enter the start time in nano seconds (typically 0): 0
Enter the end time in nano seconds: 1
Available folders:
1) protein1/
2) protein2/
3) Protein
4) MDP
5) Files/
6) Done
Select folders (enter numbers separated by spaces, press 'd' when done): 1 2 d

   ```
   3. And Script will do the rest. Happy Simulating !!

- I have kept two seperate bash script just in case anyone wants to do experiment on scripts. 

## Directory Structure

```
├── Parent Folder
│   ├── step1_1_groAndTopolFileGeneration.sh
│   ├── step1_2_gromacsauto.sh
│   ├── autoGromacs.sh
├── Protein1_folder
│   ├── protein1.pdb
│   ├── *.mdp
├── Protein2_folder
│   ├── protein2.pdb
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
