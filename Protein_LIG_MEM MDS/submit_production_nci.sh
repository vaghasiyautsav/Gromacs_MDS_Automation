#!/bin/bash

# Generated by CHARMM-GUI (http://www.charmm-gui.org) v3.7

# This folder contains GROMACS formatted force fields, a pre-optimized PDB structure, and GROMACS inputs.
# All input files were optimized for GROMACS 2019.2 or above, so lower version of GROMACS can cause some errors.
# We adopted the Verlet cut-off scheme for all minimization, equilibration, and production steps because it is 
# faster and more accurate than the group scheme. If you have a trouble with a performance of Verlet scheme while 
# running parallelized simulation, you should check if you are using appropriate command line.
# For MPI parallelizing, we recommend the following command:
# mpirun -np $NUM_CPU gmx mdrun -ntomp 1

# Variables for different steps
equi_prefix="step6.%d_equilibration"
prod_prefix="step7_production"
prod_step="step7"

# Production
cnt=1
cntmax=10

while [ ${cnt} -le ${cntmax} ]; do
    pcnt=$((cnt - 1))
    istep=${prod_step}_${cnt}
    pstep=${prod_step}_${pcnt}

    if [ ${cnt} -eq 1 ]; then
        pstep=$(printf ${equi_prefix} 6)
        gmx grompp -f ${prod_prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -p topol.top -n index.ndx
    else
        gmx grompp -f ${prod_prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -t ${pstep}.cpt -p topol.top -n index.ndx
    fi

    # Submit production job to NCI cluster
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=md_production_${cnt}
#SBATCH --output=md_production_${cnt}.out
#SBATCH --error=md_production_${cnt}.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=compute

module load gromacs/2019.2

gmx mdrun -v -deffnm ${istep}
EOF

    cnt=$((cnt + 1))
done