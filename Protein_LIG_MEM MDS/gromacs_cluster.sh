#!/bin/bash
#PBS -N GROMACS_Job
#PBS -l walltime=48:00:00
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -j oe
#PBS -M your.email@example.com
#PBS -m abe

# Load GROMACS module or source GROMACS environment
module load gromacs/2021.4  # Adjust the module name and version as needed
# Alternatively, if GROMACS is installed in a custom location, source the GMXRC file
# source /path/to/gromacs/bin/GMXRC

# Change to the directory from which the job was submitted
cd $PBS_O_WORKDIR

init="step5_input"
rest_prefix="step5_input"
mini_prefix="step6.0_minimization"
equi_prefix="step6.%d_equilibration"
prod_prefix="step7_production"
prod_step="step7"

# Function to check if a file exists
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "Error: File '$1' does not exist or is not accessible."
        exit 1
    fi
}

# Minimization
check_file_exists "${mini_prefix}.mdp"
check_file_exists "${init}.gro"
check_file_exists "${rest_prefix}.gro"
check_file_exists "topol.top"
check_file_exists "index.ndx"

gmx grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${init}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm ${mini_prefix}

# Equilibration
cnt=1
cntmax=6

while [ ${cnt} -le ${cntmax} ]
do
    pcnt=$((cnt - 1))
    istep=$(printf ${equi_prefix} ${cnt})
    pstep=$(printf ${equi_prefix} ${pcnt})
    if [ ${cnt} -eq 1 ]; then
        pstep=${mini_prefix}
    fi

    check_file_exists "${istep}.mdp"
    check_file_exists "${pstep}.gro"
    check_file_exists "${rest_prefix}.gro"
    check_file_exists "topol.top"
    check_file_exists "index.ndx"

    gmx grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx -maxwarn 1
    gmx mdrun -v -deffnm ${istep}
    cnt=$((cnt + 1))
done

# Production
prod_mdp="step7_production.mdp"
final_equi_gro="step6.6_equilibration.gro"

check_file_exists "${prod_mdp}"
check_file_exists "${final_equi_gro}"
check_file_exists "${rest_prefix}.gro"
check_file_exists "topol.top"
check_file_exists "index.ndx"

gmx grompp -f ${prod_mdp} -o ${prod_prefix}.tpr -c ${final_equi_gro} -r ${rest_prefix}.gro -p topol.top -n index.ndx -maxwarn 1
gmx mdrun -v -deffnm ${prod_prefix}
