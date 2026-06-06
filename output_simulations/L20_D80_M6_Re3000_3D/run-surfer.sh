#!/bin/sh
#PBS -N L20_D80_M6_Re3000_3D_surfer
#PBS -q mafat14_512_q
#PBS -M felipe.aguir@campus.technion.ac.il
#PBS -mbea
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l place=scatter

cd $PBS_O_WORKDIR # change to directory where qsub is issued

# Display the hostname of the MPI head node
echo "Job is submitted from \"Head Node\": $(hostname)"

echo running on nodes
cat $PBS_NODEFILE | sort -u > machines
echo "Hostfile content (unique nodes):"
cat machines

# source bashrc (for charles_launch.sh)
source ~/.bashrc

# Use with Intel MPI 2021:
source /usr/local/intel21/setup.sh

# -np should be select*ncpus

# surfer STL->sbin
charles_launch.sh -np 1 -mpi impi -mpi_dir /usr/local/intel21/mpi/2021.3.0 -target cpu surfer.exe -i surfer.in > ./logs/surfer.log

# Clean up
rm -f machines

echo "CharLES Launch in $(hostname) FINISHED!!!"
