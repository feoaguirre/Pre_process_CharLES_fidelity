#!/bin/sh
#PBS -N L20_D80_M6_Re3000_3D_steady
#PBS -q mafat14_512_q
#PBS -M felipe.aguir@campus.technion.ac.il
#PBS -mbea
#PBS -l select=1:ncpus=512:mpiprocs=512
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

# charles
charles_launch.sh -np 512 -mpi impi -mpi_dir /usr/local/intel21/mpi/2021.3.0 -target cpu --hostfile machines -perhost 512 charles_ig.exe -i steady_charles_ig.in > ./logs/charles_ig_out_$PBS_JOBID

# Clean up
rm -f machines

echo "CharLES Launch in $(hostname) FINISHED!!!"
