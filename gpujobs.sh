#!/bin/bash
#SBATCH -N 1####nodes
#SBATCH -J aagaa_drude_remd
#SBATCH -p gpu
#SBATCH -N 4####procode
#SBATCH -n 8####procode
#SBATCH -w node[18,19,20]
#SBATCH -c 4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=2
#SBATCH -o /u/shiyi/4huang_remd/test.out
#SBATCH -e /u/shiyi/4huang_remd/test.err


module purge
module load mathlib/cuda/9.1.85_387.26

module load intel/composer_xe_2019.0.117
module load intel/parallel_studio_xe_2019.0.117


export PATH=/home/shiyi/mympich/bin:$PATH:/home/shiyi/anaconda3/bin/
export LD_LIBRARY_PATH=/home/shiyi/mympich/lib:$LD_LIBRARY_PATH
#source /home/shiyi/amber18/amber.sh

export PREFIX="test5"
mpiexec -hostfile test5.hostfile -configfile test5.configfile ./run_remd.py

## request the resource, but not use it (for debug)
#sleep infinity

#build_mpirun_configfile --hostfilepath $PREFIX.hostfile --configfilepath $PREFIX.configfile "./test_book.py"
#build_mpirun_configfile --hostfilepath $PREFIX.hostfile --configfilepath $PREFIX.configfile "./my_remd2.py"
#build_mpirun_configfile --hostfilepath $PREFIX.hostfile --configfilepath $PREFIX.configfile "./run2021_1_syremd3.py"
#build_mpirun_configfile --hostfilepath $PREFIX.hostfile --configfilepath $PREFIX.configfile "./run2021_1_syremd3.py"
#build_mpirun_configfile --hostfilepath $PREFIX.hostfile --configfilepath $PREFIX.configfile "./run2021_1_syremd3_june.py"
#
#
