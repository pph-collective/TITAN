#!/bin/bash
#SBATCH -J MODEL_NAME
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=WALL_TIME
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=NoMail
#SBATCH --mem=MEMORYGB

if [ -z "$SLURM_NPROCS" ] ; then
  if [ -z "$SLURM_NTASKS_PER_NODE" ] ; then
    SLURM_NTASKS_PER_NODE=1
  fi
  SLURM_NPROCS=$(( $SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE ))
fi
#Conduct Dual analysis for 100 runs

#!/bin/bash

cd $PWD

echo Master process running on `hostname`
echo Directory is `pwd`
echo PBS has allocated the following nodes:
echo `cat $PBS_NODEFILE`
echo Starting execution at `date`
NCPU=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NCPU CPUs

# execute an MPI program
# $SLURM_NPROCS = nodes x ppn
# Change global N_MC in MPI_simulation.py to $SLURM_NPROCS
# and PROCESSES to a multiple of $SLURM_NPROCS for optimal distribution
python3 run_titan.py
