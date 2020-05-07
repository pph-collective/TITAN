#!/bin/bash
#SBATCH -J MODEL_NAME
#SBATCH -c NCORES
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

#!/bin/bash
module load graphviz

setting=""
paramPath=""
nMC=""
useBase=""
forceFlag=""
sweepDefs=""
sweepfile=""
rows=""

while getopts S:n:b:w:f:p: option
do
    case "${option}"
        in
	S) setting=${OPTARG};;
	n) nMC=${OPTARG};;
	b) useBase=${OPTARG};;
	w) sweepDefs+="-w ${OPTARG} ";;
	W) sweepfile="-W ${OPTARG}";;
	r) rows="-r ${OPTARG}";;
	F) forceFlag="-F";;
	p) paramPath=${OPTARG};;
    esac
done

cd $PWD

echo Master process running on `hostname`
echo Directory is `pwd`
echo PBS has allocated the following nodes:
echo `cat $PBS_NODEFILE`
echo Starting execution at `date`
NCPU=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NCPU CPUs

./run_titan.py -S $setting -n $nMC -p $paramPath -b $useBase $forceFlag $sweepDefs $sweepfile $rows
