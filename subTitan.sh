#!/bin/bash

#Read in source code path, then shift for optargs
srcCode="$1"
shift

date=`date +%Y-%m-%d`
parentPath="Module_$srcCode"
jobname=Analysis_${srcCode%/}_$date
dirPath="$HOME/scratch/$parentPath"
user=${USER}
jobid="JA";
cores=1
walltime=12:00:00
memory=12g
outfile="Jobname.o"
nJobs=10
nMC=100
nPop=100000
seed=0
simT=120
burn=36
repeats=1
model=${PWD##*/}
basePath=$PWD

usage() {
echo "
usage: subtitan {SourceFolder} [-r repeats] [-n iterations] [-T walltime] [-m memory] [-s seed]
                [-o outfile] [-N population] [-t timerange] [-b burntime] [-j jobname]

Starts a TITAN simulation in ~/scratch/{SourceFolder}/{jobname}

options:
  -j jobname	  name of analysis for organization (default: {SourceFolder}_date)
  -r repeats      number of times to repeat the analysis (default: $repeats)
  -n iterations   number of mode iterations per job (default: $nMC)
  -N population	  number of agents in population (default: $nPop)
  -T walltime     as hh:mm:ss, max compute time (default: $walltime)
  -m memory       as #[k|m|g] (default: $memory)
  -o outfile      save a copy of the session's output to outfile (default: off)
  -s seed         random seed for model [0 random, -1 stepwise] (default: $seed)
  -t timerange	  number of time steps per iteration in (default: $simT)
  -b burntime	  number of time steps to burn for equilibration (default: $burn)     
"
echo $dirPath
exit 0
}

updateParams() {
echo "

    Updating params:
	curPath		$PWD
	sourceCode	$srcCode
	jobname: 	$jobname
	iterations: 	$nMC
	population: 	$nPop
	seed:		$seed
	time:		$simT
	burn:		$burn
	
	walltime	$walltime
	memory		$memory
"

#TITAN params
sed -i "10s/"0"/$seed/g" params.py
sed -i "11s/"100"/$nMC/g" params.py
sed -i "12s/"100000"/$nPop/g" params.py
sed -i "13s/"120"/$simT/g" params.py
sed -i "14s/"36"/$burn/g" params.py

#Submit script params
sed -i "s/MODEL_NAME/$jobname/g" bs_Core.sh
sed -i "s/WALL_TIME/$walltime/g" bs_Core.sh

}

prepSubmit() {

    #Copy source code into parent path
    echo -e "\n\tCopying $srcCode to $finalPath"
    cp -rT $srcCode $finalPath

    #Move into new source code folder
    echo -e "\n\tMoving to model folder directory"
    cd $finalPath
    echo -e "\t$PWD"
    updateParams;

    #Submit job to cluster
    #sbatch bs_Core.sh

    #Move back to base directory
    cd $basePath
}

while getopts j:n:N:t:m:s:c:T:b:r: option
do
    case "${option}"
        in
	N) nPop=${OPTARG};;
	m) mem=${OPTARG};;
    n) nMC=${OPTARG};; 
    s) seed=${OPTARG};;
    j) jobname=${OPTARG};;
	t) simT=${OPTARG};;
	b) burn=${OPTARG};;
	T) walltime=${OPTARG};;
	r) repeats=${OPTARG};;
    esac
done


# User and Date will be ignored if job ID is specified

if [ ! $srcCode ]; then
    usage;
fi

if [ ! -d $srcCode ]; then
    echo -e "\n\n$srcCode is not a directory! Source code must be provided as a directory\n"
    exit 0
fi

if [ -d $dirPath$jobname ]; then
    echo -e "\n\n!! WARNING !!\nThe folder $jobname already exists and will be OVERWRITTEN!\n"
    read -p "Continue (y/n)?" choice
    case "$choice" in 
      y|Y ) echo "Proceeding";;
      n|N|* ) echo "Aborting"
	    exit 0;;
    esac
fi

if [ $srcCode ]; then

    echo "
        jobname     $jobname
        dirpath	    $dirPath
        user	    $user
        date        $date
        cores       $cores
        walltime    $walltime
        memory      $memory
        outfile     $outfile
        nMc         $nMC
        seed        $seed
        model       $model"
    echo -e "\n" 
    
    echo -e "\tMaking parent directory in scratch"
    mkdir -p $dirPath
    echo -e "\t $dirPath"

    if [ $repeats -gt 1 ]; then
        mkdir -p $dirPath$jobname
        basejobname=$jobname
        for ((i=1; i<=repeats; i++)); do
            echo -e "\n\nWorking on repeat $i"
            jobname=$basejobname"_"$i
            finalPath=$dirPath$basejobname"/"$jobname
            prepSubmit;
        done
    else
        finalPath=$dirPath$jobname
        prepSubmit;
    fi
    
else
    echo -e "\nSOMETHING WENT WRONG!!! Abort"
    exit 1;
fi

