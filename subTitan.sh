#!/bin/bash

#Read in source code path, then shift for optargs
version="0.1.0"
titanPath="/gpfs/data/bm8/TITAN/TITAN/"
settingPath="$1"
shift

if [ $settingPath ]; then
    setting=$(basename ${settingPath%.py})
fi

date=`date +%Y-%m-%d-T%H-%M-%S`
srcCode="${titanPath}titan/"
parentPath="Module_$setting/"
jobname=Analysis_$setting_$date
outPath="$HOME/scratch/$parentPath"
user=${USER}
jobid="JA";
cores=1
walltime=12:00:00
memory=12g
outfile="Jobname.o"
nJobs=1
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
usage: subtitan {Parameter file} [-r repeats] [-n iterations] [-T walltime] [-m memory] [-s seed]
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
echo "TITAN ver: "$version
exit 0
}

updateParams() {
echo "

    Updating params:
	savePath	$PWD
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
# sed -i "s/\(rSeed = \)\([0-9]*\)/\1${seed}/g" titan/params.py
sed -i "s/\(N_MC = \)\([0-9]*\)/\1${nMC}/g" titan/params.py
sed -i "s/\(N_POP = \)\([0-9]*\)/\1${nPop}/g" titan/params.py
sed -i "s/\(TIME_RANGE = \)\([0-9]*\)/\1${simT}/g" titan/params.py
sed -i "s/\(burnDuration = \)\([0-9]*\)/\1${burn}/g" titan/params.py

#Submit script params
sed -i "s/MODEL_NAME/$jobname/g" scripts/bs_Core.sh
sed -i "s/WALL_TIME/$walltime/g" scripts/bs_Core.sh
sed -i "s/MEMORY/$memory/g" scripts/bs_Core.sh

}

prepSubmit() {

    #Copy source code into parent path
    #echo -e "\n\tMoving setting $setting into $srcCode"
    #cp $settingPath $srcCode/params.py
    echo -e "\n\tCopying $srcCode to $finalPath"
    mkdir -p $finalPath
    cp $titanPath/run_titan.py $finalPath
    cp -rT $titanPath/titan $finalPath/titan
    cp -rT $titanPath/scripts $finalPath/scripts
    mkdir -p $finalPath/results/network
    cp $settingPath $finalPath/titan/params.py
    #Move into new source code folder
    echo -e "\n\tMoving to model folder directory"
    cd $finalPath
    echo -e "\t$PWD"
    updateParams;

    #Submit job to cluster
    sbatch scripts/bs_Core.sh

    #Move back to base directory
    cd $basePath
}

while getopts j:n:N:t:m:s:c:T:b:r: option
do
    case "${option}"
        in
	N) nPop=${OPTARG};;
	m) memory=${OPTARG};;
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

if [ ! $settingPath ]; then
    usage;
fi

if [ ! -d $srcCode ]; then
    echo -e "\n\n$srcCode is not a directory! Source code must be provided as a directory\n"
    exit 0
fi

if [ -d $outPath$jobname ]; then
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
        outPath	    $outPath
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
    mkdir -p $outPath
    echo -e "\t $outPath"

    if [ $repeats -gt 1 ]; then
        mkdir -p $outPath$jobname
        basejobname=$jobname
        for ((i=1; i<=repeats; i++)); do
            echo -e "\n\nWorking on repeat $i"
            jobname=$basejobname"_"$i
            finalPath=$outPath$basejobname"/"$jobname
            prepSubmit;
        done
    else
        finalPath=$outPath$jobname
        prepSubmit;
    fi

else
    echo -e "\nSOMETHING WENT WRONG!!! Abort"
    exit 1;
fi
