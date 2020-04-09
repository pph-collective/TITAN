#!/bin/bash

#Read in source code path, then shift for optargs
titanPath=$PWD
paramPath="$1"
shift

setting="custom"
date=`date +%Y-%m-%d-T%H-%M-%S`
user=${USER}
walltime=12:00:00
memory=12g
outfile="Jobname.o"
repeats=1
nMC=1
model=${PWD##*/}
basePath=$PWD
useBase="True"
jobname=""
folderName=""
sweepDefs=""
force=false

while getopts m:S:T:j:r:n:f:w:F:t: option
do
    case "${option}"
        in
	m) memory=${OPTARG};;
  j) jobname=${OPTARG};;
	T) walltime=${OPTARG};;
	S) setting=${OPTARG};;
	r) repeats=${OPTARG};;
	n) nMC=${OPTARG};;
	b) useBase=${OPTARG};;
	f) folderName=${OPTARG};;
	w) sweepDefs+="-w ${OPTARG} ";;
	F) force=true;;
	f) titanPath=${OPTARG};;
    esac
done

if [[ $jobname == "" ]]; then
	jobname="Analysis_$setting_$date"
fi

if [[ $folderName == "" ]]; then
	folderName="$setting/"
fi

srcCode="${titanPath}titan/"
outPath="$HOME/scratch/$folderName"

usage() {
echo "
usage: subtitan {Parameter file or directory}[-T walltime] [-m memory] [-S setting] [-j jobname] [-r repeats] [-n iterations] [-b use_base] [-f folder_name] [-w sweep_defs] [-F force] [-t titanPath ]

Starts a TITAN simulation in ~/scratch/{SourceFolder}/{jobname}

options:
  -j jobname	  name of analysis for organization (default: {SourceFolder}_date)
  -T walltime     as hh:mm:ss, max compute time (default: $walltime)
  -m memory       as #[k|m|g] (default: $memory)
	-S setting      name of setting for this model
	-r repeats      number of times to repeat the analysis (default: $repeats)
  -n iterations   number of mode iterations per job (default: $nMC)
	-b use_base     whether to use the base setting as True or False (default: $useBase)
	-f folder_name	What the parent folder for the model run outputs should be called (default: <setting>)
	-w sweep_defs   Optionally, definitions of sweep parameters in the format param:start:stop[:step]
	-F force				If the number of sweep combinations exceeds 100, run anyway
	-t titanPath		where the code is
"
exit 0
}

updateParams() {
echo "

    Updating params:
	savePath	$PWD
	sourceCode	$srcCode
	jobname: 	$jobname
	walltime	$walltime
	memory		$memory
"

#Submit script params
sed -i "s/MODEL_NAME/$jobname/g" scripts/bs_Core.sh
sed -i "s/WALL_TIME/$walltime/g" scripts/bs_Core.sh
sed -i "s/MEMORY/$memory/g" scripts/bs_Core.sh

}

prepSubmit() {

    #Copy source code into parent path
    echo -e "\n\tCopying $srcCode to $finalPath"
    mkdir -p $finalPath
    cp $titanPath/run_titan.py $finalPath
    cp -rT $titanPath/titan $finalPath/titan
    cp -rT $titanPath/scripts $finalPath/scripts
	  cp -rT $titanPath/settings $finalPath/settings
    #Move into new source code folder
    echo -e "\n\tMoving to model folder directory"
    cd $finalPath
    echo -e "\t$PWD"
    updateParams;

		# set up sweeping flags
		forceFlag=""
		if [ $force = true ]; then
			forceFlag=" -F"
		fi


    #Submit job to cluster
    sbatch scripts/bs_Core.sh -S $setting -p $paramPath -n $nMC -b $useBase $forceFlag $sweepDefs

    #Move back to base directory
    cd $basePath
}


# User and Date will be ignored if job ID is specified

if [ ! $paramPath ]; then
    usage;
fi

if [[ ${paramPath:0:1} != "/" ]] || [[ ${paramPath:0:1} == "~" ]]; then
	paramPath="${pwd}/$paramPath"
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
				paramPath   $paramPath
        user	    $user
        date        $date
        walltime    $walltime
        memory      $memory
        outfile     $outfile
        model       $model"
    echo -e "\n"

    echo -e "\tMaking parent directory in scratch"
    mkdir -p $outPath
    echo -e "\t $outPath"

    if [ $repeats -gt 1 ]; then
        # mkdir -p $outPath$jobname
        basejobname=$jobname
        for ((i=1; i<=repeats; i++)); do
            echo -e "\n\nWorking on repeat $i"
            jobname=$basejobname"_"$i
            finalPath=$outPath"/"$basejobname"/"$jobname
            prepSubmit;
        done
    else
        finalPath=$outPath"/"$jobname
        prepSubmit;
    fi

else
    echo -e "\nSOMETHING WENT WRONG!!! Abort"
    exit 1;
fi
