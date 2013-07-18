#!/bin/bash -
if [ $# -eq 0 ]; then
   echo Usage
   echo "   $0 <input_dir> <pref:input_file.root> <output_dir> <nperjob>"
   exit 1
fi

# input parameters
cfg=custom_TrkAnaSkim_data_7TeV_cfg.py
inputDir=$1
inputFile=$2
outputDir=$3
nperjob=$4
swenv=/net/hisrv0001/home/frankma/work/job/ana/trkAnaCMSSW_3_6_2/env.sh

# check
if [ ! -e $cfg ]; then
   echo cfg: $cfg does not exist!
   exit 1
fi

# setup job parameters
inputFileName=${inputFile##*/}
echo $inputFileName
outputFileName=trkhists_$inputFileName
echo $outputFileName
outputFile=$outputDir/$outputFileName
echo $outputFile
logDir=$outputDir/log
echo $logDir
logFile=$logDir/${outputFileName%.root}.txt


# setup env
. /net/hisrv0001/home/frankma/UserCode/SavedFMa/bin/cgate/basic/osg_cmssw_set_basic.sh
. $swenv

mkdir -p $logDir
cmd="cmsRun $cfg maxEvents=$nperjob files=dcache:$inputDir/$inputFile output=$outputFile >& $logFile"
#for itry in {1..5}; do
#  if [ -e /tmp/$inputFile ]; then
#    echo Copied $inputDir/$inputFile
#    break
#  else
#    sleep 10
#    echo "dccp $inputDir/$inputFile /tmp/"
#    dccp $inputDir/$inputFile /tmp/
#  fi
#done
#cmd="cmsRun $cfg maxEvents=$nperjob files=file:/tmp/$inputFile output=$outputFile >& $logFile"

# review
which cmsRun
echo $cmd

# run!
eval $cmd
#rm /tmp/$inputFile
