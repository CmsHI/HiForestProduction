#!/bin/bash -
if [ $# -eq 0 ]; then
   echo Usage
   echo "   $0 <cfg.py> <out_dir> <outputfile> <maxEvents> <input_file1,input_file2,...>"
   exit 1
fi

# input parameters
cfg=$1
outputDir=$2
outputFile=$3
maxEvents=$4
inputFiles=$5
swenv=/net/hisrv0001/home/frankma/scratch01/work/job/HLT/hi2011hlt_440/env.sh

# check
if [ ! -e $cfg ]; then
   echo cfg: $cfg does not exist!
   exit 1
fi

# setup job parameters
echo $outputFile
outputFileName=${outputFile##*/}
echo $outputFileName
logDir=$outputDir/log
logFile=${outputFileName%.root}.txt
echo $logFile


# setup env
. /net/hisrv0001/home/frankma/UserCode/SavedFMa/bin/cgate/basic/osg_cmssw_set_basic2.sh
. $swenv

mkdir -p $logDir
cmd="cmsRun $cfg maxEvents=$maxEvents files=$inputFiles output=$outputFile >& $logFile"

# review
which cmsRun
echo $cmd

# run!
eval $cmd
mv $logFile $logDir/
fsize=`stat -c%s $outputFile`
echo file size: $fsize
if [ $fsize -lt 50 ]; then
  echo bad file: $outputFile
  exit 1
fi

grep 'Event  Summary' $logDir/$logFile
if [ $? -eq 0 ]; then
  mv $outputFile $outputDir/
fi
