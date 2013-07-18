#!/bin/bash -
if [ $# -lt 2 ]; then
  echo Usage:
  echo "  $0 <input_list> <output_dir>"
  exit 1
fi

#inputTopDir=$1
inputList=$1
output_dir=$2

tag=trkAnaSkimAOD
#prefix=dcache:
#prefix=file:
njobs=1000000

echo begin
#for i in `ls $inputTopDir | grep $tag | sort | head -n $njobs`
for i in `cat $inputList | sort | head -n $njobs`
do
  ifile=${i##*/}
  inputTopDir=/pnfs/cmsaf.mit.edu/t2bat/cms${i%/*.root}
  cmd="./condor64.sh runOne.sh $inputTopDir $ifile $output_dir -1"
  echo $cmd
  eval $cmd
done

