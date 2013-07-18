#!/bin/bash -
if [ $# -eq 0 ]; then
  echo Usage
  echo "  $0 <indir or inputlist> [fromList=0] [tag=trkAnaSkimAOD]"
  exit 1
fi
#inDir=/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/MinimumBias/MB-C10-PR9-TRKANASKIM-v2/38a0275e6a2a62980cf9b1868edbd493
inDir=$1
fromList=0
tag=trkAnaSkimAOD
if [ $# -ge 2 ]; then
  fromList=$2
elif [ $# -ge 3 ]; then
  tag=$3
fi

if [ $fromList -eq 1 ]; then
  inputs=$inDir
else
  inputs=log/dirinputs.txt
  find $inDir | grep $tag | sort > $inputs
fi

for file in `cat $inputs` ; do
  base=`echo $file | sed 's/_[0-9]\.root$//'`
  #echo base: $base
  for i in `seq 1 5`; do
    this=${base}_${i}.root
    grep -q $this $inputs
    thisExists=$?
    next=${base}_$((i+1)).root
    grep -q $next $inputs
    nextExists=$?
    if [ $((thisExists+nextExists)) -eq 0 ]; then
      #echo next: $next
      echo $this
    fi
  done
done
