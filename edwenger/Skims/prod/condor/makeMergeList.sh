#!/bin/bash -
#inDir=/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/MinimumBias/MB-C10-PR9-TRKANASKIM-v3/ae98f896d123ace1f592d26e790fa90c
inDir=$1
tag=trkhists
curDir=`pwd`

printf "hadd -T all/${tag}_histOnly_all.root " > $curDir/log/mergeList_${tag}.txt
ls -lh $inDir/${tag}_* | awk '{print $NF}' | awk '{printf "dcache:%s ", $0}' >> $curDir/log/mergeList_${tag}.txt

