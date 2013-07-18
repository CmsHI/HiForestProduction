#~/bin/bash -
if [ $# -lt 3 ]; then
  echo Usage:
  echo "  $0 <inDir> <tag> <outDir>"
  echo "  * inDir is a dir containing the root files you want to merge"
  exit 1
fi

script=mergeFiles.sh
inDir=$1
tag=$2
outdir=$3

# initialize
ct=0
cmd=
nPerMerge=500
N=`ls $inDir/*.root | grep $tag | wc -l`
startDir=`pwd`
echo $inDir: $N files

# clean indir
./cleanJob2.sh $inDir

for job in `ls $inDir/*.root | grep $tag | head -n $N`; do
  subi=`expr $ct % $nPerMerge`
  #echo $job, subi: $subi
  # -- build subgroup --
  if [ $subi -eq 0 ]; then
    cmd="./condor64Basic.sh $script $outdir trkhists_sub$ct.root $job"
    #cmd="$./script $outdir trkhists_sub$ct.root $job"
  else
    #echo $ct: $job
    cmd="$cmd $job"
  fi
  # -- merge subgroup --
  if [ $subi -eq $((nPerMerge-1)) -o $ct -eq $((N-1)) ]; then
    echo $cmd
    eval $cmd
  fi
  ct=$((ct+1))
done
