#~/bin/bash -
if [ $# -lt 3 ]; then
  echo Usage:
  echo "  $0 <cfg> <input_list> <output_dir>"
  exit 1
fi

cfg=$1
inputList=$2
output_dir=$3

# initialize
ct=0
cmd=
nPerMerge=10
N=`cat $inputList | wc -l`
startDir=`pwd`
echo $inputList: $N files

for job in `cat $inputList | head -n 3000`; do
  subi=`expr $ct % $nPerMerge`
  #echo $job, subi: $subi
  # -- build subgroup --
  if [ $subi -eq 0 ]; then
    cmd="./condor64Data.sh runMany.sh $cfg $output_dir tr_trigana_$ct.root -1 $job"
    echo $ct
    sleep 1
  else
    cmd="$cmd,$job"
  fi
  # -- merge subgroup --
  if [ $subi -eq $((nPerMerge-1)) -o $ct -eq $((N-1)) ]; then
    echo $cmd
    eval $cmd
  fi
  ct=$((ct+1))
done
