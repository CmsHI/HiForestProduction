#!/bin/bash -
if [ $# -lt 3 ]; then
  echo "Usage:"
  echo "$0 <outdir> <outfile> <file1.root file2.root ...>"
  exit 1
fi

swenv=/net/hisrv0001/home/frankma/work/job/ana/trkAnaCMSSW_3_6_2/env.sh

outdir=$1
shift
outfile=$1
shift
infiles=$@

echo "outdir: $outdir"
echo "outfile: $outfile"
#echo "infiles: $infiles"

# setup env
. /net/hisrv0001/home/frankma/UserCode/SavedFMa/bin/cgate/basic/osg_cmssw_set_basic.sh
. $swenv
#which hadd

mkdir -p $outdir/log
logfile=${outfile%.root}.txt
mergecmd="hadd -T $outdir/$outfile $infiles"
#mergecmd="echo hadd -T $outdir/$outfile $infiles"
runcmd="(time $mergecmd >& $outdir/log/$logfile) 2>> $outdir/log/$logfile"
#echo $runcmd
eval $runcmd
