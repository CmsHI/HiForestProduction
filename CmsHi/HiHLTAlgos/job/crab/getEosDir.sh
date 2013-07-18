if [ $# -ge 2 ]; then
  indir=$1
  outdir=$2
  for i in `cmsLs $indir | grep root | awk '{print $NF}'`; do
    if [ -f $outdir/$i ]; then echo $i already exits will not copy
    else cmsStage $i $outdir
    fi
  done
else
  echo "Usage: $0 <inputdir> <ouputdir>"
fi
