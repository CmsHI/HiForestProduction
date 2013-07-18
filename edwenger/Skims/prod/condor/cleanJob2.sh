indir=$1
outdir=$indir/bad

mkdir -p $outdir
for i in `find $indir -maxdepth 1 -type f -size -10 | grep root$`; do mv $i $outdir; done
