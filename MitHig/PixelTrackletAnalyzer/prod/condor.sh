#!/bin/bash

#data=pion
data=pp900GeV
#cfg=count_cfg.py
#cfg=count_from_reco_cfg.py
cfg=$2
tag=1001
dir=`pwd`
ldir=$dir/joblogs
outfile=$2-output-`date +"%s"`

subfile=subfile_$outfile

echo $subfile
mkdir -p $ldir
mkdir -p /tmp/yetkin

cat > $subfile <<EOF

Universe     = vanilla

Notification = Error
Executable   = $dir/$1
Arguments    = $2 $3 $4 $5 $6 $7 $8 $9
Requirements = (Mips > 900) && (ARCH=="X86_64")
Rank         = Mips
GetEnv       = True

Initialdir   = $dir
Input        = /dev/null
Output       = /tmp/yetkin/$outfile.out
Error        = /tmp/yetkin/$outfile.err
Log          = /tmp/yetkin/$outfile.log

transfer_input_files = $dir/$cfg
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT

Queue
EOF

sleep 0
cat $subfile

echo Executable   = $dir/$2
echo Arguments    = $3 $4 $5 $6 $7 $8 $9


condor_submit $subfile
sleep 0
rm $subfile
