#!/bin/bash

dir=`pwd`
ldir=~/condor
outfile=$1-output-`date +"%s"`

subfile=subfile_$outfile

echo $subfile
mkdir -p $ldir
mkdir -p /tmp/frankma
exe=$dir/$1
shift

cat > $subfile <<EOF

Universe     = vanilla

Notification = Error
Executable   = $exe
Arguments    = $@
Requirements = (Mips > 900) && (ARCH=="X86_64")
Image_Size = 2400 Meg
Rank         = Mips
GetEnv       = True

Initialdir   = $dir
Input        = /dev/null
Output       = $ldir/$outfile.out
Error        = $ldir/$outfile.err
Log          = /tmp/frankma/$outfile.log
Notify_user  = frankma@mit.edu

+AccountingGroup = "group_cmshi.frankma"
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT

Queue
EOF

sleep 1
cat $subfile

echo Executable   = $exe
echo Arguments    = $@


condor_submit $subfile
sleep 1
rm $subfile
