logdir=$1
cd $logdir
grep -l 'FileOpenError' trkhists_trkAnaSkimAOD_*.txt | sed 's/txt/root/' > badFileList.txt
cd -
echo $logdir/badFileList.txt
