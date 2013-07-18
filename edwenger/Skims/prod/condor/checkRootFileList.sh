inputs=$1
output=log/corruptRootFiles.txt
printf "" > $output
for i in `ls $inputs`; do
  echo $i;
  $CMSSW_BASE/src/edwenger/Skims/prod/condor/checkRootFile.sh $i
  if [ $? -ne 0 ]; then echo $i >> $output; fi
done
