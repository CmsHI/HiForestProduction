result=0
root -b -q $1 >& /tmp/openrootfile.txt
cat /tmp/openrootfile.txt | grep Error
if [ $? -eq 0 ]; then result=1; fi
exit $result
