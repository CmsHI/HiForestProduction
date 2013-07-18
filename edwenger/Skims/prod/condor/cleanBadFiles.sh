badList=$1
for i in `cat $badList`; do ls | grep $i | sed 's/^/rm /'; done
