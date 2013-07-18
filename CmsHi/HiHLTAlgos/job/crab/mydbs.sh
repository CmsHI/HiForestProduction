#!/bin/bash


while ( true )
do

  dbs search --query='find run where dataset like /HIExpressPhysics/HIRun2011-Express-v1/FEVT' > currentdbs.txt

  cat currentdbs.txt | grep -v DBS | grep -v \- | grep -v run | sort | uniq > tmp.txt
  cat tmp.txt > currentdbs.txt
  cat completed_runs.txt | grep -v DBS | grep -v \- | grep -v run | sort | uniq > tmp.txt
  cat tmp.txt > completed_runs.txt
  
  diff currentdbs.txt completed_runs.txt > diff.txt

  cat diff.txt
  export lastrun=`cat lastrun.txt`
  

  for line in `cat diff.txt`
  do 
    if [[ $line != *,* && $line -gt 100000 && $line -lt 999999 ]] # valid run that is not yet completed
    then
      if [[ $line == $lastrun ]]
      then
        echo is last run $line , extend and submit
        multicrab -extend
        multicrab -submit
      else
        echo not last run $line , create and submit
        sed s@_run_selection_flag_@$line@ autohlt.cfg > multicrab.cfg
        multicrab -create
        multicrab -submit
        
        echo $lastrun >> completed_runs.txt
        export lastrun=$line
        echo $lastrun > lastrun.txt
        
        grep -o hiexp*$line*autohlt multicrab.cfg >> completed_runs_names.txt
        tail -5 completed_runs_names.txt > tmp.txt
        mv tmp.txt completed_runs_names.txt
      fi
    else
      echo not valid run $line
    fi
  done 
  sleep 1800
done

      

      # fi



