#!/bin/bash

while true; do
  for jobname in `tail -5 completed_runs_names.txt`; do
    echo $jobname
    mkdir -p /tmp/velicanu/$jobname

    crab -c $jobname -status > teststatus.txt

    nstaged=`ls -l /tmp/velicanu/$jobname/openhlt_data_*.root | wc -l`
    ntot=`grep "Total Jobs" teststatus.txt | awk '{print $2}'`
    
    echo $nstaged $ntot
    if [ $nstaged -lt $ntot ]; then
       ./getEosDir.sh /store/caf/user/velicanu/HIExpressPhysics/$jobname /tmp/velicanu/$jobname
       hadd -f /tmp/velicanu/$jobname/$jobname.root /tmp/velicanu/$jobname/openhlt_data_*.root
       rfcp /tmp/velicanu/$jobname/$jobname.root /castor/cern.ch/user/v/velicanu/openHLT/
      fi
  done
  # sleep 1
  sleep 900
done


    # if [ $ndone -eq $ntot ]; then
    # fi