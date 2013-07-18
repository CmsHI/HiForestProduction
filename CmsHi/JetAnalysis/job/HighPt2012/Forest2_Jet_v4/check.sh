#!/bin/bash
for (( ; ; ))
do
  for job in \
   promptskim-hihighpt-hltjet80-pt90-v4-lumi1 \
   promptskim-hihighpt-hltjet80-pt90-v4-lumi2 \
   promptskim-hihighpt-hltjet80-pt90-v4-lumi3 \
   promptskim-hihighpt-hltjet80-pt90-v4-lumi4 \
   promptskim-hihighpt-hltjet80-pt90-v4-lumi5 \
   promptskim-hihighpt-hltjet80-pt90-v4-lumi6
   do
    crab -c $job -status>& crablog.tmp
    cat crablog.tmp
    cat crablog.tmp |grep 'Jobs Aborted'|awk '{print "./beep"}'|bash
    cat crablog.tmp |awk '{where = match($0,"Jobs Abort")}{if (where) start=1}{gotit = match($0,"List of jobs:")}{if (start&&gotit) {print "crab -c '$job' -resubmit "$4;start=false}}'|bash
    cat crablog.tmp |awk '{where = match($0,"Jobs with Wrap")}{Zero = match($0,"Jobs with Wrapper Exit Code : 0")}{if (Zero) where=false}{if (where) start=1}{gotit = match($0,"List of jobs:")}{if (start&&gotit) {print "crab -c '$job' -resubmit "$4;start=false}}'|bash
    cat crablog.tmp |awk '{where = match($0,"Jobs Done")}{if (where) start=1}{gotit = match($0,"List of jobs:")}{if (start&&gotit) {print "crab -c '$job' -get "$4;start=false}}'|bash
    cat crablog.tmp |awk '{where = match($0,"Jobs Done")}{if (where) print "crab -c '$job' -getoutput"}'|bash   
    cat crablog.tmp |awk '{where = match($0,"Jobs Created")}{if (where) print "crab -c '$job' -submit 200"}'|bash
  done

  echo "Waiting for the next chance..."
  sleep 999
done
