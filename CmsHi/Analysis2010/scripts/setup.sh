#!/bin/sh

indir=`pwd` 
cd $CMSSW_BASE/src/


cvs co CondFormats/HIObjects
cvs co DataFormats/HeavyIonEvent
cvs co RecoHI/HiCentralityAlgos
#cvs co -d CmsHi/Analysis2010 UserCode/CmsHi/Analysis2010

# Optional for other pings than high-pt:
cvs co -d CmsHi/JetAnalysis UserCode/CmsHi/JetAnalysis

scram b

cd $indir


