#!/bin/sh

cd $CMSSW_BASE/src

# skiming code
cvs co UserCode/edwenger/Skims
cvs co -r 1.4 UserCode/edwenger/Skims/python/patAna_cff.py

# analyzers
cvs co UserCode/edwenger/EvtSelAnalyzer
cvs co UserCode/edwenger/VertexAnalyzer
cvs co UserCode/edwenger/TrackSpectraAnalyzer
cvs co UserCode/edwenger/TrkEffAnalyzer
#cvs co -d Saved UserCode/SavedFMa/Saved

# agglomerative vertex
cvs co -r V00-00-07 UserCode/FerencSiklerVertexing

# vertex constraint for refitting
cvs co UserCode/ASYoon/PPTrackingTools/VertexConstraintProducer

# cleanup
mv UserCode/edwenger .
mv UserCode/ASYoon/PPTrackingTools/ .
rm -r UserCode/ASYoon

scram b
