#!/bin/sh

cd $CMSSW_BASE/src

# Centrality Fix
cvs co CondFormats/HIObjects
cvs co DataFormats/HeavyIonEvent
cvs co RecoHI/HiCentralityAlgos
cvs co -d CmsHi/Analysis2010 UserCode/CmsHi/Analysis2010

# Event Selection
cvs co -r V01-02-09 HeavyIonsAnalysis/Configuration

# Full Track Selection Packages
cvs co -d edwenger/HiTrkEffAnalyzer UserCode/edwenger/HiTrkEffAnalyzer
cvs co -d edwenger/HiVertexAnalyzer UserCode/edwenger/HiVertexAnalyzer 

scram b
