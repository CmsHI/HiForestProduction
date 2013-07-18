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

# Fix for pixel tracking cluster shape crash

cvs co -r $CMSSW_VERSION RecoPixelVertexing/PixelLowPtUtilities
cvs co -r $CMSSW_VERSION RecoPixelVertexing/PixelTrackFitting
cvs co -r $CMSSW_VERSION RecoHI/HiTracking

cp Appeltel/PixelTracksRun2010/scripts/PixelTrackReconstruction.ccfix \
   RecoPixelVertexing/PixelTrackFitting/src/PixelTrackReconstruction.cc

cp Appeltel/PixelTracksRun2010/scripts/PixelTrackFilter.hfix \
   RecoPixelVertexing/PixelTrackFitting/interface/PixelTrackFilter.h

cp Appeltel/PixelTracksRun2010/scripts/ClusterShapeTrackFilter.hfix \
   RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeTrackFilter.h

cp Appeltel/PixelTracksRun2010/scripts/ClusterShapeTrackFilter.ccfix \
   RecoPixelVertexing/PixelLowPtUtilities/src/ClusterShapeTrackFilter.cc

cp Appeltel/PixelTracksRun2010/scripts/HIPixelTrackFilter.ccfix \
   RecoHI/HiTracking/src/HIPixelTrackFilter.cc

cp Appeltel/PixelTracksRun2010/scripts/HIPixelTrackFilter.hfix \
   RecoHI/HiTracking/interface/HIPixelTrackFilter.h

# Calorimeter Cleaning

addpkg RecoLocalCalo/HcalRecAlgos
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/src/HBHETimingShapedFlag.cc
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/interface/HBHETimingShapedFlag.h
cvs co -r V00-00-18 JetMETAnalysis/HcalReflagging
cvs co -d HcalFilter UserCode/MTonjes/HcalFilter
#ecal requires some code from the photon analyzer, 'rm' commands are to truncate the endless chain of packages
# obviously, this should be cleaned up
cvs co -d CmsHi/PhotonAnalysis UserCode/CmsHi/PhotonAnalysis
rm CmsHi/PhotonAnalysis/plugins/MultiPhotonAnalyzer.*
rm CmsHi/PhotonAnalysis/plugins/SinglePhotonAnalyzer.*
cvs co -d UserCode/HafHistogram UserCode/CmsHi/HafHistogram

# dE/dx estimators for pixel tracks

cvs co -r V00-04-28 RecoTracker/DeDx

scram b
