#!/bin/sh                                                                                                                                                                    

cd $CMSSW_BASE/src

# trak analyzer for spectra
cvs co UserCode/ASYoon/SpectraAna/TrackSpectraAnalyzer
mv UserCode/ASYoon/SpectraAna/ .

# MTV hack-up to include 2D histogams 
cvs co -r CMSSW_3_5_6 Validation/RecoTrack/
cvs co UserCode/ASYoon/PPTrackingTools/MTVPlugIn
cp UserCode/ASYoon/PPTrackingTools/MTVPlugIn/interface/MultiTrackValidatorBase.h Validation/RecoTrack/interface
cp UserCode/ASYoon/PPTrackingTools/MTVPlugIn/plugins/MultiTrackValidator.cc Validation/RecoTrack/plugins

# to run MTV with HitPixelLayersTPSelector.cc 
cvs co UserCode/CmsHi/TrackAnalysis
mv UserCode/CmsHi/ .
