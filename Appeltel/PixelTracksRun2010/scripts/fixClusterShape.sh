#!/bin/sh

cd $CMSSW_BASE/src



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


scram b
