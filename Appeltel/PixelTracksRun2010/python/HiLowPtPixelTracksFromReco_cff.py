import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtPixelTracks_cfi import *
from Appeltel.PixelTracksRun2010.HIGoodPixelTracks_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtWideCutPixelTracks_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtConformalPixelTracks_cfi import *
from Appeltel.PixelTracksRun2010.HIInvertedLowPtPixelTracks_cfi import *

lowPtPixelTrackReco = cms.Sequence(
    siPixelRecHits * 
    hiLowPtPixelTracks *
    hiGoodPixelTracks *
    hiLowPtWideCutPixelTracks *
    hiInvertedLowPtPixelTracks 
)

goodPixelTrackReco = cms.Sequence(
    siPixelRecHits *
    hiGoodPixelTracks
)

conformalPixelTrackReco = cms.Sequence(
    siPixelRecHits *
    hiConformalPixelTracks 
)

#
# This code is for secondary pixel tracking by 
# cluster removal. Please leave in the file as
# comments as we may choose to use it later.
#
# If making conformal pixel tracks, make sure that
# you set the quality bit for hiGoodTracks!
#
# You will also need to rerun heavy ion tracking
# reconstruction in order to propagate the trajectories
#
#from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi import *
#from RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi import *
#
#firstStepFilter = cms.EDProducer("QualityFilter",
#    TrackQuality = cms.string('tight'),
#    recTracks = cms.InputTag("hiGoodTracks")
#)
#
#secClusters = cms.EDProducer("TrackClusterRemover",
#    trajectories = cms.InputTag("firstStepFilter"),
#    pixelClusters = cms.InputTag("siPixelClusters"),
#    stripClusters = cms.InputTag("siStripClusters"),
#    Common = cms.PSet(
#        maxChi2 = cms.double(1000.0)
#    )
#)
#
#secPixelRecHits = siPixelRecHits.clone(
#    src = cms.InputTag("secClusters")
#    )
#
#secStripRecHits = siStripMatchedRecHits.clone(
#    ClusterProducer = 'secClusters'
#    )
#
#secpixellayertriplets = pixellayertriplets.clone()
#secpixellayertriplets.ComponentName = cms.string('SecPixelLayerTriplets')
#secpixellayertriplets.BPix.HitProducer = cms.string('secPixelRecHits')
#secpixellayertriplets.FPix.HitProducer = cms.string('secPixelRecHits')
#
#hiConformalPixelTracks.OrderedHitsFactoryPSet.SeedingLayers = cms.string("SecPixelLayerTriplets")
#
#conformalPixelTrackReco = cms.Sequence(
#    siPixelRecHits *
#    siStripMatchedRecHits *
#    firstStepFilter *
#    secClusters *
#    secPixelRecHits *
#    hiConformalPixelTracks 
#)
