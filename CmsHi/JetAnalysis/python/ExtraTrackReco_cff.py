import FWCore.ParameterSet.Config as cms

# pixel triplet tracking (HI Tracking)
from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
from RecoHI.Configuration.Reconstruction_HI_cff import *

#Track Reco
rechits = cms.Sequence(siPixelRecHits * siStripMatchedRecHits)
hiTrackReReco = cms.Sequence(rechits * heavyIonTracking)

# good track selection
hiCaloCompatibleGeneralTracksQuality = cms.EDFilter("TrackSelector",
                                                       src = cms.InputTag("hiGeneralCaloMatchedTracks"),
                                                       cut = cms.string(
       'quality("highPuritySetWithPV")')
                                                    )

hiGeneralTracksQuality = cms.EDFilter("TrackSelector",
                                         src = cms.InputTag("hiGeneralCaloMatchedTracks"),
                                         cut = cms.string(
       'quality("highPurity")')
                                      )

hiSelectedTrackQuality = cms.EDFilter("TrackSelector",
                                         src = cms.InputTag("hiGeneralCaloMatchedTracks"),
                                         cut = cms.string(
       'quality("highPurity")&&algo()==4')
                                      )




from Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff import *
from Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff import *

hiMergedTracks = hiGoodMergedTracks.clone(
    TrackProducer1  = "hiCaloTracks",
    TrackProducer2  = "hiConformalPixelTracks")

hiCaloTracks = cms.EDFilter("TrackSelector",
                            src = cms.InputTag("hiGeneralCaloMatchedTracks"),
                            cut = cms.string(
    'quality("highPuritySetWithPV")')
                            )

hiTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("hiGeneralCaloMatchedTracks"),
                                cut = cms.string(
    'quality("highPurity")')
                                )

##################################################



hiMergedTracksSelcted = hiGoodMergedTracks.clone(
        TrackProducer1  = "hiSelectedTrackQuality",
            TrackProducer2  = "hiConformalPixelTracks")

hiMergedTracksGeneral = hiGoodMergedTracks.clone(
        TrackProducer1  = "hiGeneralTracksQuality",
            TrackProducer2  = "hiConformalPixelTracks")

hiMergedTracksGeneralCalo = hiGoodMergedTracks.clone(
        TrackProducer1  = "hiCaloCompatibleGeneralTracksQuality",
            TrackProducer2  = "hiConformalPixelTracks")

##################################################

rechits = cms.Sequence(siPixelRecHits * siStripMatchedRecHits)
hiTrackReco = cms.Sequence(hiTracks * hiCaloTracks* hiMergedTracks)
hiTrackDebug = cms.Sequence( hiCaloCompatibleGeneralTracksQuality * hiMergedTracksGeneralCalo *
                                  hiSelectedTrackQuality * hiMergedTracksSelcted *
                                  hiGeneralTracksQuality * hiMergedTracksGeneral
                                  )

