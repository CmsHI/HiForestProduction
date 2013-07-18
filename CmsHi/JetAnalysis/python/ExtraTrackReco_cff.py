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




hiTracks = cms.EDFilter("TrackSelector",
                        src = cms.InputTag("hiGeneralCaloMatchedTracks"),
                        cut = cms.string(
    'quality("highPurity")')
                        )



