import FWCore.ParameterSet.Config as cms

hiTrackListMerger = cms.EDProducer("HiTrackListMerger",
    # minimum shared fraction to be called duplicate
    ShareFrac = cms.double(0.19),
    # minimum pT in GeV/c
    MinPT = cms.double(0.05),
    # minimum difference in rechit position in cm
    # negative Epsilon uses sharedInput for comparison
    Epsilon = cms.double(-0.001),
    # maximum chisq/dof
    MaxNormalizedChisq = cms.double(1000.0),
    # minimum number of RecHits used in fit
    MinFound = cms.int32(3),
    # module laber of RS Tracks from KF with material propagator
    TrackProducer2 = cms.string(''),
    # module laber of CTF Tracks from KF with material propagator
    TrackProducer1 = cms.string(''),
    # set new quality for confirmed tracks
    promoteTrackQuality = cms.bool(False),
    allowFirstHitShare = cms.bool(False),
    newQuality = cms.string('confirmed'),
    copyExtras = cms.untracked.bool(False),
    trackCollection1MinPt = cms.double(0.),
    trackCollection1MaxPt = cms.double(1000000.),
    trackCollection2MinPt = cms.double(0.),
    trackCollection2MaxPt = cms.double(1000000.)
)

