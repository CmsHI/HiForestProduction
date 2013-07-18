import FWCore.ParameterSet.Config as cms

# global tracking region for lowpt pixel tracks
HiLowPtTrackingRegionWithVertexBlock = cms.PSet(
    ptMin         = cms.double(0.2),	  
    originRadius  = cms.double(0.2),
    nSigmaZ       = cms.double(3.0),		
    beamSpot      = cms.InputTag("offlineBeamSpot"),
    precise       = cms.bool(True),		
    useFoundVertices = cms.bool(True),
    VertexCollection = cms.InputTag("hiSelectedVertex"),		
    useFixedError = cms.bool(True),
    fixedError    = cms.double(0.2),
    sigmaZVertex  = cms.double(3.0)		
    )

