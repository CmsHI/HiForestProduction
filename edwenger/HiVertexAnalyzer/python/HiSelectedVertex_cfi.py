import FWCore.ParameterSet.Config as cms

# sort by number of tracks and keep the best
hiselectedVertex = cms.EDFilter("HiBestVertexSelection",
    src = cms.InputTag("pixel3Vertices"),
    maxNumber = cms.uint32(1),
    filter = cms.bool(True) 
)
