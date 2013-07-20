import FWCore.ParameterSet.Config as cms

# keep only vertices that pass string cut
goodVertices = cms.EDFilter("GoodVertexSelection",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("ndof>4")
)

pileupEventFilter = cms.EDFilter("VertexCountFilter",
    src = cms.InputTag("goodVertices"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    filter = cms.bool(True)
 )                                
