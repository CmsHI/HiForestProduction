import FWCore.ParameterSet.Config as cms

# keep only vertices that pass string cut
goodVertices = cms.EDFilter("GoodVertexSelection",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("ndof>4 && abs(z)<=15.0 && abs(position.rho)<=2.0")
)

sortedPreGoodVertices = cms.EDFilter("BestVertexSelection",
       src = cms.InputTag("goodVertices"),
       maxNumber = cms.uint32(100),
       filter = cms.bool(False)
)

sortedGoodVertices = cms.EDFilter("BestVertexSelection",
    src = cms.InputTag("goodVertices"),
    maxNumber = cms.uint32(100),
    filter = cms.bool(False)
)

sortedSumPtOrdVertices = cms.EDFilter("SumPtOrdVertexSelector",
     src = cms.InputTag("goodVertices"),
     maxNumber = cms.uint32(100),
     filter = cms.bool(False)
)

