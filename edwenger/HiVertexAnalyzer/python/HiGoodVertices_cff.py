import FWCore.ParameterSet.Config as cms

# keep only vertices that pass string cut
goodVertices = cms.EDFilter("HiGoodVertexSelection",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("ndof>4 && abs(z)<=15.0 && abs(position.rho)<=2.0")
)

sortedGoodVertices = cms.EDFilter("HiBestVertexSelection",
    src = cms.InputTag("goodVertices"),
    maxNumber = cms.uint32(100),
    filter = cms.bool(False)
)
