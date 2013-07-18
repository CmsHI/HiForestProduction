import FWCore.ParameterSet.Config as cms

vtxconstraintprod = cms.EDAnalyzer('VertexConstraintProducer',
                                   traksrc = cms.untracked.InputTag("hiSelectedTracks"),
                                   vtxsrc = cms.untracked.InputTag("hiSelectedVertex")
                                   )
