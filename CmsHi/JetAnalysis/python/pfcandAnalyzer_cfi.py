import FWCore.ParameterSet.Config as cms

pfcandAnalyzer = cms.EDAnalyzer('HiPFCandAnalyzer',
                                pfCandidateLabel = cms.InputTag("particleFlowTmp"),
                                jetLabel = cms.InputTag("ak5patJets"),
                                genLabel = cms.InputTag("hiGenParticles"),
                                # debug
                                pfPtMin = cms.double(0.5),
                                genPtMin = cms.double(0.5),
                                jetPtMin = cms.double(20.0),                                
                                verbosity = cms.untracked.int32(0),
                                skipCharged = cms.untracked.bool(False)
                                )

