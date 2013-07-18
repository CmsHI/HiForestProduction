import FWCore.ParameterSet.Config as cms

pi0Analyzer = cms.EDAnalyzer("pi0Analyzer",
                             pfCandidatesTag  = cms.InputTag("particleFlowTmp",""),
                             isMC = cms.untracked.bool(False), 
                             useCentrality = cms.untracked.bool(True)
                             )


