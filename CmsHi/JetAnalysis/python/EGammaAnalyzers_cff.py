import FWCore.ParameterSet.Config as cms

from CmsHi.PhotonAnalysis.MultiPhotonAnalyzer_cfi import *
multiPhotonAnalyzer.GenParticleProducer = cms.InputTag("hiGenParticles")
multiPhotonAnalyzer.PhotonProducer = cms.InputTag("selectedPatPhotons")
multiPhotonAnalyzer.VertexProducer = cms.InputTag("hiSelectedVertex")
multiPhotonAnalyzer.OutputFile = cms.string('mpa.root')
multiPhotonAnalyzer.doStoreCompCone = cms.untracked.bool(False)
multiPhotonAnalyzer.doStoreJets = cms.untracked.bool(False)

