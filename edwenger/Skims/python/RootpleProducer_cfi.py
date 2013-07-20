import FWCore.ParameterSet.Config as cms

## ROOTple
rootpleProducer = cms.EDProducer("AnalysisRootpleProducer",
                                 TracksCollectionName      = cms.InputTag("goodTracks"),
                                 VertexCollectionName      = cms.InputTag("offlinePrimaryVertices"),
                                 MCEvent                   = cms.InputTag("generator"),
                                 triggerEvent              = cms.InputTag("hltTriggerSummaryAOD"),
                                 ChgGenPartCollectionName  = cms.InputTag("chargeParticles"),
                                 RecoCaloJetCollectionName = cms.InputTag("ak7CaloJets"),
                                 ChgGenJetCollectionName   = cms.InputTag("ak7GenJets"),
                                 TracksJetCollectionName   = cms.InputTag("ak7PFJets"),
                                 GenJetCollectionName      = cms.InputTag("ak7GenJets"),
                                 OnlyRECO                  = cms.bool(True),
                                 triggerResults            = cms.InputTag("TriggerResults","","HLT"),
                                 genEventScale             = cms.InputTag("generator")
                                 )


## genparticle and genjet selections for filling in MC
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
from RecoJets.Configuration.GenJetParticles_cff import *
from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *
goodParticles = cms.EDFilter("GenParticleSelector",
    filter = cms.bool(False),
    src = cms.InputTag("genParticles"),
    cut = cms.string('pt > 0.0 & status == 1')
)

chargeParticles = cms.EDFilter("GenParticleSelector",
    filter = cms.bool(False),
    src = cms.InputTag("genParticles"),
    cut = cms.string('charge != 0 & pt > 0.0 & status == 1')
)

AnalysisParticles = cms.Sequence(genParticles*genJetParticles*goodParticles*chargeParticles)

rootpleProducerMC = cms.Sequence(AnalysisParticles * rootpleProducer)
