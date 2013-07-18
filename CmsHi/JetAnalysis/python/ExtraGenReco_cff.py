import FWCore.ParameterSet.Config as cms

# Gen
from PhysicsTools.HepMCCandAlgos.HiGenParticles_cfi import *
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoHI.HiJetAlgos.HiGenJets_cff import *
from RecoHI.HiJetAlgos.HiGenCleaner_cff import *
from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonProducer_cfi import heavyIon

hiGenParticlesForJets.ignoreParticleIDs += cms.vuint32( 12,14,16)

hiGenParticles.srcVector = cms.vstring('hiSignal')

ak1HiGenJets = ak3HiGenJets.clone(rParam = 0.1)
ak2HiGenJets = ak3HiGenJets.clone(rParam = 0.2)
ak6HiGenJets = ak3HiGenJets.clone(rParam = 0.6)

hiAnaGenJets = cms.Sequence(iterativeCone5HiGenJets +
                            ak1HiGenJets +
                            ak2HiGenJets +
                            ak3HiGenJets +
                            ak4HiGenJets +
                            ak5HiGenJets +
                            ak6HiGenJets
                            )
                            

hiGen = cms.Sequence(
  heavyIon * # GenHIEventProducer
  hiGenParticles *
  hiGenParticlesForJets *
  hiAnaGenJets *
  genPartons *
  hiPartons
  )
