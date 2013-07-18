import FWCore.ParameterSet.Config as cms

# reco jet with russian pileup subtraction
from RecoHI.HiJetAlgos.HiRecoJets_cff import *
from RecoHI.HiJetAlgos.HiRecoPFJets_cff import *

akPu3CaloJets = cms.EDProducer(
    "FastjetJetProducer",
    HiCaloJetParameters,
    AnomalousCellParameters,
    MultipleAlgoIteratorBlock,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.3)
    )
akPu5CaloJets.radiusPU = 0.5

iterativeConePu5CaloJets.doPVCorrection = cms.bool(True)
iterativeConePu5CaloJets.srcPVs = 'hiSelectedVertex'

iterativeConePu5CaloJets.jetPtMin = 1.0
ak5PFJets.jetPtMin = 1.0

akPu5PFJets = ak5PFJets.clone()
akPu5PFJets.src = 'PFTowers'
akPu5PFJets.jetType = 'BasicJet'
akPu5PFJets.doPUOffsetCorr = True
akPu5PFJets.sumRecHits = False

akPu3PFJets = akPu5PFJets.clone()
akPu3PFJets.rParam = cms.double(0.3)

# pileup subtraction jet exclusion pt min
iterativeConePu5CaloJets.puPtMin = cms.double(10.0)
akPu5PFJets.puPtMin = cms.double(25.0)
akPu3PFJets.puPtMin = cms.double(15.0)
akPu5CaloJets.puPtMin = cms.double(10.0)
akPu3CaloJets.puPtMin = cms.double(10.0)

akPu5PFJets.doRhoFastjet = False
akPu5PFJets.doAreaFastjet = False

akPu3PFJets.doRhoFastjet = False
akPu3PFJets.doAreaFastjet = False

akPu5CaloJets.doRhoFastjet = False
akPu5CaloJets.doAreaFastjet = False
akPu5CaloJets.doPUOffsetCorr = True

akPu3CaloJets.doRhoFastjet = False
akPu3CaloJets.doAreaFastjet = False
akPu3CaloJets.doPUOffsetCorr = True


### Extra extended algos & sequence
akPu6PFJets = akPu3PFJets.clone(rParam = 0.6)
akPu6CaloJets = akPu3CaloJets.clone(rParam = 0.6)
akPu4PFJets = akPu3PFJets.clone(rParam = 0.4)
akPu4CaloJets = akPu3CaloJets.clone(rParam = 0.4)
akPu2PFJets = akPu3PFJets.clone(rParam = 0.2)
akPu2CaloJets = akPu3CaloJets.clone(rParam = 0.2)
akPu1PFJets = akPu3PFJets.clone(rParam = 0.1)
akPu1CaloJets = akPu3CaloJets.clone(rParam = 0.1)

iterativeCone5CaloJets = iterativeConePu5CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak6PFJets = akPu6PFJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak6CaloJets = akPu6CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak5PFJets = akPu5PFJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak5CaloJets = akPu5CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak4PFJets = akPu4PFJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak4CaloJets = akPu4CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak3PFJets = akPu3PFJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak3CaloJets = akPu3CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 10)
ak2PFJets = akPu2PFJets.clone(doPUOffsetCorr = False, jetPtMin = 5)
ak2CaloJets = akPu2CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 5)
ak1PFJets = akPu1PFJets.clone(doPUOffsetCorr = False, jetPtMin = 5)
ak1CaloJets = akPu1CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 5)


akPu1PFJets.puPtMin = cms.double(5.0)
akPu2PFJets.puPtMin = cms.double(10.0)
akPu3PFJets.puPtMin = cms.double(15.0)
akPu4PFJets.puPtMin = cms.double(20.0)
akPu5PFJets.puPtMin = cms.double(25.0)
akPu6PFJets.puPtMin = cms.double(30.0)

akPu1CaloJets.puPtMin = cms.double(2.0)
akPu2CaloJets.puPtMin = cms.double(4.0)
akPu3CaloJets.puPtMin = cms.double(6.0)
akPu4CaloJets.puPtMin = cms.double(8.0)
akPu5CaloJets.puPtMin = cms.double(10.0)
akPu6CaloJets.puPtMin = cms.double(12.0)


recoAk1to6 = cms.Sequence( akPu1PFJets * akPu2PFJets *akPu3PFJets * akPu4PFJets * akPu5PFJets * akPu6PFJets *
                           ak1PFJets * ak2PFJets *ak3PFJets * ak4PFJets * ak5PFJets * ak6PFJets *
                           akPu1CaloJets * akPu2CaloJets *akPu3CaloJets * akPu4CaloJets * akPu5CaloJets * akPu6CaloJets *
                           ak1CaloJets * ak2CaloJets *ak3CaloJets * ak4CaloJets * ak5CaloJets * ak6CaloJets
                           )








