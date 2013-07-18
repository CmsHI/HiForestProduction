from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoHI.HiJetAlgos.HiPFJetParameters_cff import *
from RecoHI.HiJetAlgos.HiCaloJetParameters_cff import *


akPu3PFConesAna = cms.EDProducer(
       "JetAlgorithmAnalyzer",
          HiPFJetParameters,
          AnomalousCellParameters,
          MultipleAlgoIteratorBlock,
          jetAlgorithm = cms.string("AntiKt"),
          rParam       = cms.double(0.3),
          useInputJets = cms.untracked.bool(True),
          inputJetSrc = cms.InputTag("akPu3PFpatJets"),
          useTowersForBkg = cms.bool(True),
          centralityTag = cms.InputTag("hiCentrality"),
          evtPlaneTag = cms.InputTag("hiEvtPlane","recoLevel"),
          avoidNegative = cms.bool(False),
          patJetSrc = cms.untracked.InputTag("akPu3PFpatJets"),
       evtPlaneIndex = cms.untracked.int32(21),
       doBackToBack  = cms.untracked.bool(True),
       centrality  = cms.untracked.int32(-1)
          )

akPu3PFConesAna.doPUOffsetCorr = True

akPu3PFConesAna.jetType = 'BasicJet'
akPu3PFConesAna.src = cms.InputTag("PFTowers")

akPu5PFConesAna = akPu3PFConesAna.clone(
          rParam       = cms.double(0.5),
          )

akPu5CaloConesAna = akPu3PFConesAna.clone(
    src = cms.InputTag("towerMaker"),
    rParam       = cms.double(0.5)
    )

akPu3CaloConesAna = akPu3PFConesAna.clone(
    src = cms.InputTag("towerMaker"),
    rParam       = cms.double(0.3)
    )

icPu5CaloConesAna = cms.EDProducer(
       "JetAlgorithmAnalyzer",
          HiCaloJetParameters,
          AnomalousCellParameters,
          MultipleAlgoIteratorBlock,
          jetAlgorithm = cms.string("IterativeCone"),
          rParam       = cms.double(0.5),
          useInputJets = cms.untracked.bool(True),
          inputJetSrc = cms.InputTag("iterativeConePu5CaloJets"),
          useTowersForBkg = cms.bool(True),
          centralityTag = cms.InputTag("hiCentrality"),
          evtPlaneTag = cms.InputTag("hiEvtPlane","recoLevel"),

          avoidNegative = cms.bool(False),
          patJetSrc = cms.untracked.InputTag("icPu5patJets")

          )


akPu5PFConesAna.puPtMin = cms.double(25)
akPu3PFConesAna.puPtMin = cms.double(15)

akPu3CaloConesAna.puPtMin = cms.double(10)
akPu5CaloConesAna.puPtMin = cms.double(10)

icPu5CaloConesAna.puPtMin = cms.double(10)

icPu5CaloConesAna.doPUOffsetCorr = True
icPu5CaloConesAna.jetType = 'BasicJet'
icPu5CaloConesAna.src = cms.InputTag("towerMaker")



randomCones = cms.Sequence(
    akPu3PFConesAna+
    akPu5PFConesAna+
    akPu3CaloConesAna+
    akPu5CaloConesAna+
    icPu5CaloConesAna
    )

