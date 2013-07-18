from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoHI.HiJetAlgos.HiPFJetParameters_cff import *
from RecoHI.HiJetAlgos.HiCaloJetParameters_cff import *

akPu1PFCones = cms.EDProducer(
    "JetAlgorithmAnalyzer",
    HiPFJetParameters,
    AnomalousCellParameters,
    MultipleAlgoIteratorBlock,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.1),
    useInputJets = cms.untracked.bool(True),
    useTowersForBkg = cms.bool(True),
    centralityTag = cms.InputTag("pACentrality"),
    evtPlaneTag = cms.InputTag("hiEvtPlane","recoLevel"),
    avoidNegative = cms.bool(False),
    patJetSrc = cms.untracked.InputTag("akPu1PFpatJets"),
    evtPlaneIndex = cms.untracked.int32(21),
    doBackToBack  = cms.untracked.bool(True),
    centrality  = cms.untracked.int32(-1),
    doRecoEvtPlane  = cms.untracked.bool(False),        
    doAnalysis  = cms.untracked.bool(False),
    puPtMin = cms.double(15.0)    
    )

akPu1PFCones.doPUOffsetCorr = True
akPu1PFCones.jetType = 'BasicJet'
akPu1PFCones.src = cms.InputTag("PFTowers")

akPu2PFCones = akPu1PFCones.clone(rParam = 0.2, patJetSrc = cms.untracked.InputTag("akPu2PFpatJets"))
akPu3PFCones = akPu1PFCones.clone(rParam = 0.3, patJetSrc = cms.untracked.InputTag("akPu3PFpatJets"))
akPu4PFCones = akPu1PFCones.clone(rParam = 0.4, patJetSrc = cms.untracked.InputTag("akPu4PFpatJets"))
akPu5PFCones = akPu1PFCones.clone(rParam = 0.5, patJetSrc = cms.untracked.InputTag("akPu5PFpatJets"))
akPu6PFCones = akPu1PFCones.clone(rParam = 0.6, patJetSrc = cms.untracked.InputTag("akPu6PFpatJets"))

akPu1CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = 0.1, patJetSrc = cms.untracked.InputTag("akPu1CalopatJets"))
akPu2CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = 0.2, patJetSrc = cms.untracked.InputTag("akPu2CalopatJets"))
akPu3CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = 0.3, patJetSrc = cms.untracked.InputTag("akPu3CalopatJets"))
akPu4CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = 0.4, patJetSrc = cms.untracked.InputTag("akPu4CalopatJets"))
akPu5CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = 0.5, patJetSrc = cms.untracked.InputTag("akPu5CalopatJets"))
akPu6CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = 0.6, patJetSrc = cms.untracked.InputTag("akPu6CalopatJets"))


akPu1PFCones.puPtMin = 5
akPu2PFCones.puPtMin = 10
akPu3PFCones.puPtMin = 15
akPu4PFCones.puPtMin = 20
akPu5PFCones.puPtMin = 25
akPu6PFCones.puPtMin = 30

akPu1CaloCones.puPtMin = 2
akPu2CaloCones.puPtMin = 4
akPu3CaloCones.puPtMin = 6
akPu4CaloCones.puPtMin = 8
akPu5CaloCones.puPtMin = 10
akPu6CaloCones.puPtMin = 12

randomCones = cms.Sequence(
    akPu2PFCones+
    akPu3PFCones+
    akPu4PFCones+
    akPu5PFCones+
    akPu6PFCones+
    akPu2CaloCones+
    akPu3CaloCones+
    akPu4CaloCones+
    akPu5CaloCones+
    akPu6CaloCones
    )

randomCones2to5 = cms.Sequence(
    akPu2PFCones+
    akPu3PFCones+
    akPu4PFCones+
    akPu5PFCones+
    akPu2CaloCones+
    akPu3CaloCones+
    akPu4CaloCones+
    akPu5CaloCones
    )

