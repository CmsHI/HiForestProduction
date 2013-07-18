import FWCore.ParameterSet.Config as cms

towersAboveThreshold = cms.EDProducer("CaloTowerCandidateCreator",
    src = cms.InputTag("towerMaker"),
    verbose = cms.untracked.int32(0),
    minimumE = cms.double(3.0),
    minimumEt = cms.double(0.0),
)

hfPosTowers = cms.EDFilter("EtaPtMinCandSelector",
    src = cms.InputTag("towersAboveThreshold"),
    ptMin   = cms.double(0),
    etaMin = cms.double(3.0),
    etaMax = cms.double(6.0)
)

hfNegTowers = cms.EDFilter("EtaPtMinCandSelector",
    src = cms.InputTag("towersAboveThreshold"),
    ptMin   = cms.double(0),
    etaMin = cms.double(-6.0),
    etaMax = cms.double(-3.0)
)


hfPosFilter = cms.EDFilter("CandCountFilter",
    src = cms.InputTag("hfPosTowers"),
    minNumber = cms.uint32(1)
)

hfNegFilter = cms.EDFilter("CandCountFilter",
    src = cms.InputTag("hfNegTowers"),
    minNumber = cms.uint32(1)
)

hfEitherTowers = cms.EDProducer("CandMerger",
    src = cms.VInputTag( "hfPosTowers", "hfNegTowers" )
)

hfEitherFilter = cms.EDFilter("CandCountFilter",
    src = cms.InputTag("hfEitherTowers"),
    minNumber = cms.uint32(1) 
)

hfCoincFilter = cms.Sequence(
        towersAboveThreshold *
        hfPosTowers *
        hfNegTowers *
        hfPosFilter *
        hfNegFilter)


hfEitherOfTwoFilter = cms.Sequence(
    towersAboveThreshold *
    hfPosTowers *
    hfNegTowers *
    hfEitherTowers* 
    hfEitherFilter)

