import FWCore.ParameterSet.Config as cms

hiRecHitSpikeCleaner = cms.EDProducer("HiRecHitSpikeCleaner",
                                recHitProducerBarrel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                ebOutputColl  = cms.string( "" ),
                                etCut          = cms.double(20),
                                TimingCut    = cms.untracked.double(4.0),
                                swissCutThr    = cms.untracked.double(0.90),
                                avoidIeta85 = cms.untracked.bool(False)
                                )
