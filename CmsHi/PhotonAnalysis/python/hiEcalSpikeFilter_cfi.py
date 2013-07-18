import FWCore.ParameterSet.Config as cms

hiEcalSpikeFilter = cms.EDFilter("HiEcalSpikeFilter",
                                 photonProducer = cms.InputTag("photons"),
                                 ebReducedRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                 eeReducedRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
                                 minEt = cms.double(10),  # 10GeV rechits
                                 swissThreshold = cms.double(0.85),
                                 avoidIeta85 = cms.bool(True), 
                                 timeThreshold = cms.double(4) # 4 ns
                                 )




