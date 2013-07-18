import FWCore.ParameterSet.Config as cms

SClusterComparison = cms.EDAnalyzer("SClusterComparison",
                                    cleanedHybridSuperClusters = cms.InputTag("cleanedHybridSuperClusters"),
                                    correctedHybridSuperClusters = cms.InputTag("correctedHybridSuperClusters"),
                                    hybridSuperClusters = cms.InputTag("hybridSuperClusters"),
                                    uncleanedHybridSuperClusters = cms.InputTag("uncleanedHybridSuperClusters"),
                                    #                                    islandBarrelSuperClusters = cms.InputTag("islandSuperClusters","islandBarrelSuperClusters"),
                                    cleanedIslandBarrelSuperClusters = cms.InputTag("spikeCleanFilter"),
                                    correctedIslandBarrelSuperClusters = cms.InputTag("correctedIslandBarrelSuperClusters")
                                    )

