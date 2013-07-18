import FWCore.ParameterSet.Config as cms

jetfilter = cms.EDFilter('JetEtaPhiFilter',
                         jsrc = cms.untracked.InputTag("iterativeCone5CaloJets"),
                         minEta = cms.untracked.double(-6.5),
                         maxEta = cms.untracked.double(6.5),
                         minPhi = cms.untracked.double(-3.2),
                         maxPhi = cms.untracked.double(3.2)
                         )
