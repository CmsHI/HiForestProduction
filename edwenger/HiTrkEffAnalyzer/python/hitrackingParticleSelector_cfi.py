import FWCore.ParameterSet.Config as cms

hitrackingParticleSelector = cms.EDFilter("HackedTrackingParticleSelector",
                             src = cms.InputTag("mergedtruth","MergedTrackTruth"),
                             chargedOnly = cms.bool(True),
                             primaryOnly = cms.bool(False),
                             pdgId = cms.vint32(),
                             tip = cms.double(999.0), # completely open 
                             signalOnly = cms.bool(False),
                             minRapidity = cms.double(-2.5),
                             lip = cms.double(999.0), # completely open
                             ptMin = cms.double(0.8),
                             maxRapidity = cms.double(2.5),
                             minHit = cms.int32(0)
                             )
