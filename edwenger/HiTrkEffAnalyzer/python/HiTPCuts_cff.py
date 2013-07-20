import FWCore.ParameterSet.Config as cms

import edwenger.HiTrkEffAnalyzer.hitrackingParticleSelector_cfi

cutsTPForFak = edwenger.HiTrkEffAnalyzer.hitrackingParticleSelector_cfi.hitrackingParticleSelector.clone()
cutsTPForEff = cutsTPForFak.clone(primaryOnly = cms.bool(True))

# Higher pt threshold 
cutsTPForFakHigh = edwenger.HiTrkEffAnalyzer.hitrackingParticleSelector_cfi.hitrackingParticleSelector.clone(ptMin = cms.double(4.0))

cutsTPForFakPxl = edwenger.HiTrkEffAnalyzer.hitrackingParticleSelector_cfi.hitrackingParticleSelector.clone(ptMin = cms.double(0.1))
cutsTPForEffPxl = cutsTPForFakPxl.clone(primaryOnly = cms.bool(True))

