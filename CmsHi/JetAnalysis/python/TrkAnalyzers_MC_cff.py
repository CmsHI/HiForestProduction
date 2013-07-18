import FWCore.ParameterSet.Config as cms

from CmsHi.JetAnalysis.TrkAnalyzers_cff import *
from CmsHi.JetAnalysis.TrkEfficiency_cff import *

anaTrack.doSimTrack = True
anaTrack.fillSimTrack = cms.untracked.bool(False)
anaTrack.simTrackPtMin = 1

pixelTrack.doSimTrack = True
pixelTrack.simTrackPtMin = 1
pixelTrack.fillSimTrack = cms.untracked.bool(False)

mergedTrack.doSimTrack = True
mergedTrack.simTrackPtMin = 1
mergedTrack.fillSimTrack = cms.untracked.bool(True)

anaTrack.tpFakeSrc = cms.untracked.InputTag("cutsTPForFak")
anaTrack.tpEffSrc = cms.untracked.InputTag("cutsTPForEff")
pixelTrack.tpFakeSrc = cms.untracked.InputTag("cutsTPForFak")
pixelTrack.tpEffSrc = cms.untracked.InputTag("cutsTPForEff")
anaTrack.tpFakeSrc = cms.untracked.InputTag("cutsTPForFak")
anaTrack.tpEffSrc = cms.untracked.InputTag("cutsTPForEff")
pixelTrack.tpFakeSrc = cms.untracked.InputTag("cutsTPForFak")
pixelTrack.tpEffSrc = cms.untracked.InputTag("cutsTPForEff")


cutsTPForFak.ptMin = 0.45
cutsTPForEff.ptMin = 0.45


                                                                          
