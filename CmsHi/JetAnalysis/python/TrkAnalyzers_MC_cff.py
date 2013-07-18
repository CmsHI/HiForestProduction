import FWCore.ParameterSet.Config as cms

from CmsHi.JetAnalysis.TrkAnalyzers_cff import *
from CmsHi.JetAnalysis.TrkEfficiency_cff import *

anaTrack.doSimVertex = True
anaTrack.doSimTrack = True
anaTrack.fillSimTrack = cms.untracked.bool(True)
anaTrack.simTrackPtMin = 1

pixelTrack.doSimVertex = True
pixelTrack.doSimTrack = True
pixelTrack.simTrackPtMin = 1
pixelTrack.fillSimTrack = cms.untracked.bool(True)

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
                                                                         
ppTrack.doSimVertex = True
ppTrack.doSimTrack = True
ppTrack.fillSimTrack = cms.untracked.bool(True)
ppTrack.simTrackPtMin = 1
ppTrack.tpFakeSrc = cms.untracked.InputTag("cutsTPForFak")
ppTrack.tpEffSrc = cms.untracked.InputTag("cutsTPForEff")
ppTrack.tpFakeSrc = cms.untracked.InputTag("cutsTPForFak")
ppTrack.tpEffSrc = cms.untracked.InputTag("cutsTPForEff")
