import FWCore.ParameterSet.Config as cms

hivertexanalyzer = cms.EDAnalyzer('HiVertexAnalyzer',
      vtxlabel=cms.untracked.InputTag("hiSelectedVertex"),
      nTracksBins=cms.untracked.int32(120),
      neededCentBins = cms.untracked.vint32(0,1,3,11,19,39)
)
