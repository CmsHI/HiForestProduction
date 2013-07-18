import FWCore.ParameterSet.Config as cms

vertexanalyzer = cms.EDAnalyzer('VertexAnalyzer',
      vtxlabel=cms.untracked.InputTag("pixel3Vertices"),
      trklabel=cms.untracked.InputTag("selectTracks"),
      jetlabel=cms.untracked.InputTag("cleanedPatJets"),
      dzcut=cms.untracked.double(0.9), # |dz|<dz_cut, which should be less then min distance between two vertices
      minJetEt=cms.untracked.double(0.0),
      nthleadingMax=cms.untracked.int32(3), # nth leading tracks (n<nthleadingMax)
      nTracksBins=cms.untracked.int32(50),
      jetTrkVerticesCorr=cms.untracked.bool(False)
)
