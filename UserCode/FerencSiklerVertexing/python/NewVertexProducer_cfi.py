import FWCore.ParameterSet.Config as cms

newVertices = cms.EDProducer("NewVertexProducer",
    TrackCollection = cms.string('generalTracks'),
    PtMin    = cms.double(0.1),  # GeV/c
    nSigma   = cms.double(3.0),  # dt < nSigma * sigma, transverse impact
    Method   = cms.int32(0),     # 0 = fPNN, 1 = KMeans, 2 = GaussMix
    dMax     = cms.double(12.0), # for fPNN, was 8.0 FIXME
    ProbMin  = cms.double(1e-3), # for kMeans and GaussMix 
    TrackMin = cms.int32(2)      # minimum number of tracks for a vertex
)
