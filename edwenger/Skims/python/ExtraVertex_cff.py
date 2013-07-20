import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Reconstruction_cff import *
from RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi import *
import RecoPixelVertexing.PixelLowPtUtilities.AllPixelTracks_cfi

############################
# Pixel-3 proto tracks
pixel3ProtoTracks = RecoPixelVertexing.PixelLowPtUtilities.AllPixelTracks_cfi.allPixelTracks.clone()
pixel3ProtoTracks.passLabel = 'Pixel triplet tracks for vertexing'
pixel3ProtoTracks.RegionFactoryPSet.ptMin = cms.double(0.2)

############################
# Vertexing
import UserCode.FerencSiklerVertexing.NewVertexProducer_cfi
pixel3Vertices = UserCode.FerencSiklerVertexing.NewVertexProducer_cfi.newVertices.clone()
pixel3Vertices.TrackCollection = 'pixel3ProtoTracks'
pixel3Vertices.dMax = cms.double(12.0)

extraVertex = cms.Sequence(
    siPixelRecHits *
    pixel3ProtoTracks *
    pixel3Vertices
)
