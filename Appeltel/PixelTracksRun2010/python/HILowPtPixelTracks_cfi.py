import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
from RecoHI.HiTracking.HIPixelTrackFilter_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtPixelTrackFilter_cfi import *
from RecoHI.HiTracking.HITrackingRegionProducer_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtTrackingRegionProducer_cfi import *

hiLowPtPixelTracks = cms.EDProducer("PixelTrackProducer",

    passLabel  = cms.string('Pixel triplet low-pt tracks with vertex constraint'),

    # Region
    RegionFactoryPSet = cms.PSet(
	  ComponentName = cms.string("GlobalTrackingRegionWithVerticesProducer"),
	  RegionPSet = cms.PSet(
                HiLowPtTrackingRegionWithVertexBlock
	  )
    ),
     
    # Ordered Hits
    OrderedHitsFactoryPSet = cms.PSet( 
          ComponentName = cms.string( "StandardHitTripletGenerator" ),
	  SeedingLayers = cms.string( "PixelLayerTriplets" ),
          GeneratorPSet = cms.PSet( 
		PixelTripletHLTGenerator
          )
    ),
	
    # Fitter
    FitterPSet = cms.PSet( 
	  ComponentName = cms.string('PixelFitterByHelixProjections'),
	  TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets')
    ),
	
    # Filter
    useFilterWithES = cms.bool( True ),
    FilterPSet = cms.PSet( 
          HiLowPtFilterBlock
    ),
	
    # Cleaner
    CleanerPSet = cms.PSet(  
          ComponentName = cms.string( "TrackCleaner" )
    )
)

# increase threshold for triplets in generation step (default: 10000)
hiLowPtPixelTracks.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000

