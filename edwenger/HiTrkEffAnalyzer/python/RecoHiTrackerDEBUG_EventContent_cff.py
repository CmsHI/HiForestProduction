import FWCore.ParameterSet.Config as cms


#Tracking DEBUG content for efficiency study
RecoHiTrackerDEBUG = cms.PSet(
         outputCommands = cms.untracked.vstring('keep *_hiGlobalPrimTracks_*_*',
                                                'keep *_hiSelectedTracks_*_*',
                                                'keep *_hiRefitPxlTracks_*_*',
                                                'keep *_hiPixel3PrimTracks_*_*',
                                                'keep *_hiPixel3ProtoTracks_*_*',
                                                'keep *_hiLowPtPixelTracks_*_*',
                                                'keep recoVertexs_hiPixelMedianVertex_*_*',
                                                'keep recoVertexs_hiPixelAdaptiveVertex_*_*',
                                                'keep recoVertexs_hiSelectedVertex_*_*',
                                                'keep *_mergedtruth_*_*',
                                                'keep *_hiGenParticles_*_*',
                                                'keep *_simSiPixelDigis_*_*',
                                                'keep *_simSiStripDigis_*_*',
                                                'keep *_siStripClusters_*_*',
                                                'keep *_siPixelClusters_*_*'
                                                )
         )
