import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
process.Timing = cms.Service("Timing") 

# Pixel Track and Pixel Vertexing
from RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi import *
#from RecoPixelVertexing.PixelLowPtUtilities.MinBiasTracking_cff import *
#from RecoPixelVertexing.PixelLowPtUtilities.firstStep_cff import *

import RecoPixelVertexing.PixelLowPtUtilities.AllPixelTracks_cfi
pixel3ProtoTracks = RecoPixelVertexing.PixelLowPtUtilities.AllPixelTracks_cfi.allPixelTracks.clone()

pixelVertices.TrackCollection = 'pixel3ProtoTracks'
pixelVertices.UseError = True
pixelVertices.WtAverage = True
pixelVertices.ZOffset = 5.
pixelVertices.ZSeparation = 0.3
pixelVertices.NTrkMin = 3
pixelVertices.PtMin = 0.15

process.pixel3ProtoTracks = pixel3ProtoTracks
process.pixelVertices = pixelVertices

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'dcache:__INPUT__'
    )
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   moduleSeeds = cms.PSet(simMuonRPCDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simEcalUnsuppressedDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simSiStripDigis = cms.untracked.uint32(__RANDOM__),
                                                                          mix = cms.untracked.uint32(__RANDOM__),
                                                                          simHcalUnsuppressedDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simMuonCSCDigis = cms.untracked.uint32(__RANDOM__),
                                                                          VtxSmeared = cms.untracked.uint32(__RANDOM__),
                                                                          g4SimHits = cms.untracked.uint32(__RANDOM__),
                                                                          simMuonDTDigis = cms.untracked.uint32(__RANDOM__),
                                                                          simSiPixelDigis = cms.untracked.uint32(__RANDOM__)
                                                                          ),
                                                   sourceSeed = cms.untracked.uint32(__RANDOM__)
                                                   )

process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertices')
                             )

process.ecalloc = cms.Sequence(process.ecalWeightUncalibRecHit*process.ecalRecHit)

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('__OUTPUT__.hist')
                                   )

process.output = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('__OUTPUT__'),
                                  outputCommands = cms.untracked.vstring('keep *',
                                                                         'drop *Digi*_*_*_*',
                                                                         'drop *_*Digi*_*_*',
                                                                         'drop *CrossingFrame*_*_*_*',
                                                                         'drop *Raw*_*_*_*',
                                                                         'keep *_siPixelDigis_*_RECO',
                                                                         'keep *_siPixelRawData_*_RECO',
                                                                         'keep *_simSiPixelDigis_*_RECO',
                                                                         'keep *CrossingFrame_*_g4SimHitsTrackerHits*_RECO'
                                                                         )
                                  )

process.reconstruct = cms.Path(process.mix*process.trackingParticles*process.simSiPixelDigis*process.siPixelRawData*process.siPixelDigis*process.pixeltrackerlocalreco*process.offlineBeamSpot*process.pixel3ProtoTracks*process.pixelVertices)

process.analyze = cms.Path(process.ana)
#process.save = cms.EndPath(process.output)


