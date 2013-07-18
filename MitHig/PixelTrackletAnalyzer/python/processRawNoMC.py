import FWCore.ParameterSet.Config as cms

process = cms.Process("myRECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
#process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
process.Timing = cms.Service("Timing") 

process.GlobalTag.globaltag = 'MC_31X_V9::All'

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
"file:Pythia_MinBias_900GeV_HLT8e29_d20091119_000001.root"
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   moduleSeeds = cms.PSet(simMuonRPCDigis = cms.untracked.uint32(49835),
                                                                          simEcalUnsuppressedDigis = cms.untracked.uint32(49835),
                                                                          simSiStripDigis = cms.untracked.uint32(49835),
                                                                          mix = cms.untracked.uint32(49835),
                                                                          simHcalUnsuppressedDigis = cms.untracked.uint32(49835),
                                                                          simMuonCSCDigis = cms.untracked.uint32(49835),
                                                                          VtxSmeared = cms.untracked.uint32(49835),
                                                                          g4SimHits = cms.untracked.uint32(49835),
                                                                          simMuonDTDigis = cms.untracked.uint32(49835),
                                                                          simSiPixelDigis = cms.untracked.uint32(49835)
                                                                          ),
                                                   sourceSeed = cms.untracked.uint32(49835)
                                                   )

#process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertices'),
                             trackSrc = cms.InputTag('pixel3ProtoTracks'),
                             doTrackingParticle = cms.untracked.bool(False)
                             )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('testPythia7TeV.root.hist')
                                   )

process.output = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('000001_2.root'),
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

#process.siPixelRecHits = cms.EDFilter("SiPixelRecHitConverter",
#    src = cms.InputTag("siPixelClusters"),
#    CPE = cms.string('PixelCPEGeneric'),
#    #    CPE = cms.string('PixelCPETemplateReco'),
#    VerboseLevel = cms.untracked.int32(0),
#)
            
#process.reconstruct = cms.Path(process.mix*process.trackingParticles*process.simSiPixelDigis*process.siPixelRawData*process.siPixelDigis*process.pixeltrackerlocalreco*process.offlineBeamSpot*process.pixel3ProtoTracks*process.pixelVertices)
process.reconstruct = cms.Path(process.mix*process.simSiPixelDigis*process.siPixelRawData*process.siPixelDigis*process.pixeltrackerlocalreco*process.offlineBeamSpot*process.pixel3ProtoTracks*process.pixelVertices)

process.analyze = cms.Path(process.ana)
#process.save = cms.EndPath(process.output)


