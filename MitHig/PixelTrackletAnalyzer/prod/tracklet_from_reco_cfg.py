import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
            "dcache:///pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/reco/__DATATAG__/__DATATAG___r__INPUT__.root"
#            "dcache:///pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/sim/pythia_900GeV_d20080526/pythia_900GeV_d20080526_r__INPUT__.root"

    )
)

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

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

process.ana1 = cms.EDAnalyzer('PixelTrackletAnalyzer',useRecoVertex = cms.untracked.bool(True),
                             skipBestMatch = cms.untracked.int32(1),
                             useDeltaPhi = cms.untracked.bool(True),
                             AssociatorParameters = cms.PSet(associatePixel = cms.bool(True),
                                                             associateStrip = cms.bool(False),
                                                             associateRecoTracks = cms.bool(False),
                                                             ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                                                                                   'g4SimHitsTrackerHitsPixelBarrelHighTof',
                                                                                   'g4SimHitsTrackerHitsPixelEndcapLowTof',
                                                                                   'g4SimHitsTrackerHitsPixelEndcapHighTof'
                                                                                   )
                                                             )
                             
                             )

process.anasim = cms.EDAnalyzer('SimTrackAnalyzer',
                                useRecoVertex = cms.untracked.bool(True),
                                vertexSrc = cms.untracked.string('pixelVertices'),
                                investigateSimTracks = cms.untracked.bool(True),
                                AssociatorParameters = cms.PSet(associatePixel = cms.bool(True),
                                                                associateStrip = cms.bool(False),
                                                                associateRecoTracks = cms.bool(False),
                                                                ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                                                                                      'g4SimHitsTrackerHitsPixelBarrelHighTof',
                                                                                      'g4SimHitsTrackerHitsPixelEndcapLowTof',
                                                                                      'g4SimHitsTrackerHitsPixelEndcapHighTof'
                                                                                      )
                                                                )
                                
                                )



process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                              vertexSrc = cms.vstring('pixelVertices')
                              )

process.ana2 = cms.EDAnalyzer('PixelTrackletAnalyzer',useRecoVertex = cms.untracked.bool(True),
                              skipBestMatch = cms.untracked.int32(1),
                              useDeltaPhi = cms.untracked.bool(True),
                              vertexSrc = cms.untracked.string('pixelVertexFromClusters'),
                              AssociatorParameters = cms.PSet(associatePixel = cms.bool(True),
                                                              associateStrip = cms.bool(False),
                                                              associateRecoTracks = cms.bool(False),
                                                              ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                                                                                    'g4SimHitsTrackerHitsPixelBarrelHighTof',
                                                                                    'g4SimHitsTrackerHitsPixelEndcapLowTof',
                                                                                    'g4SimHitsTrackerHitsPixelEndcapHighTof'
                                                                                    )
                                                              )
                              
                              )


process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('__OUTPUT__.hist.root')
                                   )

#process.reconstruct = cms.Path(process.mix*process.trackingParticles*process.doAllDigi*process.L1Emulator*process.DigiToRaw*process.RawToDigi*process.trackerlocalreco*process.ecalloc*process.offlineBeamSpot*process.pixel3ProtoTracks*process.pixelVertices*process.pixelVertexFromClusters)
process.reconstruct = cms.Path(process.siPixelClusters*process.siPixelRecHits)#*process.pixelVertexFromClusters)
process.analyze = cms.Path(process.ana)
#process.save = cms.EndPath(process.output)


