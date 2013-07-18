import FWCore.ParameterSet.Config as cms

process = cms.Process("myRECO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# Timing service
process.Timing = cms.Service("Timing") 

process.GlobalTag.globaltag = 'START311_V2::All'

process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/kimy/AllPhysics2760/ppRun_2760GeV_photon10Trig_promptRECO-v3/dd648448e93dc1238d7d4e1d0b5a9841/___outf_edm____99_1_rCe.root'
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


process.ana = cms.EDAnalyzer('PixelHitAnalyzer',
                             vertexSrc = cms.vstring('pixelVertices'),
                             trackSrc = cms.InputTag('pixelTracks'),
                             doTracking = cms.untracked.bool(False),
                             doMC = cms.untracked.bool(False),
                             L1gtReadout = cms.InputTag('gtDigis'),
                             )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('PixelTree-4.root')
                                   )


process.analyze = cms.Path(process.siPixelRecHits*process.ana)


