import FWCore.ParameterSet.Config as cms

process = cms.Process("PIXELTRACKVALIDATION")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Validation.RecoHI.TrackValidationHeavyIons_cff')

process.GlobalTag.globaltag = 'MC_38Y_V12::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/himc/Fall10/Hydjet_Quenched_MinBias_2760GeV/GEN-SIM-RECO/MC_38Y_V12-v1/0003/BED16369-95D6-DF11-B2D5-001EC9AAD5DA.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('validation.root')
)

process.Timing = cms.Service("Timing")


# reconstruct pixel tracks
process.hiLowPtPixelTracks = process.hiPixel3PrimTracks.clone()
process.hiLowPtPixelTracks.RegionFactoryPSet.RegionPSet.ptMin = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.ptMin = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.nSigmaTipMaxTolerance = cms.double(4.0)
process.hiLowPtPixelTracks.FilterPSet.tipMax = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.nSigmaLipMaxTolerance = cms.double(4.0)
process.hiLowPtPixelTracks.FilterPSet.lipmax = cms.double(0.2)
process.hiLowPtPixelTracks.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000


#customize validation for pixel tracks
process.findableSimTracks.ptMin = 0.2
process.findableSimTracks.minHit = 3
process.hiTrackValidator.label = ['recoTrackSelector']
process.hiTrackValidator.histoProducerAlgoBlock.minPt = 0.0
process.hiTrackValidator.histoProducerAlgoBlock.maxPt = 3.0
process.hiTrackValidator.ptMinTP = 0.2
process.hiTrackValidator.histoProducerAlgoBlock.nintPt = 60
process.hiTrackValidator.histoProducerAlgoBlock.useLogPt = False
process.hiTrackValidator.label_tp_effic = cms.InputTag("mergedtruth","MergedTrackTruth")
process.hiTrackValidator.label_tp_fake = cms.InputTag("mergedtruth","MergedTrackTruth")
process.hiTrackValidator.skipHistoFit = False
process.hiTrackValidator.outputFile = 'validation.root'

process.recoTrackSelector = cms.EDFilter("RecoTrackSelector",
    src = cms.InputTag("hiLowPtPixelTracks"),
    maxChi2 = cms.double(100000.0),
    tip = cms.double(9999.0),
    minRapidity = cms.double(-2.5),
    lip = cms.double(9999.0),
    ptMin = cms.double(-9999.0),
    maxRapidity = cms.double(2.5),
    quality = cms.vstring(),
    algorithm = cms.vstring(),
    minHit = cms.int32(0),
    min3DHit = cms.int32(0),
    beamSpot = cms.InputTag("offlineBeamSpot")
)


# change pseudorapidity range
process.hiTrackValidator.minRapidityTP = -1.0
process.hiTrackValidator.maxRapidityTP = 1.0
process.recoTrackSelector.minRapidity = -1.0
process.recoTrackSelector.maxRapidity = 1.0

process.rechit_step = cms.Path(process.siPixelRecHits)
process.pixeltrack_step = cms.Path(process.hiLowPtPixelTracks)
process.validation_step = cms.Path(process.recoTrackSelector * process.hiTrackValidation)

process.schedule = cms.Schedule(process.rechit_step, process.pixeltrack_step, process.validation_step)
