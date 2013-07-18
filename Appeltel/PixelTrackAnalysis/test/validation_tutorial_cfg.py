import FWCore.ParameterSet.Config as cms

process = cms.Process("PIXELTRACKVALIDATION")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Validation.RecoHI.TrackValidationHeavyIons_cff')

process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")

process.GlobalTag.globaltag = 'START39_V4HI::All'


from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    centralitySrc = cms.InputTag("hiCentrality")
    )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 'file:association.root'),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('validation3D.root')
)

process.Timing = cms.Service("Timing")


process.hitrkPixelEffAnalyzer = cms.EDAnalyzer('HiPixelTrkEffAnalyzer',
    tracks = cms.untracked.InputTag('hiMergedPtSplit12Tracks'),
    label_tp_effic = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),
    label_tp_fake = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),                 
    associatormap = cms.untracked.InputTag('tpRecoAssocMerged12Tracks'),
    vertices = cms.untracked.InputTag('hiSelectedVertex',''),
    beamspot = cms.untracked.InputTag('offlineBeamSpot'),                         
    fillHistograms = cms.bool(True),
    fillNtuples = cms.bool(False),
    constPtBins = cms.bool(False),
    lowPtMode = cms.bool(False),
    doAssociation = cms.untracked.bool(False),
    hasSimInfo = cms.untracked.bool(True)
)


process.validation_step = cms.Path(process.hitrkPixelEffAnalyzer)

process.schedule = cms.Schedule( process.validation_step)
