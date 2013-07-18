import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trackAnalysis.root')
)

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowers"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("hiCentrality")
#  pPbRunFlip = cms.uint32(99999999)
  )

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'file:/path/to/myfile.root'
  )
)


process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.GlobalTag.globaltag = 'STARTHI44_V12::All'

process.trkAna_quality.trackSrc = cms.InputTag("hiGeneralTracks")
process.trkAna_quality.vertexSrc = cms.InputTag("hiSelectedVertex")
process.trkAna_quality.dxyErrMax = cms.double(999.0)
process.trkAna_quality.dzErrMax = cms.double(999.0)
process.trkAna_quality.ptErrMax = cms.double(999.0)

process.minBias = cms.Path( 
                      process.PAcollisionEventSelection *
                      process.siPixelRecHits *
                      process.pACentrality * 
                      process.trkAna_quality
)


process.schedule = cms.Schedule(
    process.minBias
) 
