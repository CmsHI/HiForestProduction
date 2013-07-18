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
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality")
#  pPbRunFlip = cms.uint32(99999999)
  )

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames =  cms.untracked.vstring(
# private reco at MIT
'/store/user/kjung/VZhukova/EposBatch13/pPbRecoShift_f1521_e.root'
    )
)

#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltSingleTrigger = process.hltHighLevel.clone()
#process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.GlobalTag.globaltag = 'STARTHI53_V17::All'

process.p = cms.Path( process.PAcollisionEventSelection * 
                      process.pACentrality *
                      process.trkAnaMinBias 
)

