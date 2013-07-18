import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Appeltel.RpPbAnalysis.RpPbTrackingCorrections_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trackCorrections.root')
)

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality")
#  pPbRunFlip = cms.uint32(99999999)
  )

process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpRecoAssocGeneralTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("generalTracks")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames =  cms.untracked.vstring(
'/store/user/vzhukova/HIJING_GEN-SIM_YUE-SHI_Minbias_2_v1/HIJING_RECO_YUE-SHI_Minbias__2_v1/b7d33bba7673cdb1ee6f4983c0800c79/hijing_reco_fix_9_2_wnY.root'
    )
)

#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltSingleTrigger = process.hltHighLevel.clone()
#process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.GlobalTag.globaltag = 'STARTHI53_V25::All'

process.p = cms.Path( 
                      process.siPixelRecHits *
                      process.pACentrality *
                      process.tpRecoAssocGeneralTracks *
                      process.pPbTrkCorr 
)

