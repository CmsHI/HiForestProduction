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
    input = cms.untracked.int32(3)
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
'file:/scratch/appelte1/HIJING_mix_Dijet120_1000ev.root'
    )
)

# JETS

process.load('CmsHi.JetAnalysis.ExtraPfReco_cff')
process.load('CmsHi.JetAnalysis.ExtraJetReco_cff')
process.load('CmsHi.JetAnalysis.PatAna_cff')
process.PFTowers.src = cms.InputTag("particleFlow")

#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltSingleTrigger = process.hltHighLevel.clone()
#process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.GlobalTag.globaltag = 'STARTHI53_V25::All'

#process.pPbTrkCorr.occByCentrality = cms.bool(False)
#process.pPbTrkCorr.occByJetEt = cms.bool(True)
#process.pPbTrkCorr.occBins = cms.vdouble( 0., 40.,80.,120.,160.,200.,300.,500.)

process.trkCorr_akPu3PF_Et80dR03 = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3PFpatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.3),
    invertJetCuts = cms.bool(False)
)

process.trkCorr_akPu3PF_Et80dR08 = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3PFpatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.8),
    invertJetCuts = cms.bool(False)
)

process.trkCorr_akPu3PF_Et80dR03inv = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3PFpatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.3),
    invertJetCuts = cms.bool(True)
)

process.trkCorr_akPu3PF_Et80dR08inv = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3PFpatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.8),
    invertJetCuts = cms.bool(True)
)

process.trkCorr_akPu3Calo_Et80dR03 = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3CalopatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.3),
    invertJetCuts = cms.bool(False)
)

process.trkCorr_akPu3Calo_Et80dR08 = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3CalopatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.8),
    invertJetCuts = cms.bool(False)
)

process.trkCorr_akPu3Calo_Et80dR03inv = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3CalopatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.3),
    invertJetCuts = cms.bool(True)
)

process.trkCorr_akPu3Calo_Et80dR08inv = process.pPbTrkCorr.clone(
    jetSrc = cms.InputTag('akPu3CalopatJets'),
    applyJetCuts = cms.bool(True),
    jetRadius = cms.double(0.8),
    invertJetCuts = cms.bool(True)
)



process.jetReco= cms.Sequence(
    process.PFTowers *
    process.akPu3PFJets*
    process.akPu3PFcorr*
    process.akPu3PFpatJets*
    process.akPu3CaloJets*
    process.akPu3CaloJetID*
    process.akPu3Calocorr*
    process.akPu3CalopatJets
)

process.p = cms.Path( 
    process.siPixelRecHits *
    process.pACentrality *
    process.jetReco *
    process.tpRecoAssocGeneralTracks *
    process.trkCorr_akPu3PF_Et80dR03 * 
    process.trkCorr_akPu3PF_Et80dR08 * 
    process.trkCorr_akPu3PF_Et80dR03inv * 
    process.trkCorr_akPu3PF_Et80dR08inv * 
    process.trkCorr_akPu3Calo_Et80dR03 * 
    process.trkCorr_akPu3Calo_Et80dR08 * 
    process.trkCorr_akPu3Calo_Et80dR03inv * 
    process.trkCorr_akPu3Calo_Et80dR08inv  
)

