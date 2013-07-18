import FWCore.ParameterSet.Config as cms

process = cms.Process('CORRFLOWSKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")

process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
process.load('Appeltel.PixelTracksRun2010.HiTrackCandidates_cff')
process.load('Appeltel.PixelTracksRun2010.HICorrFlowSkimEventContent_cff')


process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")
process.load("RecoHI.HiTracking.LowPtTracking_PbPb_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30)
)



# And the logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring('Unknown', 
        'ProductNotFound', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ModuleFailure', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError', 
        'FatalRootError', 
        'NotFound')
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/hidata/HIRun2010/HIAllPhysics/RECO/ZS-v2/0032/00E5FF3F-A143-E011-8C93-003048CF96DA.root'
  )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.HICorrFlowEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('filter_step')
                ),
    fileName = cms.untracked.string( 'file:flowskim.root')
)

process.output.outputCommands += ['keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*']
process.output.outputCommands += ['keep *_dedxPTTru40_*_*']
process.output.outputCommands += ['keep *_hiGoodTightTracks_*_*']
process.output.outputCommands += ['keep *_hiGoodTightMergedTracks_*_*']
# Jet Stuff

process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.icPu5CaloJetsL2L3 = cms.EDProducer('CaloJetCorrectionProducer',
    src = cms.InputTag('iterativeConePu5CaloJets'),
                                   correctors = cms.vstring('ic5CaloL2L3')
                                   )
# ECAL/HCAL Noise Cleaning

#from RecoHI.HiEgammaAlgos.hiEcalSpikeFilter_cfi import *
process.load("CmsHi.PhotonAnalysis.hiEcalRecHitSpikeFilter_cfi")
process.hiEcalRecHitSpikeFilter.minEt = 20.0

#HCAL cleaning
process.load("JetMETAnalysis.HcalReflagging.hbherechitreflaggerJETMET_cfi")

process.hbheReflagNewTimeEnv = process.hbherechitreflaggerJETMET.clone()
process.hbheReflagNewTimeEnv.timingshapedcutsParameters.hbheTimingFlagBit=cms.untracked.int32(8)

# HCAL Timing
process.hcalTimingFilter = cms.EDFilter("HcalTimingFilter",
                                        hbheHits = cms.untracked.InputTag("hbheReflagNewTimeEnv")
                                )


# hcal noise filter
process.load("CommonTools.RecoAlgos.HBHENoiseFilter_cfi")
process.HBHENoiseFilter.minRatio = cms.double(-99999.0)
process.HBHENoiseFilter.maxRatio = cms.double(99999.0)
process.HBHENoiseFilter.minZeros = cms.int32(100)

# dE/dx estimators

process.dedxPTTru40 = cms.EDProducer("DeDxEstimatorProducerPixelTripplet",
    tracks                     = cms.InputTag("hiConformalPixelTracks"),
    trajectoryTrackAssociation = cms.InputTag("hiConformalPixelTracks"),

    estimator      = cms.string('truncated'),
    fraction       = cms.double(0.4),

    UseStrip       = cms.bool(False),
    UsePixel       = cms.bool(True),
    MeVperADCStrip = cms.double(3.61e-06*250),
    MeVperADCPixel = cms.double(3.61e-06),

    MisCalib_Mean      = cms.untracked.double(1.0),
    MisCalib_Sigma     = cms.untracked.double(0.00),

    UseCalibration  = cms.bool(False),
    calibrationPath = cms.string(""),
    ShapeTest       = cms.bool(False),
)

#Trigger Selection

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMinBiasHFOrBSC = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHf","HLT_HIMinBiasBSC", "HLT_HIMinBiasHfOrBSC" ]
#process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHf_Core","HLT_HIMinBiasBSC_Core", "HLT_HIMinBiasHfOrBSC_Core" ]
process.hltMinBiasHFOrBSC.andOr = cms.bool(True) # this is the default meaning either of the paths above
process.hltMinBiasHFOrBSC.throw = cms.bool(False) # don't throw exception since some runs have only one trigger or the other

# GoodTight Tracks

process.hiGoodTightMergedTracks = process.hiGoodMergedTracks.clone(
        TrackProducer1 = 'hiGoodTightTracks'
)

# Other statements
process.GlobalTag.globaltag = 'GR_R_39X_V6B::All'
#process.GlobalTag.globaltag = 'GR10_P_V12::All'

# Path and EndPath definitions

process.reco_step = cms.Path( process.hltMinBiasHFOrBSC * 
                   process.collisionEventSelection *
    process.hiEcalRecHitSpikeFilter *
    process.hbheReflagNewTimeEnv *
    process.hcalTimingFilter *
    process.HBHENoiseFilter *
                   process.hiGoodTracksSelection *
                   process.hiPostGlobalPrimTracksSelection *
                   process.hiGoodTightTracksSelection *
                   process.conformalPixelTrackReco *
#                   process.hiLowPtPixelTracks *
#                   process.hiMergedPtSplit12Tracks *
                   process.hiGoodMergedTracks *
                   process.hiGoodTightMergedTracks *
    process.icPu5CaloJetsL2L3 *
    process.dedxPTTru40
)

process.filter_step = cms.Path( process.hltMinBiasHFOrBSC * 
    process.collisionEventSelection *
    process.hiEcalRecHitSpikeFilter *
    process.hbheReflagNewTimeEnv *
    process.hcalTimingFilter *
    process.HBHENoiseFilter
)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(
    process.filter_step,
    process.reco_step,
    process.endjob_step,
    process.out_step
)


