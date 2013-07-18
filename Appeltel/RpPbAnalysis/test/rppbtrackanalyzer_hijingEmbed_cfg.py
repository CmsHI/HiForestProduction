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

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/user/vzhukova/HIJING_GEN-SIM_NEWFIX_batch9/Dijet120_mix_RECO_v/b7d33bba7673cdb1ee6f4983c0800c79/Pyquen_Dijet50_pp_Unquenched_TuneZ2_5p02TeV_RECO_10_1_vZh.root'
  )
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrackTrigger = process.hltHighLevel.clone()
process.hltSingleTrackTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]
#process.hltSingleTrigger.HLTPaths = ["HLT_PAMinBiasHfOrBHC_v1"]

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltScoutingTrigger = process.hltHighLevel.clone()
process.hltScoutingTrigger.HLTPaths = ["HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v1",
                                       "HLT_PAPixelTrackMultiplicity140_Jet80_NoJetID_v1"
                                      ]
process.hltScoutingTrigger.andOr = cms.bool(True)

process.hltJet80Trigger = process.hltHighLevel.clone()
process.hltJet80Trigger.HLTPaths = ["HLT_PAJet80_NoJetID_v1"]
#process.hltSingleTrigger.HLTPaths = ["HLT_PAMinBiasHfOrBHC_v1"]

process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.GlobalTag.globaltag = 'STARTHI53_V17::All'

process.minBias = cms.Path( process.hltSingleTrackTrigger * 
#                    ~ process.hltScoutingTrigger *
                      process.PAcollisionEventSelection *
#                      process.siPixelRecHits *
                      process.pACentrality * 
                      process.trkAnaMinBias
)

process.highPt = cms.Path( process.hltJet80Trigger *
                           process.PAcollisionEventSelection *
                           process.pACentrality *
                           process.trkAna_etaFull_highPtQuality *
                           process.trkAna_etaFull_highPt
)

process.schedule = cms.Schedule(
    process.minBias, process.highPt
) 
