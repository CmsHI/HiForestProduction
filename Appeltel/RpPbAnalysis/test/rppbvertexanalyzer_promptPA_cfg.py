import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cff')
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('vertexAnalysis.root')
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
#'/store/express/HIRun2013/ExpressPhysics/FEVT/Express-v1/000/210/498/00000/D8C0777B-1863-E211-9265-BCAEC5329719.root'
'/store/hidata/HIRun2013/PAMinBiasUPC/RECO/PromptReco-v1/000/210/759/00000/54822F4C-8468-E211-ADF9-001D09F24303.root'
  )
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrackTrigger = process.hltHighLevel.clone()
process.hltSingleTrackTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]
#process.hltSingleTrigger.HLTPaths = ["HLT_PAMinBiasHfOrBHC_v1"]
#process.hltSingleTrackTrigger.HLTPaths = ["HLT_PAPixelTracks_Multiplicity190_v1"]

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

process.vtxAna  = cms.EDAnalyzer('RpPbVertexAnalyzer',
   vertexSrc = cms.InputTag("offlinePrimaryVertices")
)
process.vtxAnaBS  = cms.EDAnalyzer('RpPbVertexAnalyzer',
   vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS")
)

process.vtxAnaCut = process.vtxAna.clone()
process.vtxAnaBSCut = process.vtxAnaBS.clone()
process.vtxAnaCutG = process.vtxAna.clone()
process.vtxAnaCutGloose = process.vtxAna.clone()
process.vtxAnaCutGtight = process.vtxAna.clone()
process.vtxAnaCutGplus = process.vtxAna.clone()
process.vtxAnaCutGplusNV = process.vtxAna.clone()
process.vtxAnaCutE = process.vtxAna.clone()
process.vtxAnaCutEandG = process.vtxAna.clone()
process.vtxAnaCutW = process.vtxAna.clone()
process.vtxAnaCutWplus = process.vtxAna.clone()
process.vtxAnaCutGplusplus = process.vtxAna.clone()


# pilot run
#process.GlobalTag.globaltag = 'GR_P_V41_AN2::All'

# 2013 prompt
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'

process.minBias = cms.Path( process.hltSingleTrackTrigger * 
                      process.siPixelRecHits *
#                    ~ process.hltScoutingTrigger *
#                      process.PAcollisionEventSelection *
                      process.pACentrality * 
                      process.vtxAna *
                      process.vtxAnaBS
)

process.minBiasCutG = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutG *
                           process.pACentrality *
                           process.vtxAnaCutG 
)
process.minBiasCutGloose = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#			   process.PAcollisionEventSelection *
			   process.pileupVertexFilterCutGloose *
			   process.pACentrality *
			   process.vtxAnaCutGloose 
)

process.minBiasCutGtight = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutGtight *
                           process.pACentrality *
                           process.vtxAnaCutGtight
)

process.minBiasCutE = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutE *
                           process.pACentrality *
                           process.vtxAnaCutE
)

process.minBiasCutEandG = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutEandG *
                           process.pACentrality *
                           process.vtxAnaCutEandG
)

process.minBiasCutGplus = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.vtxAnaCutGplus
)

process.minBiasCutGplusNV = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutGplusNV *
                           process.pACentrality *
                           process.vtxAnaCutGplusNV
)

process.minBiasCutWplus = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutWplus *
                           process.pACentrality *
                           process.vtxAnaCutWplus
)

process.minBiasCutW = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutW *
                           process.pACentrality *
                           process.vtxAnaCutW
)

process.minBiasCutGplusplus = cms.Path( process.hltSingleTrackTrigger *
                      process.siPixelRecHits *
#                           process.PAcollisionEventSelection *
                           process.pileupVertexFilterCutGplusplus *
                           process.pACentrality *
                           process.vtxAnaCutGplusplus
)



process.schedule = cms.Schedule(
    process.minBias, process.minBiasCutG, process.minBiasCutGloose,
    process.minBiasCutGtight, process.minBiasCutE, process.minBiasCutEandG,
    process.minBiasCutGplus, process.minBiasCutGplusNV,
    process.minBiasCutW, process.minBiasCutWplus,
    process.minBiasCutGplusplus 
) 

#process.schedule = cms.Schedule(process.minBias)

