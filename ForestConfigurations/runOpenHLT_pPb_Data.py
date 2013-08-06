import FWCore.ParameterSet.VarParsing as VarParsing

import FWCore.ParameterSet.Config as cms

process = cms.Process('hiForestAna2011')

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)


hiTrackQuality = "highPurity"              # iterative tracks
#hiTrackQuality = "highPuritySetWithPV"    # calo-matched tracks

doElectrons = False
doRegitForBjets = False


#####################################################################################

process.load("CmsHi.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V2 for pPb",
                                          "PF : generalTracks (pp) with highPurity",
                                          "EP Flattening OFF",
                                          "Electrons OFF",
                                          "Preshower OFF"
                                          )

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
#                            fileNames = cms.untracked.vstring("/store/group/phys_heavyions/icali/PAPhysics/pAPilotRun_Run202792GoodLumis_RAWRECO_L1Em_PrescaleActiveBitsSkimNoZB_CMSSW528_V94_FinalWorkflow_2MHz_v2_v1_v2/f3394926c5028783289fd2cd57b36909/PAPhysics_RAWRECO_inRECO_9_1_8mR.root")
                            fileNames = cms.untracked.vstring("file:FE6E7D35-895E-E211-B3E4-003048D2BEA8.root")
#                            fileNames = cms.untracked.vstring("file:PAPhysics_RAWRECO_inRECO_9_1_8mR.root")
			    )

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1))


#####################################################################################
# Load some general stuff
#####################################################################################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

process.GlobalTag.globaltag = 'GR_E_V33::All'


# load centrality

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality")
)
        
#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("pPb_openHLT.root"))

# Define Analysis sequencues
process.load('CmsHi.JetAnalysis.EventSelection_cff')

# Event tree
process.load("CmsHi/HiHLTAlgos.hievtanalyzer_cfi")
process.hiEvtAnalyzer.doMC = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(True)
process.hiEvtAnalyzer.Centrality = cms.InputTag("pACentrality") 
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVerticesWithBS")
process.hiEvtAnalyzer.doEvtPlane = False

process.primaryVertexFilter.src = cms.InputTag("offlinePrimaryVerticesWithBS")

process.reco_extra =  cms.Path(
    process.siPixelRecHits *
    process.pACentrality
)    
    
process.ana_step          = cms.Path( process.hiEvtAnalyzer +
                                      process.HiForest
)


# Event selection filters ########################################################################################
process.hfPosFilter2 = cms.EDFilter("CandCountFilter",
                                      src = cms.InputTag("hfPosTowers"),
                                      minNumber = cms.uint32(2)
                                      )

process.hfNegFilter2 = cms.EDFilter("CandCountFilter",
                                    src = cms.InputTag("hfNegTowers"),
                                      minNumber = cms.uint32(2)
                                    )

process.pcollisionEventSelection = cms.Path(process.collisionEventSelection)
process.pPAcollisionEventSelectionPA = cms.Path(process.PAcollisionEventSelection)

process.pHBHENoiseFilter = cms.Path( process.HBHENoiseFilter )
process.phiEcalRecHitSpikeFilter = cms.Path(process.hiEcalRecHitSpikeFilter )

process.phfPosFilter3 = cms.Path(process.towersAboveThreshold+process.hfPosTowers+process.hfNegTowers+process.hfPosFilter3)
process.phfNegFilter3 = cms.Path(process.towersAboveThreshold+process.hfPosTowers+process.hfNegTowers+process.hfNegFilter3)
process.phfPosFilter2 = cms.Path(process.towersAboveThreshold+process.hfPosTowers+process.hfNegTowers+process.hfPosFilter2)
process.phfNegFilter2 = cms.Path(process.towersAboveThreshold+process.hfPosTowers+process.hfNegTowers+process.hfNegFilter2)
process.phfPosFilter1 = cms.Path(process.towersAboveThreshold+process.hfPosTowers+process.hfNegTowers+process.hfPosFilter)
process.phfNegFilter1 = cms.Path(process.towersAboveThreshold+process.hfPosTowers+process.hfNegTowers+process.hfNegFilter)
process.phltPixelClusterShapeFilter = cms.Path(process.siPixelRecHits+process.hltPixelClusterShapeFilter)
process.pprimaryvertexFilter = cms.Path(process.primaryVertexFilter)
 
process.load('L1Trigger.Configuration.L1Extra_cff')
process.load('CmsHi.HiHLTAlgos.hltanalysis_cff')

process.hltanalysis.dummyBranches = cms.untracked.vstring()

process.hltAna = cms.EndPath(process.hltanalysis)
process.reco_extra*=process.L1Extra

process.pAna = cms.EndPath(process.skimanalysis)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.reco_extra,
    process.ana_step,
    process.pcollisionEventSelection,process.pHBHENoiseFilter,process.phiEcalRecHitSpikeFilter,
    process.pPAcollisionEventSelectionPA,
    process.phfPosFilter3,process.phfNegFilter3,
    process.phfPosFilter2,process.phfNegFilter2,
    process.phfPosFilter1,process.phfNegFilter1,
    process.phltPixelClusterShapeFilter,process.pprimaryvertexFilter,
    process.hltAna,process.pAna
    )
