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
                            fileNames = cms.untracked.vstring("/store/hidata/HIRun2013/PAHighPt/RAW/v1/000/210/498/00000/ACC49F8C-1B63-E211-9711-001D09F24DA8.root")
#                            fileNames = cms.untracked.vstring("file:PAPhysics_RAWRECO_inRECO_9_1_8mR.root")
			    )

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(500))


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
process.load('Configuration.StandardSequences.RawToDigi_cff')

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

process.load('L1Trigger.Configuration.L1Extra_cff')
process.load('CmsHi.HiHLTAlgos.hltanalysis_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.pACentrality.produceHFtowers = cms.bool(False)
process.pACentrality.produceEcalhits = cms.bool(False)
process.pACentrality.produceBasicClusters = cms.bool(False)
process.pACentrality.produceZDChits = cms.bool(False)
process.pACentrality.produceETmidRapidity = cms.bool(False)
process.pACentrality.produceTracks = cms.bool(False)
process.pACentrality.producePixelTracks = cms.bool(False)
process.pACentrality.reUseCentrality = cms.bool(False)

# to be updated
process.pACentrality.producePixelhits = cms.bool(False)

process.load("CmsHi/HiHLTAlgos.hievtanalyzer_cfi")
process.hiEvtAnalyzer.doMC = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.hiEvtAnalyzer.doVertex = cms.bool(False)
process.hiEvtAnalyzer.doMC = cms.bool(False)

process.hltanalysis.dummyBranches = cms.untracked.vstring()

process.reco_extra=cms.Path(process.RawToDigi*process.L1Extra*process.siPixelClusters*process.siPixelRecHits*process.hfreco*process.pACentrality)

process.hltAna = cms.EndPath(process.hltanalysis+process.hiEvtAnalyzer)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.reco_extra,
    process.hltAna
    )
