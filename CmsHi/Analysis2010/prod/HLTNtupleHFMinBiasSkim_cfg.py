    ##################################################################
                            
import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Geometry_cff')
#process.load('HLTrigger.Configuration.HLT_HIon_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR10_P_V12::All"

# =============== Conditions, Centrality ===============

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )
    
# =============== input file setting =====================

# =============== 2.36 TeV MC Sample =====================

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/cms/store/caf/user/frankma/HR10Exp3/r150308HFSkim/skim_RECO_9_1_2w4.root")
                            )

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# =============== Trigger selection ====================
process.load("FWCore.Modules.preScaler_cfi")
process.preScaler.prescaleFactor = 30

# HF coincidence
process.load("L1Trigger.Skimmer.l1Filter_cfi")
process.L1HfOrBscCoinc = process.l1Filter.clone(
      algorithms = cms.vstring("L1_BscMinBiasInnerThreshold1","L1_HcalHfCoincidencePm")
      )

# HI MB Event Selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

# =============== Analysis modules =====================

# ------------------- OpenHLT --------------------------

# Define the analyzer modules
process.load("HLTrigger.HLTanalyzers.HI_HLTAnalyser_cff")
process.hltanalysis.l1GtReadoutRecord = cms.InputTag("gtDigis")
process.hltanalysis.RunParameters = cms.PSet(
    HistogramFile  = cms.untracked.string("dummy"),
    UseTFileService  = cms.untracked.bool(True),
    Monte                = cms.bool(False),
    Debug                = cms.bool(False),
    
    ### added in 2010 ###
    DoHeavyIon           = cms.untracked.bool(True),
    DoMC           = cms.untracked.bool(False),
    DoAlCa           = cms.untracked.bool(False),
    DoTracks           = cms.untracked.bool(False),
    DoVertex           = cms.untracked.bool(True),
    DoJets           = cms.untracked.bool(True),
    
    ### MCTruth
    DoParticles          = cms.untracked.bool(False),
    DoRapidity           = cms.untracked.bool(False),
    DoVerticesByParticle = cms.untracked.bool(False),
    
    ### Egamma
    DoPhotons            = cms.untracked.bool(True),
    DoElectrons          = cms.untracked.bool(False),
    DoSuperClusters      = cms.untracked.bool(True),
    
    ### Muon
    DoMuons            = cms.untracked.bool(True),
    DoL1Muons            = cms.untracked.bool(False),
    DoL2Muons            = cms.untracked.bool(False),
    DoL3Muons            = cms.untracked.bool(False),
    DoOfflineMuons       = cms.untracked.bool(False),
    DoQuarkonia          = cms.untracked.bool(False)
    )

# rec objects
process.hltanalysis.recjets  = "iterativeConePu5CaloJets"
process.hltanalysis.BarrelPhoton = "correctedIslandBarrelSuperClusters"
process.hltanalysis.EndcapPhoton = "correctedIslandEndcapSuperClusters"

#========================= Photons ====================================

import HLTrigger.HLTfilters.hltHighLevel_cfi
# Minimum bias trigger selection (Runs after 150882)
process.HIMinBiasHFOrBsc = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.HIMinBiasHFOrBsc.HLTPaths = ["HLT_HIMinBiasHfOrBSC"]

process.HIphoton15 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.HIphoton15.HLTPaths = ["HLT_HIPhoton15"]

process.photon_filter_step = cms.Path(process.HIphoton15 * process.L1HfOrBscCoinc*process.collisionEventSelection )  #tilde means NOT

#===============================================================


# Centrality Objects
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = True
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceBasicClusters = False
process.hiCentrality.produceZDChits = True
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.doPixelCut = cms.bool(True)

process.rechits = cms.EDAnalyzer("RecHitTreeProducer",
                                 doFastJet = cms.untracked.bool(False),
                                 doEbyEonly = cms.untracked.bool(True)
                                 )

# =============== Output ================================

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('OpenHLT.root')
                                   )

process.HLToutput = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('hlt_filter_step')),
    fileName = cms.untracked.string("HLT_HIMinBiasHfOrBSC.root")
)

process.photonOutput = cms.OutputModule("PoolOutputModule",
        splitLevel = cms.untracked.int32(0),
        outputCommands = process.RECOEventContent.outputCommands,
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('photon_filter_step')),
        fileName = cms.untracked.string('HiPhoton15_L1HfOrBsc.root'),
        dataset = cms.untracked.PSet(
            filterName = cms.untracked.string(''),
            dataTier = cms.untracked.string('GEN-SIM-RECO')
            )
        )

process.ana_step = cms.Path(process.hiCentrality *
                            process.centralityBin *
                            process.rechits *
                            process.hltanalysis
                            )

process.hlt_filter_step = cms.Path(process.HIMinBiasHFOrBsc*process.collisionEventSelection*process.preScaler)
process.output_step = cms.EndPath(process.HLToutput*process.photonOutput)


