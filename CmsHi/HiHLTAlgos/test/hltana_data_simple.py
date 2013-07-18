import FWCore.ParameterSet.Config as cms

process = cms.Process('HltAna')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/d101/frankma/data/HIAllPhysics/ZS-v2/A2AD1439-F242-E011-A425-003048F00412.root'
    ))

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(1000))

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_39X_V6B::All'

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFhits"),
	nonDefaultGlauberModel = cms.string(""),
	centralitySrc = cms.InputTag("hiCentrality")
	)

# tree output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("hiOHTrees.root"))

# Define the analyzer modules
process.load("HLTrigger.HLTanalyzers.HI_HLTAnalyser_cff")
process.hltanalysis.RunParameters.Debug = False
process.hltanalysis.RunParameters.UseTFileService = True
process.hltanalysis.RunParameters.Monte = False
process.hltanalysis.RunParameters.DoMC = False
process.hltanalysis.RunParameters.DoJets = True
process.hltanalysis.RunParameters.DoPhotons = True
process.hltanalysis.RunParameters.DoSuperClusters = True
process.hltanalysis.recjets  = "iterativeConePu5CaloJets"
process.hltanalysis.BarrelPhoton = "correctedIslandBarrelSuperClusters"
process.hltanalysis.EndcapPhoton = "correctedIslandEndcapSuperClusters"
process.hltanalysis.l1GtReadoutRecord = cms.InputTag("gtDigis")


# =============== Final Paths =====================
process.ana_step = cms.Path( process.centralityBin * process.hltanalysis )
