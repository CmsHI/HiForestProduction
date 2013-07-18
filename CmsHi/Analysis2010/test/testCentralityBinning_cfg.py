
import FWCore.ParameterSet.Config as cms

process = cms.Process("TESTC")

process.HeavyIonGlobalParameters = cms.PSet(
    centralitySrc = cms.InputTag("hiCentrality"),
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("")
    )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_39X_V1::All'

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
"rfio:/castor/cern.ch/cms/store/caf/user/frankma/HR10Exp3/r150308HFSkim/skim_RECO_9_1_2w4.root",

)
    
                            )

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.analyze = cms.EDAnalyzer("AnalyzerWithCentrality")

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string("histogramsAndTable.root")
                                   )

process.p = cms.Path(process.centralityBin*process.analyze)

