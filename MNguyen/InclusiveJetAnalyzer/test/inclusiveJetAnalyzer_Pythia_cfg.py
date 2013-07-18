import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('initialEvent',mult=ivars.multiplicity.singleton,info="for testing")

ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/PatFilesWithPF_HIL1L2L3/Summer10_Pythia_Pt80/Pat_0.root"

ivars.output = 'jetTree.root'
ivars.maxEvents = -1
ivars.initialEvent = 1

ivars.parseArguments()

import FWCore.ParameterSet.Config as cms
process = cms.Process("ANALYSIS")

# Services
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration/StandardSequences/ReconstructionHeavyIons_cff')
#process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_39Y_V4HI::All'
process.GlobalTag.globaltag = 'START39_V6HI::All'  #39X
process.MessageLogger.cerr.FwkReport.reportEvery=1
#process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

#process.Timing = cms.Service("Timing")

#Input source

#process.load("MNguyen.InclusiveJetAnalyzer.Sources.source_Pyquen_Unquenched_cff")

process.source = cms.Source (
    "PoolSource",    
    fileNames = cms.untracked.vstring(ivars.files),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
)

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(ivars.maxEvents)
)

# put a jet filter here?

# Ecal spike filter

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("Hydjet_2760GeV"),
    centralitySrc = cms.InputTag("hiCentrality")
    )


process.inclusiveJetAnalyzer = cms.EDAnalyzer("InclusiveJetAnalyzer",
                                              jetTag = cms.InputTag("icPu5patJets"),
                                              isMC = cms.untracked.bool(True), 
                                              useCentrality = cms.untracked.bool(False)
                                              )


# ntuple output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string(ivars.output))
# Paths

#process.filter = cms.Path(process.filter_seq)
process.p = cms.Path(#process.filter_seq*
    process.inclusiveJetAnalyzer)

# Schedule
process.schedule = cms.Schedule(#process.filter,
                                process.p
                                )


