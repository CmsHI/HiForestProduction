import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('initialEvent',mult=ivars.multiplicity.singleton,info="for testing")

ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/HICorePhysics/Jet50-PromptReco-Runs_151217_RECOPAT-v0/9dff160a979f15b232d23f89936b0238/RECOPAT_10_1_MY1.root"

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
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration/StandardSequences/ReconstructionHeavyIons_cff')
#process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_39Y_V4HI::All'
process.GlobalTag.globaltag = 'GR10_P_V12::All'  #39X
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



# Ecal spike filter only
# put a jet filter here?

process.load("MNguyen.Configuration.HI_JetSkim_cff")
process.jetSkimSequence.remove(process.collisionEventSelection)

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")


process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )


process.inclusiveJetAnalyzer = cms.EDAnalyzer("InclusiveJetAnalyzer",
                                              jetTag = cms.InputTag("icPu5patJets"),
                                              isMC = cms.untracked.bool(False), 
                                              useCentrality = cms.untracked.bool(True),
                                              L1gtReadout = cms.InputTag("gtDigis"),
                                              hltTrgResults = cms.untracked.string("TriggerResults::HLT"),
                                              hltTrgNames  = cms.untracked.vstring('HLT_HIMinBiasHfOrBSC_Core',
                                                                                   'HLT_HIJet35U',
                                                                                   'HLT_HIJet35U_Core',
                                                                                   'HLT_HIJet50U_Core')
                                              )


process.ic3JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ic3JetAnalyzer.jetTag = 'ic3patJets'

process.ic4JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ic4JetAnalyzer.jetTag = 'ic4patJets'

process.ic5JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ic5JetAnalyzer.jetTag = 'ic5patJets'


process.ak3JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak3JetAnalyzer.jetTag = 'ak3patJets'

process.ak4JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak4JetAnalyzer.jetTag = 'ak4patJets'

process.ak5JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak5JetAnalyzer.jetTag = 'ak5patJets'

process.ak7JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak7JetAnalyzer.jetTag = 'ak7patJets'

process.kt4JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.kt4JetAnalyzer.jetTag = 'kt4patJets'

process.ic3PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ic3PFJetAnalyzer.jetTag = 'ic3PFpatJets'

process.ic4PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ic4PFJetAnalyzer.jetTag = 'ic4PFpatJets'

process.ic5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ic5PFJetAnalyzer.jetTag = 'ic5PFpatJets'


process.ak3PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak3PFJetAnalyzer.jetTag = 'ak3PFpatJets'

process.ak4PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak4PFJetAnalyzer.jetTag = 'ak4PFpatJets'

process.ak5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak5PFJetAnalyzer.jetTag = 'ak5PFpatJets'

process.ak7PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak7PFJetAnalyzer.jetTag = 'ak7PFpatJets'



# ntuple output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string(ivars.output))
# Paths

process.jetSkimPath = cms.Path(process.jetSkimSequence)

process.p = cms.Path(    process.jetSkimSequence
                         *process.inclusiveJetAnalyzer
                         #*process.ic3JetAnalyzer
                         #*process.ic4JetAnalyzer
                         *process.ic5JetAnalyzer
                         #*process.ak3JetAnalyzer
                         #*process.ak4JetAnalyzer
                         *process.ak5JetAnalyzer
                         *process.kt4JetAnalyzer
                         #*process.ak7JetAnalyzer
                         #*process.ic3PFJetAnalyzer
                         #*process.ic4PFJetAnalyzer
                         #*process.ic5PFJetAnalyzer
                         #*process.ak3PFJetAnalyzer
                         #*process.ak4PFJetAnalyzer
                         #*process.ak5PFJetAnalyzer
                         #*process.ak7PFJetAnalyzer
                         )

## Schedule
process.schedule = cms.Schedule(
    process.p
    )


