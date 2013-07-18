import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('initialEvent',mult=ivars.multiplicity.singleton,info="for testing")


ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/PatFilesWithPF_HIL1L2L3/Hydjet_Pyquen_UnquenchedDiJet_Pt80/Pat_FJ2Jets_0.root"

ivars.output = 'pfJetTree.root'
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


process.PFJetAnalyzer = cms.EDAnalyzer("PFJetAnalyzer",
                                       #jetTag = cms.InputTag("icPu5patJets"),
                                       jetTag = cms.InputTag("ak5PFpatJets"),
                                       recoJetTag = cms.InputTag("ak5PFJets"),
                                       pfCandidatesTag  = cms.InputTag("particleFlow",""),
                                       isMC = cms.untracked.bool(True), 
                                       useCentrality = cms.untracked.bool(True)
                                       )

process.ic3JetAnalyzer = process.PFJetAnalyzer.clone()
process.ic3JetAnalyzer.jetTag = 'ic3patJets'

process.ic4JetAnalyzer = process.PFJetAnalyzer.clone()
process.ic4JetAnalyzer.jetTag = 'ic4patJets'

process.ic5JetAnalyzer = process.PFJetAnalyzer.clone()
process.ic5JetAnalyzer.jetTag = 'ic5patJets'


process.ak3JetAnalyzer = process.PFJetAnalyzer.clone()
process.ak3JetAnalyzer.jetTag = 'ak3patJets'

process.ak4JetAnalyzer = process.PFJetAnalyzer.clone()
process.ak4JetAnalyzer.jetTag = 'ak4patJets'

process.ak5JetAnalyzer = process.PFJetAnalyzer.clone()
process.ak5JetAnalyzer.jetTag = 'ak5patJets'

process.ak7JetAnalyzer = process.PFJetAnalyzer.clone()
process.ak7JetAnalyzer.jetTag = 'ak7patJets'


process.ic3PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ic3PFJetAnalyzer.jetTag = 'ic3PFpatJets'

process.ic4PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ic4PFJetAnalyzer.jetTag = 'ic4PFpatJets'

process.ic5PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ic5PFJetAnalyzer.jetTag = 'ic5PFpatJets'


process.ak3PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ak3PFJetAnalyzer.jetTag = 'ak3PFpatJets'

process.ak4PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ak4PFJetAnalyzer.jetTag = 'ak4PFpatJets'

process.ak5PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ak5PFJetAnalyzer.jetTag = 'ak5PFpatJets'

process.ak7PFJetAnalyzer = process.PFJetAnalyzer.clone()
process.ak7PFJetAnalyzer.jetTag = 'ak7PFpatJets'


# ntuple output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string(ivars.output))
# Paths

#process.filter = cms.Path(process.filter_seq)
process.p = cms.Path(#process.filter_seq*
    #process.PFJetAnalyzer
    process.ic3PFJetAnalyzer
    *process.ic4PFJetAnalyzer
    *process.ic5PFJetAnalyzer
    *process.ak3PFJetAnalyzer
    *process.ak4PFJetAnalyzer
    *process.ak5PFJetAnalyzer
    *process.ak7PFJetAnalyzer
    )

# Schedule
process.schedule = cms.Schedule(#process.filter,
                                process.p
                                )


