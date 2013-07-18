import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('initialEvent',mult=ivars.multiplicity.singleton,info="for testing")

ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/HICorePhysics/Jet50-PromptReco-OfficialJSON_hiGoodMergedTracks_Runs_152652_to_152957_RECOPAT-v1/760bdb5085e074829f50ca45fe83f6ca/RECOPAT_255_1_HCm.root"

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


# run one jet algo, but others are now added in the source code
# Don't change the jet collections, it's not configurable for the moment

process.PFJetAnalyzer = cms.EDAnalyzer("PFJetAnalyzer",
                                       jetTag = cms.InputTag("icPu5patJets"),
                                       jetTag2 = cms.InputTag("ic5PFpatJets"),
                                       jetTag3 = cms.InputTag("ak5PFpatJets"),
                                       jetTag4 = cms.InputTag("ic5patJets"),
                                       recoJetTag = cms.InputTag("iterativeConePu5CaloJets"),
                                       recoJetTag2 = cms.InputTag("ic5PFJets"),
                                       recoJetTag3 = cms.InputTag("ak5PFJets"),
                                       recoJetTag4 = cms.InputTag("ic5Jets"),
                                       pfCandidatesTag  = cms.InputTag("particleFlow",""),
                                       trackTag  = cms.InputTag("hiGoodMergedTracks"),
                                       isMC = cms.untracked.bool(False), 
                                       useCentrality = cms.untracked.bool(True)
                                       )



# ntuple output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string(ivars.output))
# Paths

process.jetSkimPath = cms.Path(process.jetSkimSequence)

process.p = cms.Path(    process.jetSkimSequence
                         *process.PFJetAnalyzer    
                         )

## Schedule
process.schedule = cms.Schedule(
    process.p
    )


