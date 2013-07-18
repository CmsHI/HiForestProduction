import FWCore.ParameterSet.Config as cms


process = cms.Process('HIJETS')

#process.Timing = cms.Service("Timing")

#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#                                        ignoreTotal=cms.untracked.int32(0),
#                                        oncePerEventMode = cms.untracked.bool(False)
#                                        )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FED91C44-9AFA-DF11-B174-0002C90B3976.root")
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
                    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_39Y_V4::All'  #39X
process.GlobalTag.globaltag = 'START39_V6HI::All'  #39X

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')


# load centrality

process.load("CommonTools.UtilAlgos.TFileService_cfi")

process.HeavyIonGlobalParameters = cms.PSet(
           centralityVariable = cms.string("HFhits"),
           #nonDefaultGlauberModel = cms.string("AMPT_2760GeV"),
           nonDefaultGlauberModel = cms.string("Hydjet_2760GeV"),
           centralitySrc = cms.InputTag("hiCentrality")
                     )

process.makeCentralityTableTFile = cms.EDAnalyzer('CentralityTableProducer',
                                                  isMC = cms.untracked.bool(True),
                                                  makeDBFromTFile = cms.untracked.bool(False),
                                                  makeTFileFromDB = cms.untracked.bool(True)
                                                  )


# --- Pat Stuff --

# is this just setting defaults that I don't use anyway?  What is DisableMonteCarloDeps doing?  
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Summer09_7TeV_ReReco332")



# --- Gen stuff ---
# Note that we need to re-run gen jets b/c until AK-IC bug fix propagates to production
process.load('PhysicsTools.HepMCCandAlgos.HiGenParticles_cfi')

process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoHI.HiJetAlgos.HiGenJets_cff')


#This runs GenHIEventProducer
process.heavyIon

#process.hiGenParticles.srcVector = cms.vstring('hiSignal')

process.hiGen = cms.Sequence(
#Careful when using embedded samples
    #process.hiGenParticles +
    process.hiGenParticlesForJets *
    process.hiRecoGenJets +
    process.genPartons *
    process.hiPartons
    )



# Define Jet Algo parameters

process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')

process.patJets.addGenPartonMatch   = True
process.patJets.addGenJetMatch      = True

process.load("RecoHI.HiJetAlgos.HiL1Corrector_cff")

process.patJets.embedCaloTowers = cms.bool(False)

### IC5, Calo  ###

process.ic5sub = process.heavyIonL1SubtractedJets.clone()
process.ic5sub.src = 'ic5CaloJets'
process.ic5sub.rhoTag = 'kt4CaloJets'

process.icPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("iterativeConePu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_IC5Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_IC5Calo"),
                                                                          L5Flavor = cms.string("none")))
process.icPu5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone5HiGenJets'))
process.icPu5match = process.patJetGenJetMatch.clone(src = cms.InputTag("iterativeConePu5CaloJets"),
                                                     matched = cms.InputTag("icPu5clean"))
process.icPu5parton = process.patJetPartonMatch.clone(src = cms.InputTag("iterativeConePu5CaloJets"))
process.icPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("iterativeConePu5CaloJets"),
                                             genJetMatch = cms.InputTag("icPu5match"),
                                             genPartonMatch= cms.InputTag("icPu5parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5corr")))
process.icPu5patSequence = cms.Sequence(process.iterativeConePu5CaloJets*process.icPu5corr*process.icPu5clean*process.icPu5match*process.icPu5parton*process.icPu5patJets)




process.runAllJets = cms.Sequence(
    process.icPu5patSequence 
)


process.path = cms.Path(
    process.hiGen*
    process.runAllJets
    )

process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContentMC,
                                  fileName = cms.untracked.string("PAT.root")
                                  )

process.output.outputCommands.extend(["keep *_genParticles_*_*"])
process.output.outputCommands.extend(["keep *_hiGenParticles_*_*"])

process.output.outputCommands.extend(["keep recoCaloJets_*_*_*"])



process.out_step = cms.EndPath(process.output)

