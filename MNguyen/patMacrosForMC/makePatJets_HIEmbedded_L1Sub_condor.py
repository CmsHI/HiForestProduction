import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('initialEvent',mult=ivars.multiplicity.singleton,info="for testing")

#ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/yetkin/test/mix/MinBias_DijetUnquenched80_d20101124/Mixed_DIGI_DEFF3756-E4F0-DF11-A905-003048F1BF66.root"
ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/yetkin/test/reco/MinBias0to20_DijetUnquenched50_d20101124/Mixed_Mixed_DIGI_FE622F66-EBF0-DF11-8447-003048F118D4.root"


#ivars.output = 'patJets_AMPT_Pyquen_DiJet_Pt80.root'
ivars.maxEvents = -1
ivars.initialEvent = 1

ivars.parseArguments()

import FWCore.ParameterSet.Config as cms


process = cms.Process('HIJETS')

#process.Timing = cms.Service("Timing")

#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#                                        ignoreTotal=cms.untracked.int32(0),
#                                        oncePerEventMode = cms.untracked.bool(False)
#                                        )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(ivars.files)
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(ivars.maxEvents)
                    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_39Y_V4::All'  #39X
process.GlobalTag.globaltag = 'START39_V6HI::All'  #39X

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

# This stuff no longer appears to be needed, but keep commented out just in case ...
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')


# load centrality

process.load("CommonTools.UtilAlgos.TFileService_cfi")

process.HeavyIonGlobalParameters = cms.PSet(
           centralityVariable = cms.string("HFhits"),
           nonDefaultGlauberModel = cms.string(""),
           centralitySrc = cms.InputTag("hiCentrality")
                     )

process.makeCentralityTableTFile = cms.EDAnalyzer('CentralityTableProducer',
                                                  isMC = cms.untracked.bool(False),
                                                  makeDBFromTFile = cms.untracked.bool(False),
                                                  makeTFileFromDB = cms.untracked.bool(True)
                                                  )

# for data!
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)



# --- Pat Stuff --

# is this just setting defaults that I don't use anyway?  What is DisableMonteCarloDeps doing?  
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Summer09_7TeV_ReReco332")



process.load("HeavyIonsAnalysis.Configuration.analysisProducers_cff")
process.hiExtra = cms.Sequence(
        process.allTracks +
            process.heavyIon
            )
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
    process.hiGenParticles +
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

process.ic5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic5sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC5Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_IC5Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ic5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone5HiGenJets'))
process.ic5match = process.patJetGenJetMatch.clone(src = cms.InputTag("ic5sub"),
                                                                                                      matched = cms.InputTag("ic5clean"))
process.ic5parton = process.patJetPartonMatch.clone(src = cms.InputTag("ic5sub"))
process.ic5patJets = process.patJets.clone(jetSource  = cms.InputTag("ic5sub"),
                                                                                      genJetMatch = cms.InputTag("ic5match"),
                                                                                      genPartonMatch= cms.InputTag("ic5parton"),
                                                                                      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic5corr")))
process.ic5patSequence = cms.Sequence(process.ic5CaloJets*process.ic5sub*process.ic5corr*process.ic5clean*process.ic5match*process.ic5parton*process.ic5patJets)

### AK5, Calo  ###

process.ak5sub = process.heavyIonL1SubtractedJets.clone()
process.ak5sub.src = 'ak5CaloJets'
process.ak5sub.rhoTag = 'kt4CaloJets'

process.ak5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak5sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK5Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_AK5Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ak5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.ak5match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak5sub"),
                                                   matched = cms.InputTag("ak5clean"))
process.ak5parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak5sub"))
process.ak5patJets = process.patJets.clone(jetSource  = cms.InputTag("ak5sub"),
                                           genJetMatch = cms.InputTag("ak5match"),
                                           genPartonMatch= cms.InputTag("ak5parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5corr")))
process.ak5patSequence = cms.Sequence(process.ak5CaloJets*process.ak5sub*process.ak5corr*process.ak5clean*process.ak5match*process.ak5parton*process.ak5patJets)

### AK7, Calo  ###

process.ak7sub = process.heavyIonL1SubtractedJets.clone()
process.ak7sub.src = 'ak7CaloJets'
process.ak7sub.rhoTag = 'kt4CaloJets'

process.ak7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak7sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK7Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_AK7Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ak7clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak7HiGenJets'))
process.ak7match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak7sub"),
                                                   matched = cms.InputTag("ak7clean"))
process.ak7parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak7sub"))
process.ak7patJets = process.patJets.clone(jetSource  = cms.InputTag("ak7sub"),
                                           genJetMatch = cms.InputTag("ak7match"),
                                           genPartonMatch= cms.InputTag("ak7parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak7corr")))
process.ak7patSequence = cms.Sequence(process.ak7CaloJets*process.ak7sub*process.ak7corr*process.ak7clean*process.ak7match*process.ak7parton*process.ak7patJets)

### KT4, Calo  ###

process.kt4sub = process.heavyIonL1SubtractedJets.clone()
process.kt4sub.src = 'kt4CaloJets'
process.kt4sub.rhoTag = 'kt4CaloJets'

process.kt4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt4sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT4Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_KT4Calo"),
                                                                        L5Flavor = cms.string("none")))
process.kt4clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt4HiGenJets'))
process.kt4match = process.patJetGenJetMatch.clone(src = cms.InputTag("kt4sub"),
                                                   matched = cms.InputTag("kt4clean"))
process.kt4parton = process.patJetPartonMatch.clone(src = cms.InputTag("kt4sub"))
process.kt4patJets = process.patJets.clone(jetSource  = cms.InputTag("kt4sub"),
                                           genJetMatch = cms.InputTag("kt4match"),
                                           genPartonMatch= cms.InputTag("kt4parton"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt4corr")))
process.kt4patSequence = cms.Sequence(process.kt4CaloJets*process.kt4sub*process.kt4corr*process.kt4clean*process.kt4match*process.kt4parton*process.kt4patJets)

### KT6, Calo  ###

process.kt6sub = process.heavyIonL1SubtractedJets.clone()
process.kt6sub.src = 'kt6CaloJets'
process.kt6sub.rhoTag = 'kt4CaloJets'

process.kt6corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt6sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT6Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_KT6Calo"),
                                                                        L5Flavor = cms.string("none")))
process.kt6clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt6HiGenJets'))
process.kt6match = process.patJetGenJetMatch.clone(src = cms.InputTag("kt6sub"),
                                                   matched = cms.InputTag("kt6clean"))
process.kt6parton = process.patJetPartonMatch.clone(src = cms.InputTag("kt6sub"))
process.kt6patJets = process.patJets.clone(jetSource  = cms.InputTag("kt6sub"),
                                           genJetMatch = cms.InputTag("kt6match"),
                                           genPartonMatch= cms.InputTag("kt6parton"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt6corr")))
process.kt6patSequence = cms.Sequence(process.kt6CaloJets*process.kt6sub*process.kt6corr*process.kt6clean*process.kt6match*process.kt6parton*process.kt6patJets)


process.akPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK5Calo"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_AK5Calo"),
                                                                          L5Flavor = cms.string("none")))
process.akPu5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.akPu5match = process.patJetGenJetMatch.clone(src = cms.InputTag("akPu5CaloJets"),
                                                     matched = cms.InputTag("akPu5clean"))
process.akPu5parton = process.patJetPartonMatch.clone(src = cms.InputTag("akPu5CaloJets"))
process.akPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu5CaloJets"),
                                             genJetMatch = cms.InputTag("akPu5match"),
                                             genPartonMatch= cms.InputTag("akPu5parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu5corr")))
process.akPu5patSequence = cms.Sequence(process.akPu5CaloJets*process.akPu5corr*process.akPu5clean*process.akPu5match*process.akPu5parton*process.akPu5patJets)



process.akPu7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu7CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK7Calo"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_AK7Calo"),
                                                                          L5Flavor = cms.string("none")))
process.akPu7clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak7HiGenJets'))
process.akPu7match = process.patJetGenJetMatch.clone(src = cms.InputTag("akPu7CaloJets"),
                                                     matched = cms.InputTag("akPu7clean"))
process.akPu7parton = process.patJetPartonMatch.clone(src = cms.InputTag("akPu7CaloJets"))
process.akPu7patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu7CaloJets"),
                                             genJetMatch = cms.InputTag("akPu7match"),
                                             genPartonMatch= cms.InputTag("akPu7parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu7corr")))
process.akPu7patSequence = cms.Sequence(process.akPu7CaloJets*process.akPu7corr*process.akPu7clean*process.akPu7match*process.akPu7parton*process.akPu7patJets)



process.ktPu4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ktPu4CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT4Calo"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_KT4Calo"),
                                                                          L5Flavor = cms.string("none")))
process.ktPu4clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt4HiGenJets'))
process.ktPu4match = process.patJetGenJetMatch.clone(src = cms.InputTag("ktPu4CaloJets"),
                                                     matched = cms.InputTag("ktPu4clean"))
process.ktPu4parton = process.patJetPartonMatch.clone(src = cms.InputTag("ktPu4CaloJets"))
process.ktPu4patJets = process.patJets.clone(jetSource  = cms.InputTag("ktPu4CaloJets"),
                                             genJetMatch = cms.InputTag("ktPu4match"),
                                             genPartonMatch= cms.InputTag("ktPu4parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ktPu4corr")))
process.ktPu4patSequence = cms.Sequence(process.ktPu4CaloJets*process.ktPu4corr*process.ktPu4clean*process.ktPu4match*process.ktPu4parton*process.ktPu4patJets)



process.ktPu6corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ktPu6CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT6Calo"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_KT6Calo"),
                                                                          L5Flavor = cms.string("none")))
process.ktPu6clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt6HiGenJets'))
process.ktPu6match = process.patJetGenJetMatch.clone(src = cms.InputTag("ktPu6CaloJets"),
                                                     matched = cms.InputTag("ktPu6clean"))
process.ktPu6parton = process.patJetPartonMatch.clone(src = cms.InputTag("ktPu6CaloJets"))
process.ktPu6patJets = process.patJets.clone(jetSource  = cms.InputTag("ktPu6CaloJets"),
                                             genJetMatch = cms.InputTag("ktPu6match"),
                                             genPartonMatch= cms.InputTag("ktPu6parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ktPu6corr")))
process.ktPu6patSequence = cms.Sequence(process.ktPu6CaloJets*process.ktPu6corr*process.ktPu6clean*process.ktPu6match*process.ktPu6parton*process.ktPu6patJets)

process.icPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("iterativeConePu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC5Calo"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_IC5Calo"),
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


#begin small cones

process.ak3CaloJets = process.ak5CaloJets.clone()
process.ak3CaloJets.rParam       = cms.double(0.3)
process.ak3CaloJets.radiusPU = 0.3

process.ak3sub = process.heavyIonL1SubtractedJets.clone()
process.ak3sub.src = 'ak3CaloJets'
process.ak3sub.rhoTag = 'kt4CaloJets'

process.ak3corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak3sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK3Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_AK3Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ak3clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiGenJets'))
process.ak3match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak3sub"),
                                                   matched = cms.InputTag("ak3clean"))
process.ak3parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak3sub"))
process.ak3patJets = process.patJets.clone(jetSource  = cms.InputTag("ak3sub"),
                                           genJetMatch = cms.InputTag("ak3match"),
                                           genPartonMatch= cms.InputTag("ak3parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3corr")))
process.ak3patSequence = cms.Sequence(process.ak3CaloJets*process.ak3sub*process.ak3corr*process.ak3clean*process.ak3match*process.ak3parton*process.ak3patJets)

process.ak4CaloJets = process.ak5CaloJets.clone()
process.ak4CaloJets.rParam       = cms.double(0.4)
process.ak4CaloJets.radiusPU = 0.4

process.ak4sub = process.heavyIonL1SubtractedJets.clone()
process.ak4sub.src = 'ak4CaloJets'
process.ak4sub.rhoTag = 'kt4CaloJets'

process.ak4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak4sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK4Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_AK4Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ak4clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiGenJets'))
process.ak4match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak4sub"),
                                                   matched = cms.InputTag("ak4clean"))
process.ak4parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak4sub"))
process.ak4patJets = process.patJets.clone(jetSource  = cms.InputTag("ak4sub"),
                                           genJetMatch = cms.InputTag("ak4match"),
                                           genPartonMatch= cms.InputTag("ak4parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4corr")))
process.ak4patSequence = cms.Sequence(process.ak4CaloJets*process.ak4sub*process.ak4corr*process.ak4clean*process.ak4match*process.ak4parton*process.ak4patJets)

process.kt3CaloJets = process.kt4CaloJets.clone()
process.kt3CaloJets.rParam       = cms.double(0.3)
process.kt3CaloJets.radiusPU = 0.3

process.kt3sub = process.heavyIonL1SubtractedJets.clone()
process.kt3sub.src = 'kt3CaloJets'
process.kt3sub.rhoTag = 'kt4CaloJets'

process.kt3corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt3sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT3Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_KT3Calo"),
                                                                        L5Flavor = cms.string("none")))
process.kt3clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt3HiGenJets'))
process.kt3match = process.patJetGenJetMatch.clone(src = cms.InputTag("kt3sub"),
                                                   matched = cms.InputTag("kt3clean"))
process.kt3parton = process.patJetPartonMatch.clone(src = cms.InputTag("kt3sub"))
process.kt3patJets = process.patJets.clone(jetSource  = cms.InputTag("kt3sub"),
                                           genJetMatch = cms.InputTag("kt3match"),
                                           genPartonMatch= cms.InputTag("kt3parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt3corr")))
process.kt3patSequence = cms.Sequence(process.kt3CaloJets*process.kt3sub*process.kt3corr*process.kt3clean*process.kt3match*process.kt3parton*process.kt3patJets)

process.ic3CaloJets = process.ic5CaloJets.clone()
process.ic3CaloJets.rParam       = cms.double(0.3)
process.ic3CaloJets.radiusPU = 0.3

process.ic3sub = process.heavyIonL1SubtractedJets.clone()
process.ic3sub.src = 'ic3CaloJets'
process.ic3sub.rhoTag = 'kt4CaloJets'

process.ic3corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic3sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC3Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_IC3Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ic3clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone3HiGenJets'))
process.ic3match = process.patJetGenJetMatch.clone(src = cms.InputTag("ic3sub"),
                                                   matched = cms.InputTag("ic3clean"))
process.ic3parton = process.patJetPartonMatch.clone(src = cms.InputTag("ic3sub"))
process.ic3patJets = process.patJets.clone(jetSource  = cms.InputTag("ic3sub"),
                                           genJetMatch = cms.InputTag("ic3match"),
                                           genPartonMatch= cms.InputTag("ic3parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic3corr")))
process.ic3patSequence = cms.Sequence(process.ic3CaloJets*process.ic3sub*process.ic3corr*process.ic3clean*process.ic3match*process.ic3parton*process.ic3patJets)


process.ic4sub = process.heavyIonL1SubtractedJets.clone()
process.ic4sub.src = 'ic4CaloJets'
process.ic4sub.rhoTag = 'ic4CaloJets'

process.ic4CaloJets = process.ic5CaloJets.clone()
process.ic4CaloJets.rParam       = cms.double(0.4)
process.ic4CaloJets.radiusPU = 0.4

process.ic4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic4sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC4Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_IC4Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ic4clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone4HiGenJets'))
process.ic4match = process.patJetGenJetMatch.clone(src = cms.InputTag("ic4sub"),
                                                   matched = cms.InputTag("ic4clean"))
process.ic4parton = process.patJetPartonMatch.clone(src = cms.InputTag("ic4sub"))
process.ic4patJets = process.patJets.clone(jetSource  = cms.InputTag("ic4sub"),
                                           genJetMatch = cms.InputTag("ic4match"),
                                           genPartonMatch= cms.InputTag("ic4parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic4corr")))
process.ic4patSequence = cms.Sequence(process.ic4CaloJets*process.ic4sub*process.ic4corr*process.ic4clean*process.ic4match*process.ic4parton*process.ic4patJets)





process.runAllJets = cms.Sequence(
    process.icPu5patSequence +
    process.kt4patSequence +
    process.kt6patSequence +
    process.ic5patSequence +
    process.ak5patSequence +
    process.ak7patSequence +
    process.ak3patSequence +
    process.ak4patSequence +
    process.kt3patSequence +
    process.ic3patSequence +
    process.ic4patSequence 
# Too slow!
#    process.akPu5patSequence +
#    process.akPu7patSequence +
#    process.ktPu4patSequence +
#    process.ktPu6patSequence +
)

#actually only works on raw b/c of it needs random numbers which aren't stored
#process.load("Configuration.StandardSequences.ValidationHeavyIons_cff")

# put it all together

process.path = cms.Path(
    process.makeCentralityTableTFile*
    process.hiExtra*
    process.hiGen*
    process.runAllJets
    #*process.hiTrackValidation
    )

process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContentMC,
                                  fileName = cms.untracked.string(ivars.output)
                                  )

#process.output.outputCommands.extend(["drop *_towerMaker_*_*"])
process.output.outputCommands.extend(["keep *_hiSelectedTracks_*_HIJETS"])
process.output.outputCommands.extend(["keep *_particleFlow_*_*"])
process.output.outputCommands.extend(["keep *_mergedtruth_*_*"])
process.output.outputCommands.extend(["keep double*_*PF*_*_*"])
process.output.outputCommands.extend(["keep *_heavyIon_*_*"])
process.output.outputCommands.extend(["keep *_generator_*_*"])
process.output.outputCommands.extend(["keep *_hiSignal_*_*"])
process.output.outputCommands.extend(["keep *_genParticles_*_*"])
process.output.outputCommands.extend(["keep *_hiGenParticles_*_*"])
process.output.outputCommands.extend(["keep *_TriggerResults_*_*"])
process.output.outputCommands.extend(["keep *_heavyIon_*_*"])
# reco jets
process.output.outputCommands.extend(["keep recoCaloJets_*_*_*"])
process.output.outputCommands.extend(["keep recoPFJets_*_*_*"])
#particle flow
process.output.outputCommands.extend(["keep recoPFClusters_*_*_*"])
process.output.outputCommands.extend(["keep recoPFRecHits_*_*_*"])
#fast jet pf stuff
process.output.outputCommands.extend(["keep doubles_*PF*_*_*"])
#calorimeter stuff
process.output.outputCommands.extend(["keep *_towerMaker_*_*"])
process.output.outputCommands.extend(["keep *_caloTowers_*_*"])
process.output.outputCommands.extend(["keep *_hcalnoise_*_*"])
process.output.outputCommands.extend(["keep *_hbhereco_*_*"])
process.output.outputCommands.extend(["keep *_horeco_*_*"])
process.output.outputCommands.extend(["keep *_hfreco_*_*"])
process.output.outputCommands.extend(["keep *_ecalRecHit_*_*"])




process.out_step = cms.EndPath(process.output)

