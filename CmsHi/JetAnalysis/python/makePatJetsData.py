import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('initialEvent',mult=ivars.multiplicity.singleton,info="for testing")
#ivars.files = 'dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/nart/MinimumBiasHI/HighPtMinBiasAnaSkimHardEnriched_01/6edc9684fb662a4038079f7b1d3fbd82/minbiashardenriched_2_1_wUe.root'
#ivars.files = 'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_7_0_pre2/RelValPyquen_DiJet_pt80to120_2760GeV/GEN-SIM-RECO/MC_37Y_V1-v1/0018/1420F252-8453-DF11-9903-002618943882.root'

ivars.files = "dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/nart/MinimumBiasHI/HighPtHardEnrichedAnaSkimNZS_03/b651c142d334b6448e9d198ad759ddac/allenrichedNZS_4_1_OCO.root"

ivars.output = 'jets_pat_80to120.root'
ivars.maxEvents = -1
ivars.initialEvent = 1

ivars.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(ivars.files)
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(ivars.maxEvents)
                    )

process.load('Configuration.StandardSequences.Geometry_cff')
process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Summer09_7TeV_ReReco332")

process.patJets.addGenPartonMatch   = False
process.patJets.addGenJetMatch      = False

process.icPu5CaloJets = process.iterativeConePu5CaloJets.clone()
process.icPu5corr = process.patJetCorrFactors.clone()
process.icPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("icPu5CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_IC5Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_IC5Calo")))
process.icPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("icPu5CaloJets"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5corr")))

process.icPu5patSequence = cms.Sequence(process.icPu5CaloJets*process.icPu5corr*process.icPu5patJets)

process.icPu7CaloJets = process.iterativeConePu7CaloJets.clone()
process.icPu7corr = process.patJetCorrFactors.clone()
process.icPu7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("icPu7CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_IC5Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_IC5Calo")))
process.icPu7patJets = process.patJets.clone(jetSource  = cms.InputTag("icPu7CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu7corr")))

process.icPu7patSequence = cms.Sequence(process.icPu7CaloJets*process.icPu7corr*process.icPu7patJets)

process.ak5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak5CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK5Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_AK5Calo")))
process.ak5patJets = process.patJets.clone(jetSource  = cms.InputTag("ak5CaloJets"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5corr")))
process.ak5patSequence = cms.Sequence(process.ak5CaloJets*process.ak5corr*process.ak5patJets)


process.ak7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak7CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK7Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_AK7Calo")))
process.ak7patJets = process.patJets.clone(jetSource  = cms.InputTag("ak7CaloJets"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak7corr")))
process.ak7patSequence = cms.Sequence(process.ak7CaloJets*process.ak7corr*process.ak7patJets)

process.akPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK5Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_AK5Calo")))
process.akPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu5CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu5corr")))
process.akPu5patSequence = cms.Sequence(process.akPu5CaloJets*process.akPu5corr*process.akPu5patJets)


process.akPu7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu7CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK7Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_AK7Calo")))
process.akPu7patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu7CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu7corr")))
process.akPu7patSequence = cms.Sequence(process.akPu7CaloJets*process.akPu7corr*process.akPu7patJets)

process.kt4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt4CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT4Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_KT4Calo")))
process.kt4patJets = process.patJets.clone(jetSource  = cms.InputTag("kt4CaloJets"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt4corr")))
process.kt4patSequence = cms.Sequence(process.kt4CaloJets*process.kt4corr*process.kt4patJets)

process.kt6corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt6CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT6Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_KT6Calo")))
process.kt6patJets = process.patJets.clone(jetSource  = cms.InputTag("kt6CaloJets"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt6corr")))

process.kt6patSequence = cms.Sequence(process.kt6CaloJets*process.kt6corr*process.kt6patJets)

process.ktPu4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ktPu4CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT4Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_KT4Calo")))
process.ktPu4patJets = process.patJets.clone(jetSource  = cms.InputTag("ktPu4CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ktPu4corr")))
process.ktPu4patSequence = cms.Sequence(process.ktPu4CaloJets*process.ktPu4corr*process.ktPu4patJets)

process.ktPu6corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ktPu6CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT6Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_KT6Calo")))
process.ktPu6patJets = process.patJets.clone(jetSource  = cms.InputTag("ktPu6CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ktPu6corr")))

process.ktPu6patSequence = cms.Sequence(process.ktPu6CaloJets*process.ktPu6corr*process.ktPu6patJets)

process.load("HeavyIonsAnalysis.Configuration.analysisProducers_cff")
process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContent,
                                  fileName = cms.untracked.string(ivars.output)
                                  )

process.output.outputCommands.extend(["drop *_towerMaker_*_*"])

process.runAllJets = cms.Path(
    process.allTracks +
    process.icPu5patSequence +
    process.icPu7patSequence +
#    process.ktPu4patSequence +
#    process.ktPu6patSequence +
    process.kt4patSequence +
    process.kt6patSequence +
#    process.akPu5patSequence +
#    process.akPu7patSequence +
    process.ak5patSequence +
    process.ak7patSequence 
    )

process.out_step = cms.EndPath(process.output)

