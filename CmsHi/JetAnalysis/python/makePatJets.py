import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.files = 'rfio:/castor/cern.ch/user/y/yilmaz/pat/CMSSW_3_6_0/UnquenchedDijet80to120_runs1to500.root'
#ivars.files = 'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_7_0_pre2/RelValPyquen_DiJet_pt80to120_2760GeV/GEN-SIM-RECO/MC_37Y_V1-v1/0018/1420F252-8453-DF11-9903-002618943882.root'

ivars.output = 'jets_pat_80to120.root'
ivars.maxEvents = 3

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
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoHI.HiJetAlgos.HiGenJets_cff')
process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Summer09_7TeV_ReReco332")

process.ak5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak5CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK5Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_AK5Calo")))
process.ak5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.ak5match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak5CaloJets"),
                                                   matched = cms.InputTag("ak5clean"))
process.ak5parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak5CaloJets"))
process.ak5patJets = process.patJets.clone(jetSource  = cms.InputTag("ak5CaloJets"),
                                           genJetMatch = cms.InputTag("ak5match"),
                                           genPartonMatch= cms.InputTag("ak5parton"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5corr")))
process.ak5patSequence = cms.Sequence(process.ak5CaloJets*process.ak5corr*process.ak5clean*process.ak5match*process.ak5parton*process.ak5patJets)

process.ak7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak7CaloJets"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK7Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_AK7Calo")))
process.ak7clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak7HiGenJets'))
process.ak7match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak7CaloJets"),
                                                   matched = cms.InputTag("ak7clean"))
process.ak7parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak7CaloJets"))
process.ak7patJets = process.patJets.clone(jetSource  = cms.InputTag("ak7CaloJets"),
                                           genJetMatch = cms.InputTag("ak7match"),
                                           genPartonMatch= cms.InputTag("ak7parton"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak7corr")))
process.ak7patSequence = cms.Sequence(process.ak7CaloJets*process.ak7corr*process.ak7clean*process.ak7match*process.ak7parton*process.ak7patJets)

process.akPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK5Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_AK5Calo")))
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
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK7Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_AK7Calo")))
process.akPu7clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak7HiGenJets'))
process.akPu7match = process.patJetGenJetMatch.clone(src = cms.InputTag("akPu7CaloJets"),
                                                     matched = cms.InputTag("akPu7clean"))
process.akPu7parton = process.patJetPartonMatch.clone(src = cms.InputTag("akPu7CaloJets"))
process.akPu7patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu7CaloJets"),
                                             genJetMatch = cms.InputTag("akPu7match"),
                                             genPartonMatch= cms.InputTag("akPu7parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu7corr")))
process.akPu7patSequence = cms.Sequence(process.akPu7CaloJets*process.akPu7corr*process.akPu7clean*process.akPu7match*process.akPu7parton*process.akPu7patJets)

process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContent,
                                  fileName = cms.untracked.string(ivars.output)
                                  )

process.runAllJets = cms.Path(process.hiGenParticlesForJets *
                              process.hiRecoGenJets +
                              process.iterativeConePu5CaloJets +
                              process.makeHeavyIonJets +
                              process.akPu5patSequence +
                              process.akPu7patSequence +
                              process.ak5patSequence +
                              process.ak7patSequence
                              )

process.out_step = cms.EndPath(process.output)

