import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')


#process.load("HeavyIonsAnalysis.Configuration.Sources.Data2010_Run150476_Express_HLT_MinBias_cff")

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/cms/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/152/640/8604B150-C9FA-DF11-B748-0030487C912E.root")
                            #fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/e/edwenger/skims/mergeCleanDiJetall.root")
                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/cms/store/caf/user/edwenger/TestCentralSkim7/DiJet_100_1_Bce.root")

                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR10_P_V12::All'  #Data

process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


# --- Pat Stuff --
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

#  --- Track and muon reconstruction ---
# pixel triplet tracking (HI Tracking)
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")
process.load("RecoHI.HiTracking.LowPtTracking_PbPb_cff")

# Needed to produce "HcalSeverityLevelComputerRcd" used by CaloTowersCreator/towerMakerPF
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")

# Muon Reco
from RecoHI.HiMuonAlgos.HiRecoMuon_cff import * 
muons.JetExtractorPSet.JetCollectionLabel = cms.InputTag("iterativeConePu5CaloJets")

#Track Reco
process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking * muonRecoPbPb)

# for PF
process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")
#process.trackerDrivenElectronSeeds.UseQuality = cms.bool(False)
process.trackerDrivenElectronSeeds.TkColList = cms.VInputTag("hiGoodTracks")

# Load Heavy Ion "Good" Track Selection 
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")
process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')

# particle-flow stuff
process.HiParticleFlowRecoNoJets = cms.Sequence(
    process.particleFlowCluster
    * process.trackerDrivenElectronSeeds
    * process.particleFlowReco
    )

process.load("HeavyIonsAnalysis.Configuration.analysisProducers_cff")
process.hiExtra = cms.Sequence(
    process.allTracks 
    )

# Define Jet Algo parameters
process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')
process.load('RecoHI.HiJetAlgos.HiRecoPFJets_cff')
process.load("RecoJets.Configuration.RecoJPTJetsHIC_cff")

#Now the correct L1 
process.load("RecoHI.HiJetAlgos.HiL1Corrector_cff")


process.ic5sub = process.heavyIonL1SubtractedJets.clone()
process.ic5sub.src = 'ic5CaloJets'
process.ic5sub.rhoTag = 'kt4CaloJets'

process.patJets.embedCaloTowers = cms.bool(False)


# set up vertex matching
process.iterativeConePu5CaloJets.doPVCorrection = cms.bool(True)
process.iterativeConePu5CaloJets.srcPVs = 'hiSelectedVertex'

process.ic5CaloJets.doPVCorrection = cms.bool(True)
process.ic5CaloJets.srcPVs = 'hiSelectedVertex'

process.ak5CaloJets.doPVCorrection = cms.bool(True)
process.ak5CaloJets.srcPVs = 'hiSelectedVertex'

process.ak7CaloJets.doPVCorrection = cms.bool(True)
process.ak7CaloJets.srcPVs = 'hiSelectedVertex'

process.kt4CaloJets.doPVCorrection = cms.bool(True)
process.kt4CaloJets.srcPVs = 'hiSelectedVertex'

process.kt6CaloJets.doPVCorrection = cms.bool(True)
process.kt6CaloJets.srcPVs = 'hiSelectedVertex'


process.ic5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic5sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_IC5Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_IC5Calo"),
                                                                        L5Flavor = cms.string("none")))

process.ic5patJets = process.patJets.clone(jetSource  = cms.InputTag("ic5sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic5corr")))

process.ic5patSequence = cms.Sequence(process.ic5CaloJets*process.ic5sub*process.ic5corr*process.ic5patJets)


process.ak5sub = process.heavyIonL1SubtractedJets.clone()
process.ak5sub.src = 'ak5CaloJets'
process.ak5sub.rhoTag = 'kt4CaloJets'


process.ak5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak5sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK5Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_AK5Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ak5patJets = process.patJets.clone(jetSource  = cms.InputTag("ak5sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5corr")))
process.ak5patSequence = cms.Sequence(process.ak5CaloJets*process.ak5sub*process.ak5corr*process.ak5patJets)



process.ak7sub = process.heavyIonL1SubtractedJets.clone()
process.ak7sub.src = 'ak7CaloJets'
process.ak7sub.rhoTag = 'kt4CaloJets'


process.ak7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak7sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK7Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_AK7Calo"),
                                                                        L5Flavor = cms.string("none")))
process.ak7patJets = process.patJets.clone(jetSource  = cms.InputTag("ak7sub"),

                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak7corr")))
process.ak7patSequence = cms.Sequence(process.ak7CaloJets*process.ak7sub*process.ak7corr*process.ak7patJets)



process.kt4sub = process.heavyIonL1SubtractedJets.clone()
process.kt4sub.src = 'kt4CaloJets'
process.kt4sub.rhoTag = 'kt4CaloJets'


process.kt4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt4sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT4Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_KT4Calo"),
                                                                        L5Flavor = cms.string("none")))
process.kt4patJets = process.patJets.clone(jetSource  = cms.InputTag("kt4sub"),

                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt4corr")))
process.kt4patSequence = cms.Sequence(process.kt4CaloJets*process.kt4sub*process.kt4corr*process.kt4patJets)


process.kt6sub = process.heavyIonL1SubtractedJets.clone()
process.kt6sub.src = 'kt6CaloJets'
process.kt6sub.rhoTag = 'kt4CaloJets'


process.kt6corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt6sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT6Calo"),
                                                                        L3Absolute = cms.string("L3Absolute_KT6Calo"),
                                                                        L5Flavor = cms.string("none")))

process.kt6patJets = process.patJets.clone(jetSource  = cms.InputTag("kt6sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt6corr")))
process.kt6patSequence = cms.Sequence(process.kt6CaloJets*process.kt6sub*process.kt6corr*process.kt6patJets)



process.ak5PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak5PFsub.src = 'ak5PFJets'
process.ak5PFsub.rhoTag = 'kt4PFJets'
process.ak5PFsub.jetType    = 'PFJet'



process.ak5PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak5PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK5PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_AK5PF"),
                                                                          L5Flavor = cms.string("none")))

process.ak5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak5PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFcorr")))
process.ak5PFpatSequence = cms.Sequence(process.ak5PFJets*process.ak5PFsub*process.ak5PFcorr*process.ak5PFpatJets)



process.ak7PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak7PFsub.src = 'ak7PFJets'
process.ak7PFsub.rhoTag = 'kt4PFJets'
process.ak7PFsub.jetType    = 'PFJet'


process.ak7PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak7PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK7PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_AK7PF"),
                                                                          L5Flavor = cms.string("none")))
process.ak7PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak7PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak7PFcorr")))
process.ak7PFpatSequence = cms.Sequence(process.ak7PFJets*process.ak7PFsub*process.ak7PFcorr*process.ak7PFpatJets)


process.kt4PFsub = process.heavyIonL1SubtractedJets.clone()
process.kt4PFsub.src = 'kt4PFJets'
process.kt4PFsub.rhoTag = 'kt4PFJets'
process.kt4PFsub.jetType    = 'PFJet'


process.kt4PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt4PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT4PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_KT4PF"),
                                                                          L5Flavor = cms.string("none")))
process.kt4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("kt4PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt4PFcorr")))
process.kt4PFpatSequence = cms.Sequence(process.kt4PFJets*process.kt4PFsub*process.kt4PFcorr*process.kt4PFpatJets)



process.kt6PFsub = process.heavyIonL1SubtractedJets.clone()
process.kt6PFsub.src = 'kt6PFJets'
process.kt6PFsub.rhoTag = 'kt4PFJets'
process.kt6PFsub.jetType    = 'PFJet'

process.kt6PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt6PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT6PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_KT6PF"),
                                                                          L5Flavor = cms.string("none")))
process.kt6PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("kt6PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt6PFcorr")))
process.kt6PFpatSequence = cms.Sequence(process.kt6PFJets*process.kt6PFsub*process.kt6PFcorr*process.kt6PFpatJets)



process.ic5PFsub = process.heavyIonL1SubtractedJets.clone()
process.ic5PFsub.src = 'ic5PFJets'
process.ic5PFsub.rhoTag = 'kt4PFJets'
process.ic5PFsub.jetType    = 'PFJet'

process.ic5PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic5PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC5PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_IC5PF"),
                                                                          L5Flavor = cms.string("none")))
process.ic5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ic5PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic5PFcorr")))
process.ic5PFpatSequence = cms.Sequence(process.ic5PFJets*process.ic5PFsub*process.ic5PFcorr*process.ic5PFpatJets)




process.akPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK5Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_AK5Calo"),
                                                                          L5Flavor = cms.string("none")))
process.akPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu5CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu5corr")))
process.akPu5patSequence = cms.Sequence(process.akPu5CaloJets*process.akPu5corr*process.akPu5patJets)


process.akPu7corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("akPu7CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_AK7Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_AK7Calo"),
                                                                          L5Flavor = cms.string("none")))
process.akPu7patJets = process.patJets.clone(jetSource  = cms.InputTag("akPu7CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu7corr")))
process.akPu7patSequence = cms.Sequence(process.akPu7CaloJets*process.akPu7corr*process.akPu7patJets)



process.ktPu4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ktPu4CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT4Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_KT4Calo"),
                                                                          L5Flavor = cms.string("none")))
process.ktPu4patJets = process.patJets.clone(jetSource  = cms.InputTag("ktPu4CaloJets"),
                                          
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ktPu4corr")))
process.ktPu4patSequence = cms.Sequence(process.ktPu4CaloJets*process.ktPu4corr*process.ktPu4patJets)


process.ktPu6corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ktPu6CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_KT6Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_KT6Calo"),
                                                                          L5Flavor = cms.string("none")))
process.ktPu6patJets = process.patJets.clone(jetSource  = cms.InputTag("ktPu6CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ktPu6corr")))
process.ktPu6patSequence = cms.Sequence(process.ktPu6CaloJets*process.ktPu6corr*process.ktPu6patJets)

process.icPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("iterativeConePu5CaloJets"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("L2Relative_IC5Calo"),
                                                                          L3Absolute = cms.string("L3Absolute_IC5Calo"),
                                                                          L5Flavor = cms.string("none")))
process.icPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("iterativeConePu5CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5corr")))
process.icPu5patSequence = cms.Sequence(process.iterativeConePu5CaloJets*process.icPu5corr*process.icPu5patJets)


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
process.ak3patJets = process.patJets.clone(jetSource  = cms.InputTag("ak3sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3corr")))
process.ak3patSequence = cms.Sequence(process.ak3CaloJets*process.ak3sub*process.ak3corr*process.ak3patJets)

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
process.ak4patJets = process.patJets.clone(jetSource  = cms.InputTag("ak4sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4corr")))
process.ak4patSequence = cms.Sequence(process.ak4CaloJets*process.ak4sub*process.ak4corr*process.ak4patJets)

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
process.kt3patJets = process.patJets.clone(jetSource  = cms.InputTag("kt3sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt3corr")))
process.kt3patSequence = cms.Sequence(process.kt3CaloJets*process.kt3sub*process.kt3corr*process.kt3patJets)

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
process.ic3patJets = process.patJets.clone(jetSource  = cms.InputTag("ic3sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic3corr")))
process.ic3patSequence = cms.Sequence(process.ic3CaloJets*process.ic3sub*process.ic3corr*process.ic3patJets)


process.ic4CaloJets = process.ic5CaloJets.clone()
process.ic4CaloJets.rParam       = cms.double(0.4)
process.ic4CaloJets.radiusPU = 0.4


process.ic4sub = process.heavyIonL1SubtractedJets.clone()
process.ic4sub.src = 'ic4CaloJets'
process.ic4sub.rhoTag = 'kt4CaloJets'


process.ic4corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic4sub"),
                                                  corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC4Calo"),
                                                                        L3Absolute = cms.string("HI_L3Absolute_IC4Calo"),
                                                                        L5Flavor = cms.string("none")))

process.ic4patJets = process.patJets.clone(jetSource  = cms.InputTag("ic4sub"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic4corr")))
process.ic4patSequence = cms.Sequence(process.ic4CaloJets*process.ic4sub*process.ic4corr*process.ic4patJets)


process.ak3PFJets = process.ak5PFJets.clone()
process.ak3PFJets.rParam       = cms.double(0.3)
process.ak3PFJets.radiusPU = 0.3


process.ak3PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak3PFsub.src = 'ak3PFJets'
process.ak3PFsub.rhoTag = 'kt4PFJets'
process.ak3PFsub.jetType    = 'PFJet'

process.ak3PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak3PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK3PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_AK3PF"),
                                                                          L5Flavor = cms.string("none")))

process.ak3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak3PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3PFcorr")))
process.ak3PFpatSequence = cms.Sequence(process.ak3PFJets*process.ak3PFsub*process.ak3PFcorr*process.ak3PFpatJets)

process.ak4PFJets = process.ak5PFJets.clone()
process.ak4PFJets.rParam       = cms.double(0.4)
process.ak4PFJets.radiusPU = 0.4


process.ak4PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak4PFsub.src = 'ak4PFJets'
process.ak4PFsub.rhoTag = 'kt4PFJets'
process.ak4PFsub.jetType    = 'PFJet'

process.ak4PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ak4PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_AK4PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_AK4PF"),
                                                                          L5Flavor = cms.string("none")))

process.ak4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak4PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4PFcorr")))
process.ak4PFpatSequence = cms.Sequence(process.ak4PFJets*process.ak4PFsub*process.ak4PFcorr*process.ak4PFpatJets)


process.ic3PFJets = process.ic5PFJets.clone()
process.ic3PFJets.rParam       = cms.double(0.3)
process.ic3PFJets.radiusPU = 0.3


process.ic3PFsub = process.heavyIonL1SubtractedJets.clone()
process.ic3PFsub.src = 'ic3PFJets'
process.ic3PFsub.rhoTag = 'kt4PFJets'
process.ic3PFsub.jetType    = 'PFJet'


process.ic3PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic3PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC3PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_IC3PF"),
                                                                          L5Flavor = cms.string("none")))
process.ic3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ic3PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic3PFcorr")))
process.ic3PFpatSequence = cms.Sequence(process.ic3PFJets*process.ic3PFsub*process.ic3PFcorr*process.ic3PFpatJets)

process.ic4PFJets = process.ic5PFJets.clone()
process.ic4PFJets.rParam       = cms.double(0.4)
process.ic4PFJets.radiusPU = 0.4


process.ic4PFsub = process.heavyIonL1SubtractedJets.clone()
process.ic4PFsub.src = 'ic4PFJets'
process.ic4PFsub.rhoTag = 'kt4PFJets'
process.ic4PFsub.jetType    = 'PFJet'

process.ic4PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("ic4PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_IC4PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_IC4PF"),
                                                                          L5Flavor = cms.string("none")))
process.ic4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ic4PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic4PFcorr")))
process.ic4PFpatSequence = cms.Sequence(process.ic4PFJets*process.ic4PFsub*process.ic4PFcorr*process.ic4PFpatJets)


process.kt3PFJets = process.kt4PFJets.clone()
process.kt3PFJets.rParam       = cms.double(0.3)
process.kt3PFJets.radiusPU = 0.3



process.kt3PFsub = process.heavyIonL1SubtractedJets.clone()
process.kt3PFsub.src = 'kt3PFJets'
process.kt3PFsub.rhoTag = 'kt4PFJets'
process.kt3PFsub.jetType    = 'PFJet'


process.kt3PFcorr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("kt3PFsub"),
                                                    corrLevels = cms.PSet(L2Relative = cms.string("HI_L2Relative_KT3PF"),
                                                                          L3Absolute = cms.string("HI_L3Absolute_KT3PF"),
                                                                          L5Flavor = cms.string("none")))
process.kt3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("kt3PFsub"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt3PFcorr")))
process.kt3PFpatSequence = cms.Sequence(process.kt3PFJets*process.kt3PFsub*process.kt3PFcorr*process.kt3PFpatJets)


# JPT

process.jpticPu5corr = process.patJetCorrFactors.clone(jetSource = cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                       corrLevels =  cms.PSet(L2Relative = cms.string("L2Relative_AK5JPT"),
                                                                              L3Absolute = cms.string("L3Absolute_AK5JPT"),
                                                                              L5Flavor = cms.string("none")))
process.jpticPu5patJets = process.patJets.clone(jetSource  =cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jpticPu5corr")))
process.icPu5JPTpatSequence = cms.Sequence(process.recoJPTJetsHIC*process.jpticPu5corr*process.jpticPu5patJets)

# set JPT to look at the right track collection:
process.JPTiterativeConePu5JetTracksAssociatorAtVertex.tracks = 'hiGoodMergedTracks'
process.JPTiterativeConePu5JetTracksAssociatorAtCaloFace.tracks = 'hiGoodMergedTracks'
process.JetPlusTrackZSPCorJetIconePu5.tracks = 'hiGoodMergedTracks'



# Run kT4 first if you want to use the rho for other algorithms
process.runAllJets = cms.Sequence(    
    process.kt4patSequence +
    process.ic5patSequence +
    #process.kt6patSequence +
    process.ak5patSequence +
    #process.ak7patSequence +
    process.kt4PFpatSequence +
    process.ak5PFpatSequence +
    #process.ak7PFpatSequence +
    #process.kt6PFpatSequence +
    process.ic5PFpatSequence +
    process.icPu5patSequence 
    #process.ak3patSequence +
    #process.ak4patSequence +
    #process.kt3patSequence +
    #process.ic3patSequence +
    #process.ic4patSequence +
    #process.ak3PFpatSequence +
    #process.ak4PFpatSequence +
    #process.kt3PFpatSequence +
    #process.ic3PFpatSequence +
    #process.ic4PFpatSequence + 
# Too slow!
#    process.akPu5patSequence +
#    process.akPu7patSequence +
#    process.ktPu4patSequence +
#    process.ktPu6patSequence +
#process.icPu5JPTpatSequence
)


# jet analysis trees
# trees need centrality
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


process.load("MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_cff")
process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_cff")
process.PFJetAnalyzer.trackTag  = cms.InputTag("hiGoodTracks")




# track efficiency anlayzer
#process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")
#for tree output
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("JetAnalysisTTrees_hiGoodTracks_TrackCaloMatching.root"))



process.load("MNguyen.Configuration.HI_JetSkim_cff")
process.hltJetHI.HLTPaths = ["HLT_HIJet50U"]


process.jetSkimPath = cms.Path(
    process.jetSkimSequence*
    process.hiTrackReco*
    process.hiGoodTracksSelection*
    process.HiParticleFlowRecoNoJets*
    process.runAllJets*
    process.PFJetAnalyzerSequence
    )




# Schedule definition
process.schedule = cms.Schedule(process.jetSkimPath)


# And the logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring('Unknown', 
        'ProductNotFound', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ModuleFailure', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError', 
        'FatalRootError', 
        'NotFound')
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100









