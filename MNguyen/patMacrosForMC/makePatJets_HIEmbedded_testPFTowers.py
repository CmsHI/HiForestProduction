import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')


# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_3_10_0/RelValHydjetQ_B0_2760GeV/GEN-SIM-RECO/START310_V1HI-v1/0050/14EAEC8E-680E-E011-AF49-00261894395A.root'
                                                              #'/store/relval/CMSSW_3_10_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/START310_V3-v1/0050/54ED7C66-680E-E011-BD9F-001A92810ADE.root'
                                                              )
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START310_V1HI::All' 

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

# load centrality

process.load("CommonTools.UtilAlgos.TFileService_cfi")

process.HeavyIonGlobalParameters = cms.PSet(
           centralityVariable = cms.string("HFhits"),
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
#configureHeavyIons(process)

process.patJets.jetSource  = cms.InputTag("iterativeConePu5CaloJets")


process.patJets.addBTagInfo         = False
process.patJets.addTagInfos         = False
process.patJets.addDiscriminators   = False
process.patJets.addAssociatedTracks = False
process.patJets.addJetCharge        = False
process.patJets.addJetID            = False
process.patJets.getJetMCFlavour     = False
process.patJets.addGenPartonMatch   = True
process.patJets.addGenJetMatch      = True
process.patJets.embedGenJetMatch    = True
process.patJets.embedGenPartonMatch   = True

process.patJetPartonMatch.matched = cms.InputTag("hiGenParticles")


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
#muons.JetExtractorPSet.JetCollectionLabel = cms.InputTag("iterativeCone5CaloJets")

#Track Reco
process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking * muonRecoPbPb)

# good track selection
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")
# merge with pixel tracks
process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')

# for PF
process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")


# particle-flow stuff

process.HiParticleFlowRecoNoJets = cms.Sequence(
    process.particleFlowCluster
    * process.trackerDrivenElectronSeeds
    * process.particleFlowReco
    )
process.trackerDrivenElectronSeeds.TkColList = cms.VInputTag("hiGoodMergedTracks")

process.load("HeavyIonsAnalysis.Configuration.analysisProducers_cff")
process.hiExtra = cms.Sequence(
            process.heavyIon
            )
#This runs GenHIEventProducer
process.heavyIon


# --- Gen stuff ---
# Note that we need to re-run gen jets b/c until AK-IC bug fix propagates to production
process.load('PhysicsTools.HepMCCandAlgos.HiGenParticles_cfi')

process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoHI.HiJetAlgos.HiGenJets_cff')



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
process.load('RecoHI.HiJetAlgos.HiRecoPFJets_cff')
process.load("RecoJets.Configuration.RecoJPTJetsHIC_cff")

process.patJets.addGenPartonMatch   = True
process.patJets.addGenJetMatch      = True

process.load("RecoHI.HiJetAlgos.HiL1Corrector_cff")

process.patJets.embedCaloTowers = cms.bool(False)

# set up vertex matching, does it work for embedded events?
process.iterativeConePu5CaloJets.doPVCorrection = cms.bool(True)
process.iterativeConePu5CaloJets.srcPVs = 'hiSelectedVertex'

process.ic5CaloJets.doPVCorrection = cms.bool(True)
process.ic5CaloJets.srcPVs = 'hiSelectedVertex'

process.ak5CaloJets.doPVCorrection = cms.bool(True)
process.ak5CaloJets.srcPVs = 'hiSelectedVertex'

process.kt4CaloJets.doPVCorrection = cms.bool(True)
process.kt4CaloJets.srcPVs = 'hiSelectedVertex'

#process.iterativeConePu5CaloJets.inputEtMin = 2.0
#process.ic5CaloJets.inputEtMin = 2.0
#process.ak5CaloJets.inputEtMin = 2.0
#process.ak7CaloJets.inputEtMin = 2.0
#process.kt4CaloJets.inputEtMin = 2.0
#process.kt6CaloJets.inputEtMin = 2.0
#process.ic5PFJets.inputEtMin = 2.0
#process.ak5PFJets.inputEtMin = 2.0
#process.ak7PFJets.inputEtMin = 2.0
#process.kt4PFJets.inputEtMin = 2.0
#process.kt6PFJets.inputEtMin = 2.0

process.iterativeConePu5CaloJets.puWidth = 0.5
process.ic5CaloJets.puWidth = 0.5
process.ak5CaloJets.puWidth = 0.5
process.ak7CaloJets.puWidth = 0.5
process.kt4CaloJets.puWidth = 0.5
process.kt6CaloJets.puWidth = 0.5
process.ic5PFJets.puWidth = 0.5
process.ak5PFJets.puWidth = 0.5
process.ak7PFJets.puWidth = 0.5
process.kt4PFJets.puWidth = 0.5
process.kt6PFJets.puWidth = 0.5


# Fist are Algorithms using FastJet-style PU subtraction with KT4 used to estimate rho

### IC5, Calo  ###

process.ic5sub = process.heavyIonL1SubtractedJets.clone()
process.ic5sub.src = 'ic5CaloJets'
process.ic5sub.rhoTag = 'kt4CaloJets'

process.ic5corr = process.patJetCorrFactors.clone(src = cms.InputTag("ic5sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
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

process.ak5corr = process.patJetCorrFactors.clone(src = cms.InputTag("ak5sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5Calo')
                                                  )
process.ak5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.ak5match = process.patJetGenJetMatch.clone(src = cms.InputTag("ak5sub"),
                                                   matched = cms.InputTag("ak5clean"))
process.ak5parton = process.patJetPartonMatch.clone(src = cms.InputTag("ak5sub"))
process.ak5patJets = process.patJets.clone(jetSource  = cms.InputTag("ak5sub"),
                                           genJetMatch = cms.InputTag("ak5match"),
                                           genPartonMatch= cms.InputTag("ak5parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5corr")))
process.ak5patSequence = cms.Sequence(process.ak5CaloJets*process.ak5sub*process.ak5corr*process.ak5clean*process.ak5match*process.ak5parton*process.ak5patJets)


### KT4, Calo  ###

process.kt4sub = process.heavyIonL1SubtractedJets.clone()
process.kt4sub.src = 'kt4CaloJets'
process.kt4sub.rhoTag = 'kt4CaloJets'

process.kt4corr = process.patJetCorrFactors.clone(src = cms.InputTag("kt4sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('KT4Calo')
                                                  )
process.kt4clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt4HiGenJets'))
process.kt4match = process.patJetGenJetMatch.clone(src = cms.InputTag("kt4sub"),
                                                   matched = cms.InputTag("kt4clean"))
process.kt4parton = process.patJetPartonMatch.clone(src = cms.InputTag("kt4sub"))
process.kt4patJets = process.patJets.clone(jetSource  = cms.InputTag("kt4sub"),
                                           genJetMatch = cms.InputTag("kt4match"),
                                           genPartonMatch= cms.InputTag("kt4parton"),
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt4corr")))
process.kt4patSequence = cms.Sequence(process.kt4CaloJets*process.kt4sub*process.kt4corr*process.kt4clean*process.kt4match*process.kt4parton*process.kt4patJets)



### AK5, PF  ###

process.ak5PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak5PFsub.src = 'ak5PFJets'
process.ak5PFsub.rhoTag = 'kt4PFJets'
process.ak5PFsub.jetType    = 'PFJet'

process.ak5PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak5PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
                                                  )
process.ak5PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.ak5PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ak5PFsub"),
                                                     matched = cms.InputTag("ak5PFclean"))
process.ak5PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ak5PFsub"))
process.ak5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak5PFsub"),
                                             genJetMatch = cms.InputTag("ak5PFmatch"),
                                             genPartonMatch= cms.InputTag("ak5PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFcorr")))
process.ak5PFpatSequence = cms.Sequence(process.ak5PFJets*process.ak5PFsub*process.ak5PFcorr*process.ak5PFclean*process.ak5PFmatch*process.ak5PFparton*process.ak5PFpatJets)




### KT4, PF  ###

process.kt4PFsub = process.heavyIonL1SubtractedJets.clone()
process.kt4PFsub.src = 'kt4PFJets'
process.kt4PFsub.rhoTag = 'kt4PFJets'
process.kt4PFsub.jetType    = 'PFJet'

process.kt4PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("kt4PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('KT4PF')
                                                  )
process.kt4PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt4HiGenJets'))
process.kt4PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("kt4PFsub"),
                                                     matched = cms.InputTag("kt4PFclean"))
process.kt4PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("kt4PFsub"))
process.kt4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("kt4PFsub"),
                                             genJetMatch = cms.InputTag("kt4PFmatch"),
                                             genPartonMatch= cms.InputTag("kt4PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt4PFcorr")))
process.kt4PFpatSequence = cms.Sequence(process.kt4PFJets*process.kt4PFsub*process.kt4PFcorr*process.kt4PFclean*process.kt4PFmatch*process.kt4PFparton*process.kt4PFpatJets)


### IC5, PF  ###

process.ic5PFsub = process.heavyIonL1SubtractedJets.clone()
process.ic5PFsub.src = 'ic5PFJets'
process.ic5PFsub.rhoTag = 'kt4PFJets'
process.ic5PFsub.jetType    = 'PFJet'

process.ic5PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ic5PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5PF')
                                                  )
process.ic5PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone5HiGenJets'))
process.ic5PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ic5PFsub"),
                                                     matched = cms.InputTag("ic5PFclean"))
process.ic5PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ic5PFsub"))
process.ic5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ic5PFsub"),
                                             genJetMatch = cms.InputTag("ic5PFmatch"),
                                             genPartonMatch= cms.InputTag("ic5PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic5PFcorr")))
process.ic5PFpatSequence = cms.Sequence(process.ic5PFJets*process.ic5PFsub*process.ic5PFcorr*process.ic5PFclean*process.ic5PFmatch*process.ic5PFparton*process.ic5PFpatJets)




#begin small cones

process.ak3CaloJets = process.ak5CaloJets.clone()
process.ak3CaloJets.rParam       = cms.double(0.3)
process.ak3CaloJets.radiusPU = 0.3

process.ak3sub = process.heavyIonL1SubtractedJets.clone()
process.ak3sub.src = 'ak3CaloJets'
process.ak3sub.rhoTag = 'kt4CaloJets'

process.ak3corr = process.patJetCorrFactors.clone(src = cms.InputTag("ak3sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5Calo')
                                                  )
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

process.ak4corr = process.patJetCorrFactors.clone(src = cms.InputTag("ak4sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5Calo')
                                                  )
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

process.kt3corr = process.patJetCorrFactors.clone(src = cms.InputTag("kt3sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('KT4Calo')
                                                  )
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

process.ic3corr = process.patJetCorrFactors.clone(src = cms.InputTag("ic3sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
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

process.ic4corr = process.patJetCorrFactors.clone(src = cms.InputTag("ic4sub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
process.ic4clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone4HiGenJets'))
process.ic4match = process.patJetGenJetMatch.clone(src = cms.InputTag("ic4sub"),
                                                   matched = cms.InputTag("ic4clean"))
process.ic4parton = process.patJetPartonMatch.clone(src = cms.InputTag("ic4sub"))
process.ic4patJets = process.patJets.clone(jetSource  = cms.InputTag("ic4sub"),
                                           genJetMatch = cms.InputTag("ic4match"),
                                           genPartonMatch= cms.InputTag("ic4parton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic4corr")))
process.ic4patSequence = cms.Sequence(process.ic4CaloJets*process.ic4sub*process.ic4corr*process.ic4clean*process.ic4match*process.ic4parton*process.ic4patJets)


process.ak3PFJets = process.ak5PFJets.clone()
process.ak3PFJets.rParam       = cms.double(0.3)
process.ak3PFJets.radiusPU = 0.3

process.ak3PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak3PFsub.src = 'ak3PFJets'
process.ak3PFsub.rhoTag = 'kt4PFJets'
process.ak3PFsub.jetType    = 'PFJet'

process.ak3PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak3PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
                                                  )
process.ak3PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiGenJets'))
process.ak3PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ak3PFsub"),
                                                   matched = cms.InputTag("ak3PFclean"))
process.ak3PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ak3PFsub"))
process.ak3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak3PFsub"),
                                           genJetMatch = cms.InputTag("ak3PFmatch"),
                                           genPartonMatch= cms.InputTag("ak3PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3PFcorr")))
process.ak3PFpatSequence = cms.Sequence(process.ak3PFJets*process.ak3PFsub*process.ak3PFcorr*process.ak3PFclean*process.ak3PFmatch*process.ak3PFparton*process.ak3PFpatJets)

process.ak4PFJets = process.ak5PFJets.clone()
process.ak4PFJets.rParam       = cms.double(0.4)
process.ak4PFJets.radiusPU = 0.4

process.ak4PFsub = process.heavyIonL1SubtractedJets.clone()
process.ak4PFsub.src = 'ak4PFJets'
process.ak4PFsub.rhoTag = 'kt4PFJets'
process.ak4PFsub.jetType    = 'PFJet'

process.ak4PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak4PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
                                                  )
process.ak4PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiGenJets'))
process.ak4PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ak4PFsub"),
                                                   matched = cms.InputTag("ak4PFclean"))
process.ak4PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ak4PFsub"))
process.ak4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak4PFsub"),
                                           genJetMatch = cms.InputTag("ak4PFmatch"),
                                           genPartonMatch= cms.InputTag("ak4PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4PFcorr")))
process.ak4PFpatSequence = cms.Sequence(process.ak4PFJets*process.ak4PFsub*process.ak4PFcorr*process.ak4PFclean*process.ak4PFmatch*process.ak4PFparton*process.ak4PFpatJets)


process.ic3PFJets = process.ic5PFJets.clone()
process.ic3PFJets.rParam       = cms.double(0.3)
process.ic3PFJets.radiusPU = 0.3

process.ic3PFsub = process.heavyIonL1SubtractedJets.clone()
process.ic3PFsub.src = 'ic3PFJets'
process.ic3PFsub.rhoTag = 'kt4PFJets'
process.ic3PFsub.jetType    = 'PFJet'

process.ic3PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ic3PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5PF')
                                                  )
process.ic3PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone3HiGenJets'))
process.ic3PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ic3PFsub"),
                                                   matched = cms.InputTag("ic3PFclean"))
process.ic3PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ic3PFsub"))
process.ic3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ic3PFsub"),
                                           genJetMatch = cms.InputTag("ic3PFmatch"),
                                           genPartonMatch= cms.InputTag("ic3PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic3PFcorr")))
process.ic3PFpatSequence = cms.Sequence(process.ic3PFJets*process.ic3PFsub*process.ic3PFcorr*process.ic3PFclean*process.ic3PFmatch*process.ic3PFparton*process.ic3PFpatJets)

process.ic4PFJets = process.ic5PFJets.clone()
process.ic4PFJets.rParam       = cms.double(0.4)
process.ic4PFJets.radiusPU = 0.4

process.ic4PFsub = process.heavyIonL1SubtractedJets.clone()
process.ic4PFsub.src = 'ic4PFJets'
process.ic4PFsub.rhoTag = 'ic4PFJets'
process.ic4PFsub.jetType    = 'PFJet'

process.ic4PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ic4PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5PF')
                                                  )
process.ic4PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone4HiGenJets'))
process.ic4PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ic4PFsub"),
                                                   matched = cms.InputTag("ic4PFclean"))
process.ic4PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ic4PFsub"))
process.ic4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ic4PFsub"),
                                           genJetMatch = cms.InputTag("ic4PFmatch"),
                                           genPartonMatch= cms.InputTag("ic4PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ic4PFcorr")))
process.ic4PFpatSequence = cms.Sequence(process.ic4PFJets*process.ic4PFsub*process.ic4PFcorr*process.ic4PFclean*process.ic4PFmatch*process.ic4PFparton*process.ic4PFpatJets)


process.kt3PFJets = process.kt4PFJets.clone()
process.kt3PFJets.rParam       = cms.double(0.3)
process.kt3PFJets.radiusPU = 0.3

process.kt3PFsub = process.heavyIonL1SubtractedJets.clone()
process.kt3PFsub.src = 'kt3PFJets'
process.kt3PFsub.rhoTag = 'kt4PFJets'
process.kt3PFsub.jetType    = 'PFJet'

process.kt3PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("kt3PFsub"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('KT4PF')
                                                  )
process.kt3PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('kt3HiGenJets'))
process.kt3PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("kt3PFsub"),
                                                   matched = cms.InputTag("kt3PFclean"))
process.kt3PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("kt3PFsub"))
process.kt3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("kt3PFsub"),
                                           genJetMatch = cms.InputTag("kt3PFmatch"),
                                           genPartonMatch= cms.InputTag("kt3PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("kt3PFcorr")))
process.kt3PFpatSequence = cms.Sequence(process.kt3PFJets*process.kt3PFsub*process.kt3PFcorr*process.kt3PFclean*process.kt3PFmatch*process.kt3PFparton*process.kt3PFpatJets)

# Now algorithms with Russian-style PU subtraction:  icPu5Calo, icPu5PF and jpticPu5

process.icPu5corr = process.patJetCorrFactors.clone(src = cms.InputTag("iterativeConePu5CaloJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
process.icPu5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone5HiGenJets'))
process.icPu5match = process.patJetGenJetMatch.clone(src = cms.InputTag("iterativeConePu5CaloJets"),
                                                     matched = cms.InputTag("icPu5clean"))
process.icPu5parton = process.patJetPartonMatch.clone(src = cms.InputTag("iterativeConePu5CaloJets"))
process.icPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("iterativeConePu5CaloJets"),
                                             genJetMatch = cms.InputTag("icPu5match"),
                                             genPartonMatch= cms.InputTag("icPu5parton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5corr")))
process.icPu5patSequence = cms.Sequence(process.iterativeConePu5CaloJets*process.icPu5corr*process.icPu5clean*process.icPu5match*process.icPu5parton*process.icPu5patJets)

# PF using a grid of pseudo-towers to emulate Russian-style subtraction

process.load("RecoHI.HiJetAlgos.ParticleTowerProducer_cff")
process.PFTowers = process.particleTowerProducer.clone()

### IC5, PF  ###
process.icPu5PFJets = process.ic5PFJets.clone()
process.icPu5PFJets.src = 'PFTowers'
process.icPu5PFJets.jetType = 'BasicJet'
process.icPu5PFJets.doPUOffsetCorr = cms.bool(True)
process.icPu5PFJets.sumRecHits = cms.bool(False)

process.icPu5PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("icPu5PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5PF')
)
process.icPu5PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone5HiGenJets'))
process.icPu5PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("icPu5PFJets"),
                                                     matched = cms.InputTag("icPu5PFclean"))
process.icPu5PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("icPu5PFJets"))
process.icPu5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("icPu5PFJets"),
                                             genJetMatch = cms.InputTag("icPu5PFmatch"),
                                             genPartonMatch= cms.InputTag("icPu5PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5PFcorr")))
process.icPu5PFpatSequence = cms.Sequence(process.PFTowers*process.icPu5PFJets*process.icPu5PFcorr*process.icPu5PFclean*process.icPu5PFmatch*process.icPu5PFparton*process.icPu5PFpatJets)


# JPT

process.jpticPu5corr = process.patJetCorrFactors.clone(src = cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
process.jpticPu5clean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('iterativeCone5HiGenJets'))
process.jpticPu5match = process.patJetGenJetMatch.clone(src = cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                        matched = cms.InputTag("jpticPu5clean"))
process.jpticPu5parton = process.patJetPartonMatch.clone(src = cms.InputTag("JetPlusTrackZSPCorJetIconePu5"))
process.jpticPu5patJets = process.patJets.clone(jetSource  =cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                genJetMatch = cms.InputTag("jpticPu5match"),
                                                genPartonMatch= cms.InputTag("jpticPu5parton"),                                       
                                                jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jpticPu5corr")))
process.icPu5JPTpatSequence = cms.Sequence(process.recoJPTJetsHIC*process.jpticPu5corr*process.jpticPu5clean*process.jpticPu5match*process.jpticPu5parton*process.jpticPu5patJets)



process.runAllJets = cms.Sequence(
    process.kt4patSequence +
    process.ic5patSequence +
    #process.ak5patSequence +
    process.kt4PFpatSequence +
    #process.ak5PFpatSequence +
    process.ic5PFpatSequence +
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
    process.icPu5patSequence +
    process.icPu5PFpatSequence #+
    #process.icPu5JPTpatSequence 
)

# jet analysis trees
process.load("MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_cff")
process.icPu5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.icPu5PFJetAnalyzer.jetTag = 'icPu5PFpatJets'


process.inclusiveJetAnalyzerSequence = cms.Sequence(                  process.inclusiveJetAnalyzer
                                                                      #*process.ic3JetAnalyzer
                                                                      #*process.ic4JetAnalyzer
                                                                      *process.ic5JetAnalyzer
                                                                      #*process.ak3JetAnalyzer
                                                                      #*process.ak4JetAnalyzer
                                                                      #*process.ak5JetAnalyzer
                                                                      #*process.kt4JetAnalyzer
                                                                      #*process.ic3PFJetAnalyzer
                                                                      #*process.ic4PFJetAnalyzer
                                                                      *process.ic5PFJetAnalyzer
                                                                      #*process.ak3PFJetAnalyzer
                                                                      #*process.ak4PFJetAnalyzer
                                                                      #*process.ak5PFJetAnalyzer
                                                                      #*process.kt4PFJetAnalyzer
                                                                      #*process.icPu5JPTJetAnalyzer
                                                                      *process.icPu5PFJetAnalyzer
                                                                      )

process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_cff")
#process.PFJetAnalyzer.hasSimInfo = cms.untracked.bool(True)
process.PFJetAnalyzer.jetTag3 = cms.InputTag("icPu5PFpatJets")
process.PFJetAnalyzer.recoJetTag3 = cms.InputTag("icPu5PFJets")

# track efficiency anlayzer
#process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")
#for tree output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("JetAnalysisTTrees_hiGoodMergedTracks_v1.root"))



# put it all together



process.path = cms.Path(
    process.makeCentralityTableTFile*
    process.hiTrackReco*
    process.hiGoodTracksSelection*
    process.conformalPixelTrackReco *
    process.hiGoodMergedTracks *
    process.HiParticleFlowRecoNoJets*
    process.hiExtra*
    process.hiGen*
    process.runAllJets*
    process.inclusiveJetAnalyzerSequence#*
    #process.PFJetAnalyzerSequence
    #process.hitrkEffAna*
    #process.hipxltrkEffAna
    )

#process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")
#
#process.output = cms.OutputModule("PoolOutputModule",
#                                  process.jetTrkSkimContentMC,
#                                  fileName = cms.untracked.string("/tmp/mnguyen/PAT.root")
#                                  )

#process.output.outputCommands.extend(["drop *_towerMaker_*_*"])
#process.output.outputCommands.extend(["keep *_PFTowers_*_*"])
#process.output.outputCommands.extend(["keep *_hiGoodMergedTracks_*_*"])
#process.output.outputCommands.extend(["keep *_particleFlow_*_*"])
#process.output.outputCommands.extend(["keep *_mergedtruth_*_*"])
#process.output.outputCommands.extend(["keep double*_*PF*_*_*"])
#process.output.outputCommands.extend(["keep *_heavyIon_*_*"])
#process.output.outputCommands.extend(["keep *_generator_*_*"])
#process.output.outputCommands.extend(["keep *_hiSignal_*_*"])
#process.output.outputCommands.extend(["keep *_genParticles_*_*"])
#process.output.outputCommands.extend(["keep *_hiGenParticles_*_*"])
#process.output.outputCommands.extend(["keep *_TriggerResults_*_*"])
#process.output.outputCommands.extend(["keep *_heavyIon_*_*"])
## reco jets
#process.output.outputCommands.extend(["keep recoCaloJets_*_*_*"])
#process.output.outputCommands.extend(["keep recoPFJets_*_*_*"])
##particle flow
#process.output.outputCommands.extend(["keep recoPFClusters_*_*_*"])
#process.output.outputCommands.extend(["keep recoPFRecHits_*_*_*"])
##fast jet pf stuff
#process.output.outputCommands.extend(["keep doubles_*PF*_*_*"])
##calorimeter stuff
#process.output.outputCommands.extend(["keep *_towerMaker_*_*"])
#process.output.outputCommands.extend(["keep *_caloTowers_*_*"])
#process.output.outputCommands.extend(["keep *_hcalnoise_*_*"])
#process.output.outputCommands.extend(["keep *_hbhereco_*_*"])
#process.output.outputCommands.extend(["keep *_horeco_*_*"])
#process.output.outputCommands.extend(["keep *_hfreco_*_*"])
#process.output.outputCommands.extend(["keep *_ecalRecHit_*_*"])
##JPT
#process.output.outputCommands.extend(["keep *_jptic*_*_*"])
#process.output.outputCommands.extend(["keep *_recoJPT*_*_*"])
#process.output.outputCommands.extend(["keep *_JetPlusTrack*_*_*"])
#
#
#
#
#process.out_step = cms.EndPath(process.output)

