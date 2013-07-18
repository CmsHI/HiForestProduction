import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')

#process.Timing = cms.Service("Timing")

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'/store/relval/CMSSW_3_10_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/START310_V3-v1/0050/54ED7C66-680E-E011-BD9F-001A92810ADE.root'
    #'/store/relval/CMSSW_3_9_9/RelValHydjetQ_B0_2760GeV/GEN-SIM-RECO/START39_V7HI-v1/0001/0618A00D-E23D-E011-A7FE-001A9281173E.root'                                                              
    '/store/user/davidlw/Hydjet_Bass_MinBias_2760GeV/Pyquen_UnquenchedDiJet_Pt50_GEN-SIM-RECO_393_v2/1243c1b8707a4e7eb28eae64e1474920/hiReco_RAW2DIGI_RECO_9_1_ZVb.root'
    )
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START39_V9::All'

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


#Track Reco
process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking)


# good track selection
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")
# merge with pixel tracks
#process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
#process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')

# for PF
process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")


# particle-flow stuff

process.HiParticleFlowRecoNoJets = cms.Sequence(
    process.particleFlowCluster
    * process.trackerDrivenElectronSeeds
    * process.particleFlowReco
    )
process.trackerDrivenElectronSeeds.TkColList = cms.VInputTag("hiGoodTracks")

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


process.ak3PFJets = process.ak5PFJets.clone()
process.ak3PFJets.rParam       = cms.double(0.3)


process.ak4PFJets = process.ak5PFJets.clone()
process.ak4PFJets.rParam       = cms.double(0.4)

    
process.iterativeConePu5CaloJets.puPtMin = cms.double(15.0)
process.ak5PFJets.puPtMin = cms.double(30.0)
process.ak4PFJets.puPtMin = cms.double(25.0)
process.ak3PFJets.puPtMin = cms.double(20.0)


# PF using a grid of pseudo-towers to emulate Russian-style subtraction

process.load("RecoHI.HiJetAlgos.ParticleTowerProducer_cff")
process.PFTowers = process.particleTowerProducer.clone()




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


### AKPU5, PF  ###
process.akPu5PFJets = process.ak5PFJets.clone()
process.akPu5PFJets.src = 'PFTowers'
process.akPu5PFJets.jetType = 'BasicJet'
process.akPu5PFJets.doPUOffsetCorr = cms.bool(True)
process.akPu5PFJets.sumRecHits = cms.bool(False)

process.akPu5PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("akPu5PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
)
process.akPu5PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.akPu5PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("akPu5PFJets"),
                                                     matched = cms.InputTag("akPu5PFclean"))
process.akPu5PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("akPu5PFJets"))
process.akPu5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("akPu5PFJets"),
                                             genJetMatch = cms.InputTag("akPu5PFmatch"),
                                             genPartonMatch= cms.InputTag("akPu5PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu5PFcorr")))
process.akPu5PFpatSequence = cms.Sequence(process.akPu5PFJets*process.akPu5PFcorr*process.akPu5PFclean*process.akPu5PFmatch*process.akPu5PFparton*process.akPu5PFpatJets)



### AKPU4, PF  ###
process.akPu4PFJets = process.ak4PFJets.clone()
process.akPu4PFJets.src = 'PFTowers'
process.akPu4PFJets.jetType = 'BasicJet'
process.akPu4PFJets.doPUOffsetCorr = cms.bool(True)
process.akPu4PFJets.sumRecHits = cms.bool(False)

process.akPu4PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("akPu4PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
)
process.akPu4PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiGenJets'))
process.akPu4PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("akPu4PFJets"),
                                                     matched = cms.InputTag("akPu4PFclean"))
process.akPu4PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("akPu4PFJets"))
process.akPu4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("akPu4PFJets"),
                                             genJetMatch = cms.InputTag("akPu4PFmatch"),
                                             genPartonMatch= cms.InputTag("akPu4PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu4PFcorr")))
process.akPu4PFpatSequence = cms.Sequence(process.akPu4PFJets*process.akPu4PFcorr*process.akPu4PFclean*process.akPu4PFmatch*process.akPu4PFparton*process.akPu4PFpatJets)



### AKPU3, PF  ###
process.akPu3PFJets = process.ak3PFJets.clone()
process.akPu3PFJets.src = 'PFTowers'
process.akPu3PFJets.jetType = 'BasicJet'
process.akPu3PFJets.doPUOffsetCorr = cms.bool(True)
process.akPu3PFJets.sumRecHits = cms.bool(False)

process.akPu3PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("akPu3PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
)
process.akPu3PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiGenJets'))
process.akPu3PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("akPu3PFJets"),
                                                     matched = cms.InputTag("akPu3PFclean"))
process.akPu3PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("akPu3PFJets"))
process.akPu3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("akPu3PFJets"),
                                             genJetMatch = cms.InputTag("akPu3PFmatch"),
                                             genPartonMatch= cms.InputTag("akPu3PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu3PFcorr")))
process.akPu3PFpatSequence = cms.Sequence(process.akPu3PFJets*process.akPu3PFcorr*process.akPu3PFclean*process.akPu3PFmatch*process.akPu3PFparton*process.akPu3PFpatJets)



# JPT


#
process.JetPlusTrackZSPCorJetIconePu5.tagName = cms.vstring('ZSP_CMSSW361_Akt_05_PU0')
process.JetPlusTrackZSPCorJetIconePu5.LeakageMap = "CondFormats/JetMETObjects/data/CMSSW_362_TrackLeakage.txt"
process.JetPlusTrackZSPCorJetIconePu5.ResponseMap = "CondFormats/JetMETObjects/data/CMSSW_362_response.txt"
process.JetPlusTrackZSPCorJetIconePu5.EfficiencyMap = "CondFormats/JetMETObjects/data/CMSSW_362_TrackNonEff.txt"
process.JetPlusTrackZSPCorJetIconePu5.UseEfficiency = True
process.JetPlusTrackZSPCorJetIconePu5.Verbose = cms.bool(False)

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

# set JPT to look at the right track collection:
process.trackExtrapolator.trackSrc = cms.InputTag("hiGoodTracks")
process.JPTiterativeConePu5JetTracksAssociatorAtVertex.tracks = 'hiGoodTracks'
process.JPTiterativeConePu5JetTracksAssociatorAtCaloFace.tracks = 'hiGoodTracks'
process.JetPlusTrackZSPCorJetIconePu5.tracks = 'hiGoodTracks'


process.runAllJets = cms.Sequence(
    process.icPu5patSequence +
    process.akPu5PFpatSequence +
    process.akPu4PFpatSequence +
    process.akPu3PFpatSequence +
    process.icPu5JPTpatSequence 
)

# jet analysis trees
process.load("MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_cff")
#process.inclusiveJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.akPu5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu5PFJetAnalyzer.jetTag = 'akPu5PFpatJets'
process.akPu5PFJetAnalyzer.genjetTag = 'ak5HiGenJets'
process.akPu4PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu4PFJetAnalyzer.jetTag = 'akPu4PFpatJets'
process.akPu4PFJetAnalyzer.genjetTag = 'ak4HiGenJets'
process.akPu3PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu3PFJetAnalyzer.jetTag = 'akPu3PFpatJets'
process.akPu3PFJetAnalyzer.genjetTag = 'ak3HiGenJets'


process.inclusiveJetAnalyzerSequence = cms.Sequence(                  process.inclusiveJetAnalyzer
                                                                      *process.akPu5PFJetAnalyzer
                                                                      *process.akPu4PFJetAnalyzer
                                                                      *process.akPu3PFJetAnalyzer
                                                                      )

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_cff")
process.PFJetAnalyzer.trackTag  = cms.InputTag("hiGoodTracks")
process.PFJetAnalyzer.hasSimInfo = cms.untracked.bool(True)
process.PFJetAnalyzer.jetTag2 = cms.InputTag("akPu5PFpatJets")
process.PFJetAnalyzer.recoJetTag2 = cms.InputTag("akPu5PFJets")
process.PFJetAnalyzer.jetTag3 = cms.InputTag("akPu4PFpatJets")
process.PFJetAnalyzer.recoJetTag3 = cms.InputTag("akPu4PFJets")
process.PFJetAnalyzer.jetTag4 = cms.InputTag("akPu3PFpatJets")
process.PFJetAnalyzer.recoJetTag4 = cms.InputTag("akPu3PFJets")

# track efficiency anlayzer
#process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")
#for tree output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("JetAnalysisTTrees_hiGoodTracks_v1.root"))



# put it all together



process.path = cms.Path(
    process.makeCentralityTableTFile*
    process.hiTrackReco*
    process.hiGoodTracksSelection*
    #process.conformalPixelTrackReco *
    #process.hiGoodMergedTracks *
    process.HiParticleFlowRecoNoJets*
    process.hiExtra*
    process.hiGen*
    process.PFTowers*
    process.runAllJets*
    process.inclusiveJetAnalyzerSequence*
    process.PFJetAnalyzerSequence
    #process.hitrkEffAna*
    #process.hipxltrkEffAna
    )

process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContentMC,
                                  fileName = cms.untracked.string("/tmp/mnguyen/PAT.root")
                                  )

process.output.outputCommands.extend(["drop *_towerMaker_*_*"])
process.output.outputCommands.extend(["keep *_PFTowers_*_*"])
process.output.outputCommands.extend(["keep *_hiGoodMergedTracks_*_*"])
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
#JPT
process.output.outputCommands.extend(["keep *_jptic*_*_*"])
process.output.outputCommands.extend(["keep *_recoJPT*_*_*"])
process.output.outputCommands.extend(["keep *_JetPlusTrack*_*_*"])




process.out_step = cms.EndPath(process.output)

