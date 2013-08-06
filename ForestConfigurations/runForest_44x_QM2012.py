import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.register ('randomNumber',
                1,
                ivars.multiplicity.singleton,
                ivars.varType.int,
                "Random Seed")

ivars.randomNumber = 1
#ivars.inputFiles = "file:/net/hisrv0001/home/icali/hadoop/Pythia/Z2/ppDijet50/reco_v0/set2_random70000_HydjetDrum_362.root"
ivars.inputFiles = "file:/mnt/hadoop/cms/store/user/yetkin/MC_Production/Pythia80_HydjetDrum_mix01/RECO/set1_random30000_HydjetDrum_12.root"

ivars.outputFile = './mixed_440_B_02.root'

ivars.parseArguments()


import FWCore.ParameterSet.Config as cms

process = cms.Process('hiForestAna2011')

process.options = cms.untracked.PSet(
   # wantSummary = cms.untracked.bool(True)
)


doElectrons = False
useIterTracking = False
doRegitForBjets = False

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(ivars.inputFiles))

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1))


#####################################################################################
# Load some general stuff
#####################################################################################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
#process.load('MitHig.PixelTrackletAnalyzer.pixelHitAnalyzer_cfi')

# Data Global Tag 44x 
#process.GlobalTag.globaltag = 'GR_P_V27::All'

# MC Global Tag 44x 
process.GlobalTag.globaltag = 'STARTHI44_V7::All'

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
	nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
	centralitySrc = cms.InputTag("hiCentrality")
	)

process.hiCentrality.pixelBarrelOnly = False
process.load("CmsHi.JetAnalysis.RandomCones_cff")

process.RandomNumberGeneratorService.generator.initialSeed = ivars.randomNumber 
process.RandomNumberGeneratorService.akPu3PFConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.icPu5CaloConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.akPu5PFConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.akPu3CaloConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.akPu5CaloConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.multiPhotonAnalyzer = process.RandomNumberGeneratorService.generator.clone()

# EcalSeverityLevel ES Producer
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load("RecoEcal.EgammaCoreTools.EcalNextToDeadChannelESProducer_cff")


# Pixel Tracks
process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')


#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string(ivars.outputFile))

#####################################################################################
# Additional Reconstruction and Analysis
#####################################################################################

# MET: Calorimeter based MET
process.load("RecoMET.METProducers.CaloMET_cfi") 

process.load("CmsHi.JetAnalysis.hcalNoise_cff")

# Define Analysis sequencues
process.load('CmsHi.JetAnalysis.EventSelection_cff')
process.load('CmsHi.JetAnalysis.ExtraGenReco_cff')
#process.load('CmsHi.JetAnalysis.ExtraTrackReco_cff')
process.load('CmsHi.JetAnalysis.ExtraPfReco_cff')
process.load('CmsHi.JetAnalysis.ExtraJetReco_cff')
process.load('CmsHi.JetAnalysis.ExtraEGammaReco_cff')
process.load('CmsHi.JetAnalysis.PatAna_cff')
process.load('CmsHi.JetAnalysis.JetAnalyzers_cff')
process.load('CmsHi.JetAnalysis.TrkAnalyzers_cff')
process.load('CmsHi.JetAnalysis.EGammaAnalyzers_cff')
process.load("MitHig.PixelTrackletAnalyzer.trackAnalyzer_cff")

process.anaTrack.trackPtMin = 0.5
process.anaTrack.useQuality = True
#print "Temporarily commented out sim tracks "
process.anaTrack.doSimTrack = True
process.anaTrack.simTrackPtMin = 2
process.anaTrack.doPFMatching = True
process.anaTrack.pfCandSrc = cms.InputTag("particleFlowTmp")
process.anaTrack.trackSrc = cms.InputTag("hiSelectedTracks")

process.pixelTrack = process.anaTrack.clone(trackSrc = cms.InputTag("hiConformalPixelTracks"))
process.pixelTrack.useQuality = False
process.pixelTrack.trackPtMin = 0.5

process.load("MitHig.PixelTrackletAnalyzer.METAnalyzer_cff")
process.load("CmsHi.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.interestingTrackEcalDetIds.TrackCollection = cms.InputTag("hiSelectedTracks")

#process.load("HiMuonAlgos.HLTMuTree.hltMuTree_cfi")

process.genpana = cms.EDAnalyzer("GenParticleCounter",
                                 src = cms.untracked.string("hiGenParticles"),
                                 doCentrality = cms.untracked.bool(False),
                                 VertexProducer = cms.untracked.string("hiSelectedVertex")
                                 )


# Muons 
process.load("MuTrig.HLTMuTree.hltMuTree_cfi")
process.muonTree = process.hltMuTree.clone()
process.muonTree.doGen = cms.untracked.bool(True)

# Event tree
process.load("CmsHi/HiHLTAlgos.hievtanalyzer_cfi")
# Not working for the moment..
#process.hiEvtAnalyzer.doMC = cms.bool(True)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(True)

process.hiGenParticles.srcVector = cms.vstring('hiSignal','generator')
process.icPu5JetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.akPu3PFJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.akPu5PFJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.akPu3CaloJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.akPu5CaloJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.multiPhotonAnalyzer.GenEventScale = cms.InputTag("hiSignal")

process.multiPhotonAnalyzer.HepMCProducer = cms.InputTag("hiSignal")

process.icPu5JetAnalyzer.hltTrgResults = cms.untracked.string('TriggerResults::RECO')
process.akPu3PFJetAnalyzer.hltTrgResults = cms.untracked.string('TriggerResults::RECO')
process.icPu5JetAnalyzer.isMC   = cms.untracked.bool(True)
process.akPu3PFJetAnalyzer.isMC = cms.untracked.bool(True)

process.iterativeCone5CaloJets = process.iterativeConePu5CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 1)
process.ak5PFJets = process.akPu5PFJets.clone(doPUOffsetCorr = False, jetPtMin = 1)
process.ak5CaloJets = process.akPu5CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 1)
process.ak3PFJets = process.akPu3PFJets.clone(doPUOffsetCorr = False, jetPtMin = 1)
process.ak3CaloJets = process.akPu3CaloJets.clone(doPUOffsetCorr = False, jetPtMin = 1)

process.ic5JetAnalyzer = process.icPu5JetAnalyzer.clone(jetTag = cms.InputTag("iterativeCone5CaloJets"))
process.ak3PFJetAnalyzer = process.akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak3PFJets"))
process.ak5PFJetAnalyzer = process.akPu5PFJetAnalyzer.clone(jetTag = cms.InputTag("ak5PFJets"))
process.ak3CaloJetAnalyzer = process.akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak3CaloJets"))
process.ak5CaloJetAnalyzer = process.akPu5CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak5CaloJets"))


#Commented by Yen-Jie
#process.hiPixelAdaptiveVertex.useBeamConstraint = False

process.HiGenParticleAna = cms.EDAnalyzer("HiGenAnalyzer")
process.HiGenParticleAna.src= cms.untracked.InputTag("hiGenParticles")
process.HiGenParticleAna.chargedOnly = cms.untracked.bool(False)

process.load("RecoHI.HiMuonAlgos.HiRecoMuon_cff")
process.muons.JetExtractorPSet.JetCollectionLabel = cms.InputTag("iterativeConePu5CaloJets")

process.hiSelectedTrackHighPurity = cms.EDFilter("TrackSelector",
                                                 src = cms.InputTag("hiSelectedTracks"),
                                                 cut = cms.string(
    'quality("highPurity")')
                                                 )

process.particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)

#process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")



process.cutsTPForFak.tracks = cms.untracked.InputTag('TrackingParticles')

#  Following is the reco before adding b-tagging -Matt
'''
if doElectrons:
    process.reco_extra        = cms.Path( process.siPixelRecHits * process.siStripMatchedRecHits *
                                          process.hiPixel3PrimTracks *
                                          process.hiPixelTrackSeeds *
                                          process.hiSelectedTrackHighPurity *
                                          process.electronGsfTrackingHi *
                                          process.hiElectronSequence *
                                          process.HiParticleFlowReco *
                                          process.iterativeConePu5CaloJets *
                                          process.PFTowers
                                          )
    
else:
    process.pfTrack.GsfTracksInEvents = cms.bool(False)
    process.HiParticleFlowReco.remove(process.electronsWithPresel)
    process.HiParticleFlowReco.remove(process.electronsCiCLoose)
    process.particleFlowTmp.usePFElectrons = cms.bool(False)
    process.particleFlowTmp.useEGammaElectrons = cms.bool(False)
    
    process.reco_extra        = cms.Path( process.hiSelectedTrackHighPurity *
                                          process.HiParticleFlowLocalReco *
                                          process.HiParticleFlowReco *
                                          process.iterativeConePu5CaloJets *
                                          process.PFTowers
                                          )
'''

###########################################
# Here it is after including b-tagging -Matt

process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
process.hiTrackReco = cms.Sequence(process.rechits *
                                   process.heavyIonTracking *
                                   process.hiSelectedTrackHighPurity *
                                   process.conformalPixelTrackReco)

# fixed necessary for muons in HI -Matt
process.particleFlowBlock.RecMuons = 'muons'
process.particleFlowTmp.postMuonCleaning = False
# should be set elsewhere, but just to make double sure...
process.particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)
# seed PF with high purity tracks only
process.pfTrack.UseQuality = True 
process.pfTrack.TrackQuality = cms.string('highPurity')


if doElectrons:
    process.reco_extra =  cms.Path(
        process.hiTrackReco
        *process.muonRecoPbPb
        *process.electronGsfTrackingHi        
        *process.HiParticleFlowLocalReco
        *process.gsfEcalDrivenElectronSequence
        *process.HiParticleFlowReco
        *process.iterativeConePu5CaloJets
        *process.PFTowers
        )        
else:
    process.pfTrack.GsfTracksInEvents = cms.bool(False)
    process.HiParticleFlowReco.remove(process.electronsWithPresel)
    process.HiParticleFlowReco.remove(process.electronsCiCLoose)
    process.particleFlowTmp.usePFElectrons = cms.bool(False)
    process.particleFlowTmp.useEGammaElectrons = cms.bool(False)         

    process.reco_extra =  cms.Path(
        process.hiTrackReco
        *process.muonRecoPbPb
        *process.HiParticleFlowLocalReco
        *process.HiParticleFlowReco
        *process.iterativeConePu5CaloJets
        *process.PFTowers
        )
    

if useIterTracking:
    # load the iterative tacks
    process.load("RecoHI.HiTracking.hiIterTracking_cff")
    process.heavyIonTracking *= process.hiIterTracking
    
    # seed the muons with iterative tracks
    process.globalMuons.TrackerCollectionLabel = "hiGeneralTracks"
    process.muons.TrackExtractorPSet.inputTrackCollection = "hiGeneralTracks"
    process.muons.inputCollectionLabels = ["hiGeneralTracks", "globalMuons", "standAloneMuons:UpdatedAtVtx", "tevMuons:firstHit", "tevMuons:picky", "tevMuons:dyt"]

    # set track collection to iterative tracking
    process.pfTrack.TkColList = cms.VInputTag("hiGeneralTracks")


# End modifications to reco sequence -Matt
##########################################




process.reco_extra_jet    = cms.Path( process.iterativeConePu5CaloJets *
                                      process.iterativeCone5CaloJets *
                                      process.ak3PFJets * process.ak5PFJets *
                                      process.ak3CaloJets * process.ak5CaloJets *     
                                      process.akPu3PFJets * process.akPu5PFJets *
                                      process.akPu3CaloJets * process.akPu5CaloJets *
                                      process.photon_extra_reco)
process.gen_step          = cms.Path( process.hiGenParticles * process.hiGenParticlesForJets * process.genPartons * process.hiPartons * process.hiRecoGenJets)

###########################################
# Do regit and b-tagging

if doRegitForBjets:
    

    process.akPu3PFSelectedJets = cms.EDFilter("CandViewSelector",
                                               src = cms.InputTag("akPu3PFJets"),
                                               cut = cms.string("pt > 40. && eta > -2. && eta < 2")
                                               )
    
    
    process.load("RecoHI.HiTracking.hiRegitTracking_cff")
    
    if useIterTracking == False:
        process.hiGeneralTrackFilter.recTracks = cms.InputTag("hiSelectedTracks") # uses hiGeneralTracks by default
        del process.hiRegitInitialStepClusters.oldClusterRemovalInfo

        
    # Region optimized for b-jets
    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 0.02
    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 0.02
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 0.02
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 0.02
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 0.015
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 0.015
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 1.5
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 1.5
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.originRadius = 1.0
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.originHalfLength = 1.0
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.originRadius = 0.5
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.originHalfLength = 0.5

    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.deltaPhiRegion = 0.3
    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.deltaEtaRegion = 0.3
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.deltaPhiRegion = 0.3
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.deltaEtaRegion = 0.3
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.deltaPhiRegion = 0.3
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.deltaEtaRegion = 0.3
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.deltaPhiRegion = 0.3
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.deltaEtaRegion = 0.3
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.deltaPhiRegion = 0.3
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.deltaEtaRegion = 0.3
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.deltaPhiRegion = 0.3
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.deltaEtaRegion = 0.3

    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")


    # merge the regit with the global tracking
    process.load("RecoHI.HiTracking.MergeRegit_cff")

    if useIterTracking:
        process.hiGeneralAndRegitTracks.TrackProducer1 = 'hiGeneralTracks'
    else:
        process.hiGeneralAndRegitTracks.TrackProducer1 = 'hiSelectedTracks'

    # redo the muons, too to get the displaced muons in jets
    process.regGlobalMuons = process.globalMuons.clone(
        TrackerCollectionLabel = "hiGeneralAndRegitTracks"
        )
    process.regGlbTrackQual = process.glbTrackQual.clone(
        InputCollection = "regGlobalMuons",
        InputLinksCollection = "regGlobalMuons"
        )
    process.regMuons = process.muons.clone()
    process.regMuons.TrackExtractorPSet.inputTrackCollection = "hiGeneralAndRegitTracks"
    process.regMuons.globalTrackQualityInputTag = "regGlbTrackQual"
    process.regMuons.inputCollectionLabels = cms.VInputTag("hiGeneralAndRegitTracks", "regGlobalMuons", "standAloneMuons:UpdatedAtVtx", "tevMuons:firstHit", "tevMuons:picky", "tevMuons:dyt")
    
    process.regMuonReco = cms.Sequence(
        process.regGlobalMuons*
        process.regGlbTrackQual*
        process.regMuons
        )
    
    # re-run particle flow to do the calo-matching
    process.pfRegTrack = process.pfTrack.clone(TkColList = cms.VInputTag("hiGeneralAndRegitTracks"))
    process.pfRegBlock = process.particleFlowBlock.clone()
    process.regParticleFlow = process.particleFlowTmp.clone()
        
    process.pfRegTrack.GsfTracksInEvents = False
    process.regParticleFlow.usePFElectrons = False
    process.regParticleFlow.muons = "regMuons"
        
    process.hiRegPF =  cms.Sequence(
        process.pfRegTrack
        *process.pfRegBlock
        *process.regParticleFlow
    )

    process.load("RecoHI.HiTracking.HICaloCompatibleTracks_cff")
    process.hiGeneralAndRegitCaloMatchedTracks = process.hiCaloCompatibleTracks.clone(    
        srcTracks = 'hiGeneralAndRegitTracks',
        srcPFCands = 'regParticleFlow'
        )

        
    process.hiCaloMatchFilteredTracks = cms.EDFilter("TrackSelector",
                                                     src = cms.InputTag("hiGeneralAndRegitCaloMatchedTracks"),
                                                     cut = cms.string(
        'quality("highPuritySetWithPV") && pt > 1')                                                                                            
                                                     )
    
    
    process.regionalTracking = cms.Path(
        process.akPu3PFSelectedJets *
        process.hiRegitTracking *
        process.hiGeneralAndRegitTracks *
        process.regMuonReco *
        process.hiRegPF *
        process.hiGeneralAndRegitCaloMatchedTracks *
        process.hiCaloMatchFilteredTracks
        )
    
    
    process.icPu5SecondaryVertexTagInfos.trackSelection.qualityClass = 'any'
    process.akPu3PFSecondaryVertexTagInfos.trackSelection.qualityClass = 'any'
    process.icPu5JetTracksAssociatorAtVertex.tracks = cms.InputTag("hiCaloMatchFilteredTracks")
    process.akPu3PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("hiCaloMatchFilteredTracks")

    process.akPu3PFpatJets.addAssociatedTracks = True
    process.akPu3PFpatJets.addTagInfos = True
    process.akPu3PFpatJets.addBTagInfo         = True
    process.akPu3PFpatJets.addDiscriminators   = True
    process.akPu3PFpatJets.getJetMCFlavour     = True


    process.akPu3PFpatJets.tagInfoSources = cms.VInputTag(
        cms.InputTag("akPu3PFImpactParameterTagInfos"),
        cms.InputTag("akPu3PFGhostTrackVertexTagInfos"),
        cms.InputTag("akPu3PFSecondaryVertexTagInfos"),
        cms.InputTag("akPu3PFSoftMuonTagInfos"),
        )


# end btagging mods, expect for pat replace just below -Matt
##########################################




process.pat_step          = cms.Path(process.icPu5patSequence +
                                     process.akPu3PFpatSequence + process.akPu5PFpatSequence +
                                     process.akPu3patSequence + process.akPu5patSequence +
                                     process.makeHeavyIonPhotons)

# change to btagging pat sequence
if doRegitForBjets:
    process.pat_step.replace(process.akPu3PFpatSequence,process.akPu3PFpatSequence_withBtagging)
# Note that for data you need: akPu3PFpatSequence_withBtagging_data
    
process.pat_step.remove(process.interestingTrackEcalDetIds)
process.photonMatch.matched = cms.InputTag("hiGenParticles")
#process.pat_step.remove(process.photonMatch)
#+ process.patPhotons)

process.patPhotons.addPhotonID = cms.bool(False)
#process.makeHeavyIonPhotons)
process.extrapatstep = cms.Path(process.selectedPatPhotons)

process.multiPhotonAnalyzer.GammaEtaMax = cms.untracked.double(100)
process.multiPhotonAnalyzer.GammaPtMin = cms.untracked.double(0)
process.ana_step          = cms.Path( process.genpana +
                                      process.hcalNoise +
                                      process.jetana_seq +                                      
                                      process.multiPhotonAnalyzer +
                                      process.HiGenParticleAna +
#                                      process.cutsTPForFak +
#                                      process.cutsTPForEff +
                                      process.anaTrack + process.pixelTrack +
                                      process.pfcandAnalyzer +
                                      process.met * process.anaMET +
				      process.muonTree +
				      process.hiEvtAnalyzer +
                                      process.randomCones +
                                      process.hltMuTree
                                      )


process.phltJetHI = cms.Path( process.hltJetHI )
process.pcollisionEventSelection = cms.Path(process.collisionEventSelection)
process.pHBHENoiseFilter = cms.Path( process.HBHENoiseFilter )
process.phiEcalRecHitSpikeFilter = cms.Path(process.hiEcalRecHitSpikeFilter )
#process.ppreTrgTest = cms.Path(process.preTrgTest )
#process.pminBiasBscFilter = cms.Path(process.minBiasBscFilter )
#process.ppostTrgTest = cms.Path(process.postTrgTest )
#process.phfCoincFilter = cms.Path(process.hfCoincFilter )
#process.ppurityFractionFilter = cms.Path(process.purityFractionFilter )

# Customization
from CmsHi.JetAnalysis.customise_cfi import *
setPhotonObject(process,"cleanPhotons")

process.load('L1Trigger.Configuration.L1Extra_cff')
process.load('CmsHi.HiHLTAlgos.hltanalysis_cff')

process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","RECO")
process.hltAna = cms.Path(process.hltanalysis)
process.pAna = cms.EndPath(process.skimanalysis)
process.reco_extra*=process.L1Extra


process.load('CmsHi.JetAnalysis.rechitanalyzer_cfi')
  
process.rechitanalyzer.HBHETreePtMin = cms.untracked.double(0.5)
process.rechitanalyzer.HFTreePtMin = cms.untracked.double(0.5)
process.rechitanalyzer.EBTreePtMin = cms.untracked.double(0.5)
process.rechitanalyzer.EETreePtMin = cms.untracked.double(0.5)
process.rechitanalyzer.TowerTreePtMin = cms.untracked.double(0.5)
process.rechitanalyzer.doHF = cms.untracked.bool(True)

process.pfTowers = process.rechitanalyzer.clone(
    doEcal  = cms.untracked.bool(False),
    doHcal  = cms.untracked.bool(False),
    hasVtx  = cms.untracked.bool(False),
    doFastJet = cms.untracked.bool(False),
    towersSrc = cms.untracked.InputTag("PFTowers")
    )

process.rechitAna = cms.Path(process.rechitanalyzer+process.pfTowers)

########### random number seed

#####################################################################################
# Edm Output
#####################################################################################

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("output.root")
#                               )
#process.save = cms.EndPath(process.out)
