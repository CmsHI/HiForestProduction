import FWCore.ParameterSet.Config as cms

process = cms.Process('hiForestAna2011')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
    'file:/mnt/hadoop/cms/store/himc/Summer11/Hydjet_Bass_MinBias_2760GeV/GEN-SIM/STARTHI44_V4-v1/0002/4C8726CB-80FA-E011-8692-0002C90B3966.root',
#    'file:/mnt/hadoop/cms/store/user/yinglu/MC_Production/photon30/RECO/edmOut_364.root',
    ))

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(5))


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
# Data Global Tag 44x 
#process.GlobalTag.globaltag = 'GR_P_V27::All'

# MC Global Tag 44x 
process.GlobalTag.globaltag = 'STARTHI44_V7::All'

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFhits"),
	nonDefaultGlauberModel = cms.string("Hydjet_2760GeV"),
	centralitySrc = cms.InputTag("hiCentrality")
	)

process.hiCentrality.pixelBarrelOnly = False
process.load("CmsHi.JetAnalysis.RandomCones_cff")

process.RandomNumberGeneratorService.akPu3PFConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.icPu5CaloConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.akPu5PFConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.akPu3CaloConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.akPu5CaloConesAna = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.multiPhotonAnalyzer = process.RandomNumberGeneratorService.generator.clone()

# EcalSeverityLevel ES Producer
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load("RecoEcal.EgammaCoreTools.EcalNextToDeadChannelESProducer_cff")


#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("HiForest.root"))

#####################################################################################
# Jet energy correction
#####################################################################################

process.jec = cms.ESSource("PoolDBESSource",
	DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
	timetype = cms.string('runnumber'),
	toGet = cms.VPSet(
		cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         tag = cms.string("JetCorrectorParametersCollection_HI_Calo_hiGoodTightTracks_D6T_413_IC5Calo"),
                         label = cms.untracked.string("IC5Calo")),
                
		cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_413_AK3PF'),
                         label = cms.untracked.string("AK3PF")),
                
		cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_413_AK4PF'),
                         label = cms.untracked.string("AK4PF")),
                
		cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_413_AK5PF'),
                         label = cms.untracked.string("AK5PF")
                         ),
	),
                           connect = cms.string("sqlite_file:JEC_HI_PFTower_413patch2_2011_v3.db"),
                           )
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#####################################################################################
# Additional Reconstruction and Analysis
#####################################################################################

# MET: Calorimeter based MET
process.load("RecoMET.METProducers.CaloMET_cfi") 

# Define Analysis sequencues
process.load('CmsHi.JetAnalysis.EventSelection_cff')
process.load('CmsHi.JetAnalysis.ExtraGenReco_cff')
#process.load('CmsHi.JetAnalysis.ExtraTrackReco_cff')
process.load('CmsHi.JetAnalysis.ExtraPfReco_cff')
process.load('CmsHi.JetAnalysis.ExtraJetReco_cff')
process.load('CmsHi.JetAnalysis.ExtraEGammaReco_cff')
process.load('CmsHi.JetAnalysis.PatAna_cff')
process.load('CmsHi.JetAnalysis.JetAnalyzers_cff')
process.load('CmsHi.JetAnalysis.EGammaAnalyzers_cff')
process.load("MitHig.PixelTrackletAnalyzer.trackAnalyzer_cff")
process.anaTrack.trackPtMin = 0
process.anaTrack.useQuality = True
process.anaTrack.doPFMatching = True
process.anaTrack.trackSrc = cms.InputTag("hiSelectedTracks")
process.load("MitHig.PixelTrackletAnalyzer.METAnalyzer_cff")
process.load("CmsHi.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.interestingTrackEcalDetIds.TrackCollection = cms.InputTag("hiSelectedTracks")

# Muons 
process.load("MuTrig.HLTMuTree.hltMuTree_cfi")
process.muonTree = process.hltMuTree.clone()
process.muonTree.doGen = cms.untracked.bool(True)

# Event tree
process.load("CmsHi/HiHLTAlgos.hievtanalyzer_cfi")
# Not working for the moment..
#process.hiEvtAnalyzer.doMC = cms.bool(True)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(True)


# Jet stuff
process.ak5CaloJets = process.akPu5CaloJets.clone(doPUOffsetCorr = False)
process.ak5corr = process.icPu5corr.clone(
        src = cms.InputTag("ak5CaloJets"),
            payload = cms.string('AK5Calo')
            )
process.ak5patJets = process.akPu5PFpatJets.clone(
        jetSource = cms.InputTag("ak5CaloJets"),
            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5corr"))
            )

process.hiGenParticles.srcVector = cms.vstring('generator')
process.icPu5JetAnalyzer.eventInfoTag = cms.InputTag("generator")
process.akPu3PFJetAnalyzer.eventInfoTag = cms.InputTag("generator")
process.multiPhotonAnalyzer.GenEventScale = cms.InputTag("generator")

process.icPu5JetAnalyzer.hltTrgResults = cms.untracked.string('TriggerResults::RECO')
process.akPu3PFJetAnalyzer.hltTrgResults = cms.untracked.string('TriggerResults::RECO')
process.icPu5JetAnalyzer.isMC   = cms.untracked.bool(True)
process.akPu3PFJetAnalyzer.isMC = cms.untracked.bool(True)
process.icPu5JetAnalyzer.useCentrality   = cms.untracked.bool(False) # doesn't fill cent info
process.akPu3PFJetAnalyzer.useCentrality = cms.untracked.bool(False) # doesn't fill cent info


process.ak5CaloJetAnalyzer = process.icPu5JetAnalyzer.clone(
    jetTag = 'ak5patJets',
    genjetTag = 'ak5HiGenJets',
    isMC = True
    )

process.ak3PFJetsX = process.akPu3PFJets.clone(doPUOffsetCorr = False)
process.ak3corrX = process.icPu5corr.clone(
    src = cms.InputTag("ak3PFJetsX"),
    payload = cms.string('AK3PF')
    )
process.ak3patJetsX = process.akPu5PFpatJets.clone(
    jetSource = cms.InputTag("ak3PFJetsX"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3corrX"))
    )
process.ak3PFJetAnalyzer = process.icPu5JetAnalyzer.clone(
    jetTag = 'ak3patJetsX',
    genjetTag = 'ak3HiGenJets',
    isMC = True
    )

process.ak5extra = cms.Sequence(process.ak5CaloJets*process.ak5corr*process.ak5patJets*process.ak5CaloJetAnalyzer)
process.ak3extra = cms.Sequence(process.ak3PFJetsX*process.ak3corrX*process.ak3patJetsX*process.ak3PFJetAnalyzer)

#process.load("edwenger.Skims.EventFilter_cff")
#from edwenger.Skims.customise_cfi import *
#run2760GeVmode(process)

# Filtering
process.hltJetHI.HLTPaths = ['HLT_HIJet35U','HLT_HIPhoton20']
process.hltJetHI.TriggerResultsTag = cms.InputTag("TriggerResults::RECO")

print "Add cleaning to analysis"
process.event_filter_seq = cms.Sequence(
  process.hltJetHI * 
  process.siPixelRecHits *
  process.collisionEventSelection *
  process.HBHENoiseFilter *
  process.hiEcalRecHitSpikeFilter 

#  process.preTrgTest *
#  process.minBiasBscFilter *
#  process.postTrgTest *
#  process.hfCoincFilter *
#  process.purityFractionFilter

  )

#Commented by Yen-Jie
#process.hiPixelAdaptiveVertex.useBeamConstraint = False

process.load("RecoHI.HiMuonAlgos.HiRecoMuon_cff")
process.muons.JetExtractorPSet.JetCollectionLabel = cms.InputTag("iterativeConePu5CaloJets")

process.hiSelectedTrackHighPurity = cms.EDFilter("TrackSelector",
                                                 src = cms.InputTag("hiSelectedTracks"),
                                                 cut = cms.string(
    'quality("highPurity")')
                                                 )

process.particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)
process.reco = cms.Path(process.pdigi * process.SimL1Emulator * process.DigiToRaw * process.RawToDigi * process.reconstructionHeavyIons_withPF)
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

process.reco_extra_jet    = cms.Path( process.iterativeConePu5CaloJets * process.akPu3PFJets * process.akPu5PFJets * process.photon_extra_reco)
process.gen_step          = cms.Path( process.hiGenParticles * process.hiGenParticlesForJets * process.genPartons * process.hiPartons * process.hiRecoGenJets)
process.pat_step          = cms.Path( process.icPu5patSequence + process.akPu3PFpatSequence + process.akPu5PFpatSequence + process.makeHeavyIonPhotons)
process.pat_step.remove(process.interestingTrackEcalDetIds)
process.photonMatch.matched = cms.InputTag("hiGenParticles")
#process.pat_step.remove(process.photonMatch)
#+ process.patPhotons)

process.patPhotons.addPhotonID = cms.bool(False)
#process.makeHeavyIonPhotons)
process.extrapatstep = cms.Path(process.selectedPatPhotons)

process.multiPhotonAnalyzer.GammaEtaMax = cms.untracked.double(100)
process.multiPhotonAnalyzer.GammaPtMin = cms.untracked.double(0)
process.ana_step          = cms.Path( process.icPu5JetAnalyzer + process.akPu3PFJetAnalyzer +
                                      process.multiPhotonAnalyzer + process.anaTrack + process.pfcandAnalyzer +
                                      process.met * process.anaMET +
				      process.muonTree +
				      process.hiEvtAnalyzer +
                                      process.randomCones
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
process.rechitAna = cms.Path(process.rechitanalyzer)

########### random number seed

#####################################################################################
# Edm Output
#####################################################################################

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("output.root")
#                               )
#process.save = cms.EndPath(process.out)
