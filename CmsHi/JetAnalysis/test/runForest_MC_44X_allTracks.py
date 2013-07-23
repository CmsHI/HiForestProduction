import FWCore.ParameterSet.Config as cms

process = cms.Process('hiForestAna')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
    'file:MCRaw/step2_RAW2DIGI_RECO_1.root',
    ))

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(100))


#####################################################################################
# Load some general stuff
#####################################################################################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
# Data Global Tag 44x 
#process.GlobalTag.globaltag = 'GR_R_44_V6C::All'

# MC Global Tag 44x 
process.GlobalTag.globaltag = 'STARTHI44_V4::All'

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFhits"),
	nonDefaultGlauberModel = cms.string(""),
	centralitySrc = cms.InputTag("hiCentrality")
	)

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

# EcalSeverityLevel ES Producer
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load("RecoEcal.EgammaCoreTools.EcalNextToDeadChannelESProducer_cff")


#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("JetAnaTrees.root"))

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

# Remove neutrinos
process.hiGenParticlesForJets.ignoreParticleIDs += cms.vuint32( 12,14,16)

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

#process.hiGoodTightTracks.src = cms.InputTag("hiGlobalPrimTracks")
#process.hiGoodTightTracksDirect = process.hiGoodTightTracks.clone(keepAllTracks = True)
#process.hiGoodTracks = process.hiGoodTightTracks.clone()

process.reco_extra        = cms.Path( process.iterativeConePu5CaloJets *
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
				      process.hiEvtAnalyzer 
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
process.load('HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi')

process.hltbitanalysis.UseTFileService			= cms.untracked.bool(True)
process.hltanalysis = process.hltbitanalysis.clone(
   	            l1GtReadoutRecord            = cms.InputTag("gtDigis"),
	 	    l1GctHFBitCounts     = cms.InputTag("gctDigis"),
	 	    l1GctHFRingSums      = cms.InputTag("gctDigis"),
	 	    l1extramu            = cms.string('l1extraParticles'),
	 	    l1extramc            = cms.string('l1extraParticles'),
	 	    hltresults           = cms.InputTag("TriggerResults","","HLT"),
	 	    HLTProcessName       = cms.string("HLT")
	 	   )
process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","RECO")
process.skimanalysis = process.hltanalysis.clone(
    HLTProcessName                  = cms.string("JetAna"),
    hltresults                      = cms.InputTag("TriggerResults::hiForestAna")
    )
process.hlt = cms.Path(process.hltanalysis)
process.pAna = cms.EndPath(process.skimanalysis)
process.reco_extra*=process.L1Extra


process.load('CmsHi.JetAnalysis.rechitanalyzer_cfi')
  
process.rechitanalyzer.HBHETreePtMin = cms.untracked.double(-1000)
process.rechitanalyzer.HFTreePtMin = cms.untracked.double(-1000)
process.rechitanalyzer.EBTreePtMin = cms.untracked.double(-1000)
process.rechitanalyzer.EETreePtMin = cms.untracked.double(-1000)
process.rechitanalyzer.TowerTreePtMin = cms.untracked.double(-1000)
process.rechitanalyzer.doHF = cms.untracked.bool(True)
process.rechitAna = cms.Path(process.rechitanalyzer)

########### random number seed
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  multiPhotonAnalyzer = cms.PSet(
    engineName = cms.untracked.string("TRandom3"),
    initialSeed = cms.untracked.uint32(98236)
    )
  )

#####################################################################################
# Edm Output
#####################################################################################

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("output.root")
#                               )
#process.save = cms.EndPath(process.out)
