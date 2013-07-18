import FWCore.ParameterSet.Config as cms

process = cms.Process('hiForestAna2012')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#####################################################################################
# High Level Setup
#####################################################################################
import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.register ('randomNumber',
                1,
                ivars.multiplicity.singleton,
                ivars.varType.int,
                "Random Seed")

ivars.register ('doSkim',
               False, # default value
               VarParsing.VarParsing.multiplicity.singleton, # singleton or list
               VarParsing.VarParsing.varType.bool,          # string, int, or float
               "whether to do skimming")
                  
ivars.randomNumber = 1
ivars.inputFiles = "file:/mnt/hadoop/cms/store/user/frankmalocal/data/reco/96A31CE6-1921-E111-B606-00A0D1E953AA.root"
ivars.outputFile = 'output_data.root'
ivars.maxEvents = -1

ivars.parseArguments()

process.load("CmsHi.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V2",
                                          "PF : iterative tracks",
                                          "EP Flattening OFF",
                                          "Electrons OFF",
                                          "Preshower OFF"
                                          )

#####################################################################################
# Input source
#####################################################################################
process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(ivars.inputFiles))

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(ivars.maxEvents))

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
# process.load('HLTrigger.Configuration.HLT_HIon_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
# Data Global Tag 44x 
process.GlobalTag.globaltag = 'GR_P_V27A::All'

# process.Timing = cms.Service("Timing")

# MC Global Tag 44x 
# process.GlobalTag.globaltag = 'STARTHI44_V7::All'

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
	centralitySrc = cms.InputTag("hiCentrality")
	)
process.hiCentrality.pixelBarrelOnly = False

# EcalSeverityLevel ES Producer
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load("RecoEcal.EgammaCoreTools.EcalNextToDeadChannelESProducer_cff")

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string(ivars.outputFile))

#####################################################################################
# Additional Reconstruction and Analysis
#####################################################################################
process.load("CmsHi.JetAnalysis.hcalNoise_cff")
process.load('CmsHi.JetAnalysis.EventSelection_cff')
process.load('CmsHi.JetAnalysis.ExtraGenReco_cff')
process.load('CmsHi.JetAnalysis.ExtraTrackReco_cff')
process.load('CmsHi.JetAnalysis.ExtraPfReco_cff')
process.load('CmsHi.JetAnalysis.ExtraJetReco_cff')
process.load('CmsHi.JetAnalysis.ExtraEGammaReco_cff')
process.load('CmsHi.JetAnalysis.PatAna_cff')
process.load('CmsHi.JetAnalysis.JetAnalyzers_cff')
process.load('CmsHi.JetAnalysis.TrkAnalyzers_cff')
process.load('CmsHi.JetAnalysis.EGammaAnalyzers_cff')

##################### Track Related
process.load("RecoHI.HiTracking.hiIterTracking_cff")
process.heavyIonTracking *= process.hiIterTracking

hiTrackQuality = "highPurity"              # iterative tracks
#hiTrackQuality = "highPuritySetWithPV"    # calo-matched tracks
process.hiTracks.src = "hiGeneralTracks"
process.hiTracks.cut = cms.string('quality("' + hiTrackQuality+  '")')
process.heavyIonTracking *= process.hiTracks

##################### Jet Related
# Remove neutrinos
process.hiGenParticlesForJets.ignoreParticleIDs += cms.vuint32( 12,14,16)

##################### Photon Related
# EcalSeverityLevel ES Producer
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load("RecoEcal.EgammaCoreTools.EcalNextToDeadChannelESProducer_cff")
process.hiSelectedTrackHighPurity = process.hiSelectedTrackQuality.clone()

##################### Particle Flow
process.particleFlowTmp.postMuonCleaning = False
process.particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)
process.pfTrack.TrackQuality = cms.string('loose') # To recover tracks before calo matching
process.pfTrack.TkColList = cms.VInputTag("hiGeneralTracks") # needed to run the calo matching on hiGeneralTracks

##################### Pat
process.interestingTrackEcalDetIds.TrackCollection = cms.InputTag("hiGeneralTracks")
process.patPhotons.addPhotonID = cms.bool(False)
process.photonMatch.matched = cms.InputTag("hiGenParticles")

################################################################
# Analyzers
################################################################
##################### Event Tree Analyzers
process.load("CmsHi/HiHLTAlgos.hievtanalyzer_cfi")
process.hiEvtAnalyzer.doEvtPlane = cms.bool(True)

##################### Jet Tree Analyzers
for m in [
   process.icPu5JetAnalyzer,
   process.akPu1PFJetAnalyzer,
   process.akPu2PFJetAnalyzer,
   process.akPu3PFJetAnalyzer,
   process.akPu4PFJetAnalyzer,
   process.akPu5PFJetAnalyzer,
   process.akPu6PFJetAnalyzer,
   process.akPu1CaloJetAnalyzer,
   process.akPu2CaloJetAnalyzer,
   process.akPu3CaloJetAnalyzer,
   process.akPu4CaloJetAnalyzer,
   process.akPu5CaloJetAnalyzer,
   process.akPu6CaloJetAnalyzer,
   process.ak1PFJetAnalyzer,
   process.ak2PFJetAnalyzer,
   process.ak3PFJetAnalyzer,
   process.ak4PFJetAnalyzer,
   process.ak5PFJetAnalyzer,
   process.ak6PFJetAnalyzer,
   process.ak1CaloJetAnalyzer,
   process.ak2CaloJetAnalyzer,
   process.ak3CaloJetAnalyzer,
   process.ak4CaloJetAnalyzer,
   process.ak5CaloJetAnalyzer,
   process.ak6CaloJetAnalyzer
   ]:
   m.trackTag = "hiTracks"
   m.hltTrgResults = cms.untracked.string('TriggerResults::HLT')

##################### Photon Tree Analyzers
process.load('CmsHi.JetAnalysis.EGammaAnalyzers_cff')
process.multiPhotonAnalyzer.GammaEtaMax = cms.untracked.double(100)
process.multiPhotonAnalyzer.GammaPtMin = cms.untracked.double(0)
process.multiPhotonAnalyzer.gsfElectronCollection = cms.untracked.InputTag("ecalDrivenGsfElectrons")


##################### Track Analyzers
process.anaTrack.trackSrc = cms.InputTag("hiTracks")
process.anaTrack.qualityStrings = cms.untracked.vstring('highPurity')

############################################ Pf Analyzers
process.load("CmsHi.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0.5

############################################ Gen Particle Analyzers
process.HiGenParticleAna = cms.EDAnalyzer('HiGenAnalyzer',
    useHepMCProduct = cms.untracked.bool(False),
    ptMin = cms.untracked.double(0.8),
    chargedOnly = cms.untracked.bool(True),
    src = cms.untracked.InputTag("hiSignal"),
    genpSrc = cms.untracked.InputTag("hiGenParticles"), # save just the signal genp
    genHiSrc = cms.untracked.InputTag("heavyIon"),
    )
# if saveAllGenParticles:
#     process.HiGenParticleAna.genpSrc = cms.untracked.InputTag("hiGenParticles","","HISIGNAL")

process.genpana = cms.EDAnalyzer("GenParticleCounter",
                              src = cms.untracked.string("hiGenParticles"),
                              doCentrality = cms.untracked.bool(False),
                              VertexProducer = cms.untracked.string("hiSelectedVertex")
                              )

############################################ Other MC Info

##################### Final Paths
process.reco_extra =  cms.Path(
   process.hiTrackReReco *
   (process.electronGsfTrackingHi * process.hiElectronSequence) *
   process.HiParticleFlowLocalReco *
   process.HiParticleFlowReco *
   process.PFTowers
)    
    

process.reco_extra_jet    = cms.Path(process.iterativeConePu5CaloJets *
                                     process.iterativeCone5CaloJets *
                                     process.recoAk1to6 *
                                     process.photon_extra_reco)

process.pat_step          = cms.Path(process.makeHeavyIonJets +
                                     process.makeHeavyIonPhotons
                                     )    
process.pat_step.remove(process.interestingTrackEcalDetIds)

process.extrapatstep = cms.Path(process.selectedPatPhotons)

process.ana_step          = cms.Path(
                                      process.hcalNoise +
                                      process.jetAnalyzers +
                                      process.multiPhotonAnalyzer +
                                      (process.anaTrack) +
                                      process.pfcandAnalyzer +
                                      process.hiEvtAnalyzer +
                                      process.HiForest
                                      )


process.phltJetHI = cms.Path( process.hltJetHI )
process.pcollisionEventSelection = cms.Path(process.collisionEventSelection)
process.pHBHENoiseFilter = cms.Path( process.HBHENoiseFilter )
process.phiEcalRecHitSpikeFilter = cms.Path(process.hiEcalRecHitSpikeFilter )

# Customization
from CmsHi.JetAnalysis.customise_cfi import *
setPhotonObject(process,"cleanPhotons")

#############################################################################
# Event Triggering/Filtering
#############################################################################
# process.load('L1Trigger.Configuration.L1Extra_cff')
process.load('CmsHi.HiHLTAlgos.hltanalysis_cff')
process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HLT")
process.hltAna = cms.Path(process.hltanalysis)
# process.reco_extra*=process.L1Extra

process.skimanalysis.hltresults = cms.InputTag("TriggerResults","",process.name_())
process.pAna = cms.EndPath(process.skimanalysis)
process.endjob_step = cms.EndPath(process.endOfProcess)

#############################################################################
# Final Schedule
#############################################################################
# process.schedule = cms.Schedule(process.L1simulation_step,process.reco_extra, process.reco_extra_jet, process.gen_step, process.pat_step, process.extrapatstep,process.ana_step, process.phltJetHI,process.pcollisionEventSelection,process.pHBHENoiseFilter,process.phiEcalRecHitSpikeFilter,process.hltAna,process.pAna)
process.schedule = cms.Schedule(process.reco_extra, process.reco_extra_jet,process.pat_step, process.extrapatstep,process.ana_step, process.pcollisionEventSelection,process.pHBHENoiseFilter,process.phiEcalRecHitSpikeFilter,process.hltAna,process.pAna)

# process.HLTSchedule.extend(process.schedule)



########### random number seed
process.RandomNumberGeneratorService.generator.initialSeed = ivars.randomNumber 
process.RandomNumberGeneratorService.multiPhotonAnalyzer = process.RandomNumberGeneratorService.generator.clone()

#####################################################################################
# Skimming
#####################################################################################
if (ivars.doSkim):
   process.schedule.remove(process.pAna)
   import HLTrigger.HLTfilters.hltHighLevel_cfi
   process.photonFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
   process.photonFilter.HLTPaths = ["HLT_HISinglePhoton30_v*"]
   process.singlePhotonPtFilter = cms.EDFilter("PhotonSelector",
                                                src = cms.InputTag("photons"),
                                                cut = cms.string('pt > 45 && abs(eta) < 1.48 && sigmaIetaIeta > 0.002' ),
                                                filter = cms.bool(True)
                                                )
   
   process.superFilterSequence = cms.Sequence(process.photonFilter*process.collisionEventSelection*process.singlePhotonPtFilter*process.HBHENoiseFilter)
   process.superFilterPath = cms.Path(process.superFilterSequence)
   
   for path in process.paths:
       getattr(process,path)._seq = process.superFilterSequence*getattr(process,path)._seq

#####################################################################################
# Edm Output
#####################################################################################

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("output.root")
#                               )
#process.save = cms.EndPath(process.out)
