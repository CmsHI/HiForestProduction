
import FWCore.ParameterSet.Config as cms


process = cms.Process('HIJETS')

#process.Timing = cms.Service("Timing")

#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#                                        ignoreTotal=cms.untracked.int32(0),
#                                        oncePerEventMode = cms.untracked.bool(False)
#                                        )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_9_9/RelValQCD_Pt_80_120/GEN-SIM-RECO/START39_V9-v1/0000/0622F65F-623D-E011-A38D-002618943949.root'
    #'/store/relval/CMSSW_3_9_9/RelValSingleMuPt100/GEN-SIM-RECO/MC_39Y_V8-v1/0000/EEE53664-583D-E011-8433-0018F3D09692.root'
    #'/store/user/davidlw/Pyquen_UnquenchedDiJet_Pt15_START39V7HI_GEN_SIM_RAW_RECO_393_v1/Pyquen_UnquenchedDiJet_Pt15_START39V7HI_GEN_SIM_RAW_RECO_393_v1/24f3ed1e4d50df3073080480887969a1/Pyquen_UnquenchedDiJet_Pt15_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_100_1_dYJ.root'
)
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
                    )


useHighPtTrackCollection = False
useGoodTightCollection = True

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START39_V7::All'

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("CommonTools.UtilAlgos.TFileService_cfi")


process.HeavyIonGlobalParameters = cms.PSet(
           centralityVariable = cms.string("HFhits"),
           nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
           centralitySrc = cms.InputTag("hiCentrality")
                     )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)



#JEC from DB

process.load("CondCore.DBCommon.CondDBCommon_cfi")



if useGoodTightCollection:

    
    process.jec = cms.ESSource("PoolDBESSource",
                               DBParameters = cms.PSet(
                    messageLevel = cms.untracked.int32(0)
                    ),
                               timetype = cms.string('runnumber'),
                               toGet = cms.VPSet(
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_hiGoodTightTracks_D6T_399_IC5Calo"),
                 label = cms.untracked.string("IC5Calo")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_399_AK3PF"),
                 label = cms.untracked.string("AK3PF")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_399_AK4PF"),
                 label = cms.untracked.string("AK4PF")
                 ),
        
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
         tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_399_AK5PF"),
                 label = cms.untracked.string("AK5PF")
                 ),
        ),
                               connect = cms.string("sqlite_file:JEC_HI_PFTowers_hiGoodTightTracks_D6T_399.db"),
                               
)
    


elif useHighPtTrackCollection:

    process.jec = cms.ESSource("PoolDBESSource",
                               DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
                               timetype = cms.string('runnumber'),
                               toGet = cms.VPSet(
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_hiHighPtTracks_IC5Calo"),
                 label = cms.untracked.string("IC5Calo")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiHighPtTracks_AK3PF"),
                 label = cms.untracked.string("AK3PF")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiHighPtTracks_AK4PF"),
                 label = cms.untracked.string("AK4PF")
                 ),
        
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiHighPtTracks_AK5PF"),
                 label = cms.untracked.string("AK5PF")
                 ),
        ),
                               connect = cms.string("sqlite_file:JEC_HI_PFTowers_hiHighPtTracks_2011.db"),
                               
                               )
    

else:
    
      
    
    process.jec = cms.ESSource("PoolDBESSource",
                               DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
                               timetype = cms.string('runnumber'),
                               toGet = cms.VPSet(
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiSelFix_AK3PF"),
                 label = cms.untracked.string("AK3PF")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                   tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiSelFix_AK4PF"),
                 label = cms.untracked.string("AK4PF")
                 ),
        
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_hiSelFix_AK5PF"),
                 label = cms.untracked.string("AK5PF")
                 ),
        ),
                               connect = cms.string("sqlite_file:JEC_HI_2011.db"),
                               
                               )
    
      
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


# --- Pat Stuff --

# is this just setting defaults that I don't use anyway?  What is DisableMonteCarloDeps doing?  
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *

# This confusing script don't work no mo
#configureHeavyIons(process)

#let's do it by hand

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

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Summer09_7TeV_ReReco332")


#  --- Track and muon reconstruction ---
# pixel triplet tracking (HI Tracking)
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")
process.load("RecoHI.HiTracking.LowPtTracking_PbPb_cff")
# Needed to produce "HcalSeverityLevelComputerRcd" used by CaloTowersCreator/towerMakerPF
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")


# keep all the tracks for muon reco, then hack flags such that 'tight' tracks are input into final track selection and final tracks are set to 'highPurity'  (Need a cleaner way to do this)

if useGoodTightCollection:
    print "Using hiGlobalTight track collection"
    process.hiPixelAdaptiveVertex.useBeamConstraint = cms.bool(False)

elif useHighPtTrackCollection:
    print "Using hiHighPtTracks collection"

else:
    print "Using hiGoodTracks collection"
    process.hiTracksWithLooseQualityKeepAll = process.hiTracksWithLooseQuality.clone()
    process.hiTracksWithTightQualityKeepAll = process.hiTracksWithTightQuality.clone()
    process.hiSelectedTracksKeepAll = process.hiSelectedTracks.clone()    
    process.hiTracksWithTightQualityKeepAll.src = cms.InputTag("hiTracksWithLooseQualityKeepAll")
    process.hiSelectedTracksKeepAll.src = cms.InputTag("hiTracksWithTightQualityKeepAll")    
    process.hiTracksWithTightQuality.qualityBit = 'loose'
    process.hiSelectedTracks.qualityBit = 'tight'    
    process.hiTracksWithTightQualityKeepAll.qualityBit = 'loose'
    process.hiSelectedTracksKeepAll.qualityBit = 'tight'    
    process.hiTracksWithLooseQualityKeepAll.keepAllTracks = True
    process.hiTracksWithTightQualityKeepAll.keepAllTracks = True
    process.hiSelectedTracksKeepAll.keepAllTracks = True
    process.heavyIonTracking += process.hiTracksWithLooseQualityKeepAll*process.hiTracksWithTightQualityKeepAll*process.hiSelectedTracksKeepAll


# Muon Reco
from RecoHI.HiMuonAlgos.HiRecoMuon_cff import * 
muons.JetExtractorPSet.JetCollectionLabel = cms.InputTag("iterativeConePu5CaloJets")
#muons.JetExtractorPSet.JetCollectionLabel = cms.InputTag("iterativeCone5CaloJets")

if useGoodTightCollection:
    process.globalMuons.TrackerCollectionLabel = cms.InputTag("hiGoodTightTracksDirect")
    muons.TrackExtractorPSet.inputTrackCollection = cms.InputTag("hiGoodTightTracksDirect")
    muons.inputCollectionLabels = cms.VInputTag("hiGoodTightTracksDirect", "globalMuons", "standAloneMuons:UpdatedAtVtx")    
elif useHighPtTrackCollection:
    process.globalMuons.TrackerCollectionLabel = cms.InputTag("hiHighPtTracksKeepAll")
    muons.TrackExtractorPSet.inputTrackCollection = cms.InputTag("hiHighPtTracksKeepAll")
    muons.inputCollectionLabels = cms.VInputTag("hiHighPtTracksKeepAll", "globalMuons", "standAloneMuons:UpdatedAtVtx")
else:
    process.globalMuons.TrackerCollectionLabel = cms.InputTag("hiGoodTracksKeepAll")
    muons.TrackExtractorPSet.inputTrackCollection = cms.InputTag("hiGoodTracksKeepAll")
    muons.inputCollectionLabels = cms.VInputTag("hiGoodTracksKeepAll", "globalMuons", "standAloneMuons:UpdatedAtVtx")

process.muonRecoPbPb = muonRecoPbPb

#Track Reco
process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking)

# for PF
process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")

# good track selection
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")

#***********************************
#  for Efficiency analyzer need goodTight with Keep All off
# Dpt/pt = 0.05
process.hiGoodTightTracksDontKeepAll = process.hiGoodTightTracks.clone(keepAllTracks=False)
# Dpt/pt = 0.06
process.hiGoodTightTracksDirectDontKeepAll = process.hiGoodTightTracksDirect.clone(keepAllTracks=False)

process.hiGoodTightTracksSelection *= process.hiGoodTightTracksDontKeepAll
process.hiGoodTightTracksDirectSelection *= process.hiGoodTightTracksDirectDontKeepAll


process.hiGoodTracksKeepAll = process.hiGoodTracks.clone()
process.hiGoodTracksKeepAll.keepAllTracks = True
process.hiGoodTracksSelection += process.hiGoodTracksKeepAll
# "High pT track selection"
process.hiHighPtTracks.qualityBit = 'highPurity'
process.hiHighPtTracksKeepAll = process.hiHighPtTracks.clone()
process.hiHighPtTracksKeepAll.keepAllTracks = True
process.hiHighPtTrackSelection += process.hiHighPtTracksKeepAll
# merge with pixel tracks
process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
process.hiGoodMergedTracks.src = cms.InputTag("hiGoodTracks")



# particle-flow stuff

process.HiParticleFlowRecoNoJets = cms.Sequence(
    process.particleFlowCluster
    * process.trackerDrivenElectronSeeds
    * process.particleFlowReco
    )

if useGoodTightCollection:
    process.trackerDrivenElectronSeeds.TkColList = cms.VInputTag("hiGoodTightTracksDirect")
elif useHighPtTrackCollection:
    process.trackerDrivenElectronSeeds.TkColList = cms.VInputTag("hiHighPtTracksKeepAll")
else:
    process.trackerDrivenElectronSeeds.TkColList = cms.VInputTag("hiGoodTracksKeepAll")


process.load("HeavyIonsAnalysis.Configuration.analysisProducers_cff")
process.hiExtra = cms.Sequence(
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

# needed for HERWIG
process.hiGenParticles.abortOnUnknownPDGCode = cms.untracked.bool(False)
process.hiGenParticlesForJets.abortOnUnknownPDGCode = cms.untracked.bool(False)

#ignore neutrinos
process.hiGenParticlesForJets.ignoreParticleIDs += cms.vuint32(12,14,16)

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

process.patJets.addGenPartonMatch   = True
process.patJets.addGenJetMatch      = True


process.patJets.embedCaloTowers = cms.bool(False)


process.iterativeConePu5CaloJets.jetPtMin = 1.0
process.ak5PFJets.jetPtMin = 1.0


process.load("RecoHI.HiJetAlgos.ParticleTowerProducer_cff")
process.PFTowers = process.particleTowerProducer.clone()



process.iterativeConePu5CaloJets.doPVCorrection = cms.bool(True)
process.iterativeConePu5CaloJets.srcPVs = 'hiSelectedVertex'
process.iterativeConePu5CaloJets.doPUOffsetCorr = False


#Vanilla algo, but PU turned off

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


### AK5, PF  ###

process.ak5PFJets.src = 'PFTowers'
process.ak5PFJets.jetType = 'BasicJet'

process.ak5PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak5PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
                                                    )
process.ak5PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5HiGenJets'))
process.ak5PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ak5PFJets"),
                                                     matched = cms.InputTag("ak5PFclean"))
process.ak5PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ak5PFJets"))
process.ak5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak5PFJets"),
                                             genJetMatch = cms.InputTag("ak5PFmatch"),
                                             genPartonMatch= cms.InputTag("ak5PFparton"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFcorr")))
process.ak5PFpatSequence = cms.Sequence(process.ak5PFJets*process.ak5PFcorr*process.ak5PFclean*process.ak5PFmatch*process.ak5PFparton*process.ak5PFpatJets)




#begin small cones

process.ak3PFJets = process.ak5PFJets.clone()
process.ak3PFJets.rParam       = cms.double(0.3)

process.ak3PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak3PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK3PF')
                                                    )
process.ak3PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiGenJets'))
process.ak3PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ak3PFJets"),
                                                   matched = cms.InputTag("ak3PFclean"))
process.ak3PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ak3PFJets"))
process.ak3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak3PFJets"),
                                           genJetMatch = cms.InputTag("ak3PFmatch"),
                                           genPartonMatch= cms.InputTag("ak3PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3PFcorr")))
process.ak3PFpatSequence = cms.Sequence(process.ak3PFJets*process.ak3PFcorr*process.ak3PFclean*process.ak3PFmatch*process.ak3PFparton*process.ak3PFpatJets)

process.ak4PFJets = process.ak5PFJets.clone()
process.ak4PFJets.rParam       = cms.double(0.4)

process.ak4PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak4PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK4PF')
                                                    )
process.ak4PFclean = process.heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiGenJets'))
process.ak4PFmatch = process.patJetGenJetMatch.clone(src = cms.InputTag("ak4PFJets"),
                                                   matched = cms.InputTag("ak4PFclean"))
process.ak4PFparton = process.patJetPartonMatch.clone(src = cms.InputTag("ak4PFJets"))
process.ak4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak4PFJets"),
                                           genJetMatch = cms.InputTag("ak4PFmatch"),
                                           genPartonMatch= cms.InputTag("ak4PFparton"),                                       
                                           jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4PFcorr")))
process.ak4PFpatSequence = cms.Sequence(process.ak4PFJets*process.ak4PFcorr*process.ak4PFclean*process.ak4PFmatch*process.ak4PFparton*process.ak4PFpatJets)





process.runAllJets = cms.Sequence(
    process.icPu5patSequence +
    process.ak5PFpatSequence +
    process.ak3PFpatSequence +
    process.ak4PFpatSequence 
)

process.load("MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_ppMC_cff")
process.inclusiveJetAnalyzer.useCentrality = cms.untracked.bool(False)

process.inclusiveJetAnalyzerSequence = cms.Sequence(
    process.inclusiveJetAnalyzer
    *process.ak3PFJetAnalyzer
    *process.ak4PFJetAnalyzer
    *process.ak5PFJetAnalyzer
    )

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_ppMC_cff")
if useGoodTightCollection:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodTightTracksDirect")
elif useHighPtTrackCollection:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiHighPtTracks")
else:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodMergedTracks")
process.PFJetAnalyzer.hasSimInfo = cms.untracked.bool(True)
process.PFJetAnalyzer.writeGenParticles = cms.untracked.bool(True)
process.PFJetAnalyzer.jetTag2 = cms.InputTag("ak5PFpatJets")
process.PFJetAnalyzer.recoJetTag2 = cms.InputTag("ak5PFJets")
process.PFJetAnalyzer.jetTag3 = cms.InputTag("ak4PFpatJets")
process.PFJetAnalyzer.recoJetTag3 = cms.InputTag("ak4PFJets")
process.PFJetAnalyzer.jetTag4 = cms.InputTag("ak3PFpatJets")
process.PFJetAnalyzer.recoJetTag4 = cms.InputTag("ak3PFJets")


#Frank's analyzer
#import to avoid conflicts with inclusiveJetAnalyzer
from Saved.QM11Ana.Analyzers_cff import trkAnalyzer
from Saved.QM11Ana.Analyzers_cff import genpAnalyzer

from Saved.QM11Ana.Analyzers_cff import cutsTPForFak
from Saved.QM11Ana.Analyzers_cff import cutsTPForEff
from Saved.QM11Ana.Analyzers_cff import hitrkEffAnalyzer_akpu3pf_j1
from Saved.QM11Ana.Analyzers_cff import hitrkEffAnalyzer_akpu3pf_j2
from Saved.QM11Ana.Analyzers_cff import hitrkEffAnalyzer_akpu3pf

from Saved.QM11Ana.Analyzers_cff import cutsTPForFakPxl
from Saved.QM11Ana.Analyzers_cff import cutsTPForEffPxl
from Saved.QM11Ana.Analyzers_cff import hipixtrkEffAnalyzer_akpu3pf_j1
from Saved.QM11Ana.Analyzers_cff import hipixtrkEffAnalyzer_akpu3pf_j2
from Saved.QM11Ana.Analyzers_cff import hipixtrkEffAnalyzer_akpu3pf

hitrkEffAnalyzer_akpu3pf_j1.useQaulityStr = False
hitrkEffAnalyzer_akpu3pf_j2.useQaulityStr = False
hitrkEffAnalyzer_akpu3pf.useQaulityStr = False

if useGoodTightCollection:
    trkAnalyzer.trackSrc = 'hiGoodTightTracksDirectDontKeepAll'
    hitrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTightTracksDirectDontKeepAll'
    hitrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTightTracksDirectDontKeepAll'
    hitrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTightTracksDirectDontKeepAll'
    hipixtrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTightTracksDirectDontKeepAll'
    hipixtrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTightTracksDirectDontKeepAll'
    hipixtrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTightTracksDirectDontKeepAll'
elif useHighPtTrackCollection:
    trkAnalyzer.trackSrc = 'hiHighPtTracks'
    hitrkEffAnalyzer_akpu3pf_j1.tracks = 'hiHighPtTracks'
    hitrkEffAnalyzer_akpu3pf_j2.tracks = 'hiHighPtTracks'
    hitrkEffAnalyzer_akpu3pf.tracks = 'hiHighPtTracks'
    hipixtrkEffAnalyzer_akpu3pf_j1.tracks = 'hiHighPtTracks'
    hipixtrkEffAnalyzer_akpu3pf_j2.tracks = 'hiHighPtTracks'
    hipixtrkEffAnalyzer_akpu3pf.tracks = 'hiHighPtTracks'
else:
    trkAnalyzer.trackSrc = 'hiGoodMergedTracks'
    hitrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTracks'
    hitrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTracks'
    hitrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTracks'
    hipixtrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTracks'
    hipixtrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTracks'
    hipixtrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTracks'



process.hiGoodTightTrkEffAnalyzer_j1 = hitrkEffAnalyzer_akpu3pf_j1.clone(tracks = 'hiGoodTightTracksDontKeepAll')
process.hiGoodTightTrkEffAnalyzer_j2 = hitrkEffAnalyzer_akpu3pf_j2.clone(tracks = 'hiGoodTightTracksDontKeepAll')
process.hiGoodTightTrkEffAnalyzer = hitrkEffAnalyzer_akpu3pf.clone(tracks = 'hiGoodTightTracksDontKeepAll')
process.hiHighPtTrkEffAnalyzer_j1 = hitrkEffAnalyzer_akpu3pf_j1.clone(tracks = 'hiHighPtTracksDirect')
process.hiHighPtTrkEffAnalyzer_j2 = hitrkEffAnalyzer_akpu3pf_j2.clone(tracks = 'hiHighPtTracksDirect')
process.hiHighPtTrkEffAnalyzer = hitrkEffAnalyzer_akpu3pf.clone(tracks = 'hiHighPtTracksDirect')


process.hiGoodTightTrkEffAna = cms.Sequence(
    process.hiGoodTightTrkEffAnalyzer_j1*
    process.hiGoodTightTrkEffAnalyzer_j2*
    process.hiGoodTightTrkEffAnalyzer
    )

process.hiHighPtTrkEffAna = cms.Sequence(
    process.hiHighPtTrkEffAnalyzer_j1*
    process.hiHighPtTrkEffAnalyzer_j2*
    process.hiHighPtTrkEffAnalyzer
    )


process.trkAnalyzer = trkAnalyzer
process.cutsTPForFak = cutsTPForFak
process.cutsTPForEff = cutsTPForEff
process.hitrkEffAnalyzer_akpu3pf_j1 = hitrkEffAnalyzer_akpu3pf_j1
process.hitrkEffAnalyzer_akpu3pf_j2 = hitrkEffAnalyzer_akpu3pf_j2
process.hitrkEffAnalyzer_akpu3pf = hitrkEffAnalyzer_akpu3pf 
process.cutsTPForFakPxl = cutsTPForFakPxl
process.cutsTPForEffPxl = cutsTPForEffPxl
process.hipixtrkEffAnalyzer_akpu3pf_j1 = hipixtrkEffAnalyzer_akpu3pf_j1
process.hipixtrkEffAnalyzer_akpu3pf_j2 = hipixtrkEffAnalyzer_akpu3pf_j2
process.hipixtrkEffAnalyzer_akpu3pf = hipixtrkEffAnalyzer_akpu3pf 
process.genpAnalyzer = genpAnalyzer


process.trkAnalyzer.trackPtMin = 0.5
process.genpAnalyzer.ptMin = 0.5

process.hitrkEffAnalyzer_akpu3pf.jets = 'ak3PFpatJets'
process.hitrkEffAnalyzer_akpu3pf_j1.jets = 'ak3PFpatJets'
process.hitrkEffAnalyzer_akpu3pf_j2.jets = 'ak3PFpatJets'

process.hipixtrkEffAnalyzer_akpu3pf.jets = 'ak3PFpatJets'
process.hipixtrkEffAnalyzer_akpu3pf_j1.jets = 'ak3PFpatJets'
process.hipixtrkEffAnalyzer_akpu3pf_j2.jets = 'ak3PFpatJets'

process.hitrkEffAna_akpu3pf = cms.Sequence(process.cutsTPForFak*process.cutsTPForEff*process.hitrkEffAnalyzer_akpu3pf*process.hitrkEffAnalyzer_akpu3pf_j1*process.hitrkEffAnalyzer_akpu3pf_j2)
process.hipixtrkEffAna_akpu3pf = cms.Sequence(process.cutsTPForFakPxl*process.cutsTPForEffPxl*process.hipixtrkEffAnalyzer_akpu3pf*process.hipixtrkEffAnalyzer_akpu3pf_j1*process.hipixtrkEffAnalyzer_akpu3pf_j2)


#process.franksAnalyzers = cms.Sequence(process.trkAnalyzer*process.hitrkEffAna_akpu3pf*process.hipixtrkEffAna_akpu3pf*process.genpAnalyzer)
process.franksAnalyzers = cms.Sequence(process.trkAnalyzer*process.hitrkEffAna_akpu3pf*process.hiGoodTightTrkEffAna*process.hiHighPtTrkEffAna*process.genpAnalyzer)

process.load("PbPbTrackingTools.HiTrackValidator.hitrackvalidator_cfi")
if useHighPtTrackCollection:
    process.hitrkvalidator.trklabel=cms.untracked.InputTag("hiHighPtTracks")
else:
    process.hitrkvalidator.trklabel=cms.untracked.InputTag("hiGoodTracks")


process.higoodTrkVal = process.hitrkvalidator.clone(
    trklabel = "hiGoodTracks",
    neededCentBins = [0,1,3,11,19,35],
    jetlabel = "akPu3PFpatJets")
process.higoodtightTrkVal = process.higoodTrkVal.clone(
    trklabel=cms.untracked.InputTag("hiGoodTightTracksDirect")
    )
process.hihighTrkVal = process.higoodTrkVal.clone(trklabel=cms.untracked.InputTag("hiHighPtTracksDirect"),
    qualityString = "tight")

process.hiCentrality.produceHFhits = cms.bool(True)
process.hiCentrality.produceHFtowers = cms.bool(False)
process.hiCentrality.produceEcalhits = cms.bool(False)
process.hiCentrality.produceBasicClusters = cms.bool(False)
process.hiCentrality.produceZDChits = cms.bool(False)
process.hiCentrality.produceETmidRapidity = cms.bool(False)
process.hiCentrality.producePixelhits = cms.bool(False)
process.hiCentrality.produceTracks = cms.bool(False)
process.hiCentrality.producePixelTracks = cms.bool(False)
process.hiCentrality.doPixelCut = cms.bool(True)

process.hiValidatorSequence = cms.Sequence(
    process.hiCentrality*
    process.higoodTrkVal*
    process.higoodtightTrkVal*
    process.hihighTrkVal
)


#Yetkin's analyzer
from CmsHi.JetAnalysis.ak3PFTowerJetsAna_cff import *
process.yetkinsAna = akPu3PFtowerJetsAna.clone()
process.RandomNumberGeneratorService.yetkinsAna = process.RandomNumberGeneratorService.generator.clone(initialSeed = 1)
process.yetkinsAna.patJetSrc = cms.untracked.InputTag("ak3PFpatJets")
process.yetkinsAna.doRecoEvtPlane = cms.untracked.bool(False)

#for tree output
process.TFileService = cms.Service("TFileService",
                                   #fileName=cms.string("JetAnalysisTTrees_hiGoodMergedTracks_seedGoodTracks_v1.root")
                                   #fileName=cms.string("JetAnalysisTTrees_hiHighPtTracks_v2.root")
                                   fileName=cms.string("JetAnalysisTTrees_hiGoodTightTracks_v1.root")
                                   )


# put it all together
if useGoodTightCollection:
    process.trackRecoAndSelection = cms.Path(
        process.hiTrackReco*
         process.hiPostGlobalPrimTracks*
        process.hiGoodTracksSelection*
        process.hiGoodTightTracksSelection*
        process.hiGoodTightTracksDirectSelection*
        process.hiHighPtTrackDirectSelection*
        process.muonRecoPbPb 
        )
elif useHighPtTrackCollection:
    process.trackRecoAndSelection = cms.Path(
        process.hiTrackReco*
        process.hiHighPtTrackSelection*
        process.muonRecoPbPb 
        )
else:
    process.trackRecoAndSelection = cms.Path(
        process.hiTrackReco*
        process.hiGoodTracksSelection*
        process.conformalPixelTrackReco *
        process.hiGoodMergedTracks *
        process.muonRecoPbPb 
        )


process.jetReco = cms.Path(
    process.HiParticleFlowRecoNoJets *
    process.hiExtra*
    process.hiGen*
    process.PFTowers*
    process.runAllJets
    )
process.jetAna = cms.Path(
    #process.franksAnalyzers*
    process.inclusiveJetAnalyzerSequence
    #*process.PFJetAnalyzerSequence
    *process.hiValidatorSequence
    *process.yetkinsAna
    )

process.schedule = cms.Schedule(process.trackRecoAndSelection, process.jetReco, process.jetAna)

