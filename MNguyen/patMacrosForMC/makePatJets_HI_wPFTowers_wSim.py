import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')

#process.Timing = cms.Service("Timing")

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'/store/relval/CMSSW_3_10_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/START310_V3-v1/0050/54ED7C66-680E-E011-BD9F-001A92810ADE.root'
    #'/store/relval/CMSSW_3_9_9/RelValHydjetQ_B0_2760GeV/GEN-SIM-RECO/START39_V7HI-v1/0001/0618A00D-E23D-E011-A7FE-001A9281173E.root'
    'file:/d101/frankma/ana/qm11jet/hiReco_RAW2DIGI_RECO_9_1_QNo.root'
#'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_9/RelValPyquen_DiJet_pt80to120_2760GeV/GEN-SIM-RECO/START39_V7HI-v1/0000/3AD28469-FE3D-E011-A504-003048D15E24.root'
    #'/store/user/davidlw/Hydjet_Bass_MinBias_2760GeV/Pyquen_UnquenchedDiJet_Pt50_GEN-SIM-RECO_393_v2/1243c1b8707a4e7eb28eae64e1474920/hiReco_RAW2DIGI_RECO_9_1_ZVb.root'
    )
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(2)
    )


isMinBias = False
useHighPtTrackCollection = False
useGoodTightCollection = True

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START39_V7::All'

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

# load centrality

process.load("CommonTools.UtilAlgos.TFileService_cfi")

process.HeavyIonGlobalParameters = cms.PSet(
           centralityVariable = cms.string("HFhits"),
           nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
           centralitySrc = cms.InputTag("hiCentrality")
                     )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

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
                 tag = cms.string("JetCorrectorParametersCollection_HI_fromFit_IC5Calo"),
                 label = cms.untracked.string("IC5Calo")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_EscaleXcheck_v2_AK3PF"),
                 label = cms.untracked.string("AK3PF")
                 ),
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_EscaleXcheck_AK4PF"),
                 label = cms.untracked.string("AK4PF")
                 ),
        
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                 tag = cms.string("JetCorrectorParametersCollection_HI_PFTowers_EscaleXcheck_AK5PF"),
                 label = cms.untracked.string("AK5PF")
                 ),
        ),
                               connect = cms.string("sqlite_file:JEC_HI_PFTowers_hiGoodTracks_v2_2011.db"),
                               
                               )

process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

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

# keep all the tracks for muon reco, then hack flags such that 'tight' tracks are input into final track selection and final tracks are set to 'highPurity'  (Need a cleaner way to do this)

if useGoodTightCollection:
    print "Using hiGoodTight track collection"
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

process.muonRecoPbpb = muonRecoPbPb

#Track Reco
process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking)



# good track selection
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")
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



# for PF
process.load("RecoHI.Configuration.Reconstruction_hiPF_cff")


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
#This runs GenHIEventProducer
process.heavyIon


# --- Gen stuff ---
# Note that we need to re-run gen jets b/c until AK-IC bug fix propagates to production
process.load('PhysicsTools.HepMCCandAlgos.HiGenParticles_cfi')

process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoHI.HiJetAlgos.HiGenJets_cff')



if isMinBias:
    process.hiGenParticles.srcVector = cms.vstring('generator')
else:
    process.hiGenParticles.srcVector = cms.vstring('hiSignal')


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

    
process.iterativeConePu5CaloJets.puPtMin = cms.double(10.0)
process.ak5PFJets.puPtMin = cms.double(25.0)
process.ak4PFJets.puPtMin = cms.double(20.0)
process.ak3PFJets.puPtMin = cms.double(15.0)


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
                                                    payload = cms.string('AK4PF')
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
                                                    payload = cms.string('AK3PF')
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
if useGoodTightCollection:
    process.trackExtrapolator.trackSrc = cms.InputTag("hiGoodTightTracksDirect")
    process.JPTiterativeConePu5JetTracksAssociatorAtVertex.tracks = 'hiGoodTightTracksDirect'
    process.JPTiterativeConePu5JetTracksAssociatorAtCaloFace.tracks = 'hiGoodTightTracksDirect'
    process.JetPlusTrackZSPCorJetIconePu5.tracks = 'hiGoodTightTracksDirect'
elif useHighPtTrackCollection:
    process.trackExtrapolator.trackSrc = cms.InputTag("hiHighPtTracks")
    process.JPTiterativeConePu5JetTracksAssociatorAtVertex.tracks = 'hiHighPtTracks'
    process.JPTiterativeConePu5JetTracksAssociatorAtCaloFace.tracks = 'hiHighPtTracks'
    process.JetPlusTrackZSPCorJetIconePu5.tracks = 'hiHighPtTracks'
else:
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
if isMinBias:
    process.inclusiveJetAnalyzer.eventInfoTag = cms.InputTag("generator")
else:
    process.inclusiveJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
process.icPu5JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu5PFJetAnalyzer.jetTag = 'akPu5PFpatJets'
process.akPu5PFJetAnalyzer.genjetTag = 'ak5HiGenJets'
process.akPu4PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu4PFJetAnalyzer.jetTag = 'akPu4PFpatJets'
process.akPu4PFJetAnalyzer.genjetTag = 'ak4HiGenJets'
process.akPu3PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu3PFJetAnalyzer.jetTag = 'akPu3PFpatJets'
process.akPu3PFJetAnalyzer.genjetTag = 'ak3HiGenJets'
process.icPu5JPTJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.icPu5JPTJetAnalyzer.jetTag = 'jpticPu5patJets'


process.inclusiveJetAnalyzerSequence = cms.Sequence(
    #process.icPu5JetAnalyzer
    #*process.akPu5PFJetAnalyzer
    #*process.akPu4PFJetAnalyzer
    #*process.akPu3PFJetAnalyzer
    process.icPu5JPTJetAnalyzer
    )

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.load("edwenger/HiTrkEffAnalyzer/hitrackingParticleSelector_cfi")
process.hitrackingParticleSelector.ptMin = 10

process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_cff")
if useGoodTightCollection:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodTightTracksDirect")
elif useHighPtTrackCollection:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiHighPtTracks")
else:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodMergedTracks")
    
process.PFJetAnalyzer.hasSimInfo = cms.untracked.bool(True)
process.PFJetAnalyzer.SimTracks = cms.InputTag("hitrackingParticleSelector")
process.PFJetAnalyzer.doSimTracks = cms.untracked.bool(True)

if isMinBias:
    process.PFJetAnalyzer.eventInfoTag = cms.InputTag("generator")
    process.PFJetAnalyzer.genParticleTrhesh = cms.double(2.0)
else:
    process.PFJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
    process.PFJetAnalyzer.genParticleTrhesh = cms.double(0.5)
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

process.hiValidatorSequence = cms.Sequence(
    process.higoodTrkVal*
    process.higoodtightTrkVal*
    process.hihighTrkVal
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

if useGoodTightCollection:
  #trkAnalyzer.trackSrc = 'hiGoodTightTracksDirect'
    hitrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTightTracksDirect'
    hitrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTightTracksDirect'
    hitrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTightTracksDirect'
    hipixtrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTightTracksDirect'
    hipixtrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTightTracksDirect'
    hipixtrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTightTracksDirect'
elif useHighPtTrackCollection:
  #trkAnalyzer.trackSrc = 'hiHighPtTracks'
    hitrkEffAnalyzer.qualityString = cms.untracked.string("tight")
    hitrkEffAnalyzer_akpu3pf_j1.tracks = 'hiHighPtTracks'
    hitrkEffAnalyzer_akpu3pf_j2.tracks = 'hiHighPtTracks'
    hitrkEffAnalyzer_akpu3pf.tracks = 'hiHighPtTracks'
    hipixtrkEffAnalyzer_akpu3pf_j1.tracks = 'hiHighPtTracks'
    hipixtrkEffAnalyzer_akpu3pf_j2.tracks = 'hiHighPtTracks'
    hipixtrkEffAnalyzer_akpu3pf.tracks = 'hiHighPtTracks'
else:
  #trkAnalyzer.trackSrc = 'hiGoodMergedTracks'
    hitrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTracks'
    hitrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTracks'
    hitrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTracks'
    hipixtrkEffAnalyzer_akpu3pf_j1.tracks = 'hiGoodTracks'
    hipixtrkEffAnalyzer_akpu3pf_j2.tracks = 'hiGoodTracks'
    hipixtrkEffAnalyzer_akpu3pf.tracks = 'hiGoodTracks'

process.hiGoodTightTrkEffAnalyzer_j1 = hitrkEffAnalyzer_akpu3pf_j1.clone(tracks = 'hiGoodTightTracksDirect')
process.hiGoodTightTrkEffAnalyzer_j2 = hitrkEffAnalyzer_akpu3pf_j2.clone(tracks = 'hiGoodTightTracksDirect')
process.hiGoodTightTrkEffAnalyzer = hitrkEffAnalyzer_akpu3pf.clone(tracks = 'hiGoodTightTracksDirect')
process.hiGoodTightTrkEffAnalyzer_icpu5calo_j1 = process.hiGoodTightTrkEffAnalyzer_j1.clone(jets = 'icPu5patJets')
process.hiGoodTightTrkEffAnalyzer_icpu5calo_j2 = process.hiGoodTightTrkEffAnalyzer_j2.clone(jets = 'icPu5patJets')
process.hiGoodTightTrkEffAnalyzer_icpu5calo = process.hiGoodTightTrkEffAnalyzer.clone(jets = 'icPu5patJets')
process.hiHighPtTrkEffAnalyzer_j1 = hitrkEffAnalyzer_akpu3pf_j1.clone(
    qualityString = cms.untracked.string("tight"),
    tracks = 'hiHighPtTracksDirect')
process.hiHighPtTrkEffAnalyzer_j2 = hitrkEffAnalyzer_akpu3pf_j2.clone(
    qualityString = cms.untracked.string("tight"),
    tracks = 'hiHighPtTracksDirect')
process.hiHighPtTrkEffAnalyzer = hitrkEffAnalyzer_akpu3pf.clone(
    qualityString = cms.untracked.string("tight"),
    tracks = 'hiHighPtTracksDirect')


process.hiGoodTightTrkEffAna = cms.Sequence(
    process.hiGoodTightTrkEffAnalyzer_j1*
    process.hiGoodTightTrkEffAnalyzer_j2*
    process.hiGoodTightTrkEffAnalyzer *
    process.hiGoodTightTrkEffAnalyzer_icpu5calo_j1*
    process.hiGoodTightTrkEffAnalyzer_icpu5calo_j2*
    process.hiGoodTightTrkEffAnalyzer_icpu5calo
    )

process.hiHighPtTrkEffAna = cms.Sequence(
    process.hiHighPtTrkEffAnalyzer_j1*
    process.hiHighPtTrkEffAnalyzer_j2*
    process.hiHighPtTrkEffAnalyzer
    )


process.trkAnalyzer.trackPtMin = 0.5

if isMinBias:
    process.genpAnalyzer.ptMin = 0.5
else:
    process.genpAnalyzer.ptMin = 2.0



process.hitrkEffAna_akpu3pf = cms.Sequence(process.cutsTPForFak*process.cutsTPForEff*process.hitrkEffAnalyzer_akpu3pf*process.hitrkEffAnalyzer_akpu3pf_j1*process.hitrkEffAnalyzer_akpu3pf_j2)
process.hipixtrkEffAna_akpu3pf = cms.Sequence(process.cutsTPForFakPxl*process.cutsTPForEffPxl*process.hipixtrkEffAnalyzer_akpu3pf*process.hipixtrkEffAnalyzer_akpu3pf_j1*process.hipixtrkEffAnalyzer_akpu3pf_j2)


#process.franksAnalyzers = cms.Sequence(process.trkAnalyzer*process.hitrkEffAna_akpu3pf*process.hipixtrkEffAna_akpu3pf*process.genpAnalyzer)
process.franksAnalyzers = cms.Sequence(process.trkAnalyzer*process.hitrkEffAna_akpu3pf*process.hiGoodTightTrkEffAna*process.hiHighPtTrkEffAna*process.genpAnalyzer)

#Yetkin's analyzer
from CmsHi.JetAnalysis.ak3PFTowerJetsAna_cff import *
process.yetkinsAna = akPu3PFtowerJetsAna.clone()
process.RandomNumberGeneratorService.yetkinsAna = process.RandomNumberGeneratorService.generator.clone(initialSeed = 1)

#for tree output
process.TFileService = cms.Service("TFileService",
                                   #fileName=cms.string("JetAnalysisTTrees_hiGoodMergedTracks_seedGoodTracks_v1.root"))
                                   fileName=cms.string("JetAnalysisTTrees_hiGoodTightTracks_v1.root"))



# put it all together
if useGoodTightCollection:
    process.trackRecoAndSelection = cms.Path(
        process.hiTrackReco*
        process.hiGoodTracksSelection*
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
    process.franksAnalyzers*
    process.hiValidatorSequence*
    #process.inclusiveJetAnalyzerSequence*
    process.hitrackingParticleSelector *
    process.PFJetAnalyzerSequence
    *process.yetkinsAna
    )

process.schedule = cms.Schedule(process.trackRecoAndSelection, process.jetReco, process.jetAna)





#process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")
#
#process.output = cms.OutputModule("PoolOutputModule",
#                                  process.jetTrkSkimContentMC,
#                                  fileName = cms.untracked.string("/tmp/mnguyen/PAT.root")
#                                  )
#
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

