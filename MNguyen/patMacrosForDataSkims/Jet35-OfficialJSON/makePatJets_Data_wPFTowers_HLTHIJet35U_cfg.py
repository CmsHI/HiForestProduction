import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #"/store/hidata/HIRun2010/HIAllPhysics/RECO/SDmaker_3SD_1CS_PDHIAllPhysicsZSv2_SD_JetHI-v1/0055/FAA6A135-B04F-E011-BCB4-003048FEAED8.root",
    "/store/hidata/HIRun2010/HIAllPhysics/RECO/SDmaker_3SD_1CS_PDHIAllPhysicsZSv2_SD_JetHI-v1/0020/FC80C7C2-0B4C-E011-A1E8-003048FEC040.root"
    #"file:/tmp/mnguyen/FAA6A135-B04F-E011-BCB4-003048FEAED8.root",
    #"file:/tmp/mnguyen/FC80C7C2-0B4C-E011-A1E8-003048FEC040.root"
    )

                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

useHighPtTrackCollection = False
useGoodTightCollection = True

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR_R_39X_V6B::All'


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


# --- Pat Stuff --
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
process.patJets.addGenPartonMatch   = False
process.patJets.addGenJetMatch      = False
process.patJets.embedGenJetMatch    = False
process.patJets.embedGenPartonMatch   = False


#  --- Track and muon reconstruction ---
# pixel triplet tracking (HI Tracking)
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")
#process.load("RecoHI.HiTracking.LowPtTracking_PbPb_cff")
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

 #good track selection
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




process.icPu5corr = process.patJetCorrFactors.clone(src = cms.InputTag("iterativeConePu5CaloJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
process.icPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("iterativeConePu5CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5corr")))
process.icPu5patSequence = cms.Sequence(process.iterativeConePu5CaloJets*process.icPu5corr*process.icPu5patJets)




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
process.akPu5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("akPu5PFJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu5PFcorr")))
process.akPu5PFpatSequence = cms.Sequence(process.akPu5PFJets*process.akPu5PFcorr*process.akPu5PFpatJets)



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
process.akPu4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("akPu4PFJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu4PFcorr")))
process.akPu4PFpatSequence = cms.Sequence(process.akPu4PFJets*process.akPu4PFcorr*process.akPu4PFpatJets)



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
process.akPu3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("akPu3PFJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPu3PFcorr")))
process.akPu3PFpatSequence = cms.Sequence(process.akPu3PFJets*process.akPu3PFcorr*process.akPu3PFpatJets)







# JPT
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
process.jpticPu5patJets = process.patJets.clone(jetSource  =cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jpticPu5corr")))
process.icPu5JPTpatSequence = cms.Sequence(process.recoJPTJetsHIC*process.jpticPu5corr*process.jpticPu5patJets)

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
# trees need centrality
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string(""),
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

process.load("MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_cff")
process.inclusiveJetAnalyzer.isMC = cms.untracked.bool(False)
process.icPu5JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.icPu5JetAnalyzer.jetTag = 'icPu5patJets'
process.akPu5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu5PFJetAnalyzer.jetTag = 'akPu5PFpatJets'
process.akPu4PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu4PFJetAnalyzer.jetTag = 'akPu4PFpatJets'
process.akPu3PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.akPu3PFJetAnalyzer.jetTag = 'akPu3PFpatJets'
process.icPu5JPTJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.icPu5JPTJetAnalyzer.jetTag = 'jpticPu5patJets'

process.inclusiveJetAnalyzerSequence = cms.Sequence(                  #process.icPu5JetAnalyzer
                                                                      #*process.akPu5PFJetAnalyzer
                                                                      #*process.akPu4PFJetAnalyzer
                                                                      #*process.akPu3PFJetAnalyzer
                                                                      process.icPu5JPTJetAnalyzer
                                                                      )

process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_cff")
process.PFJetAnalyzer.isMC = cms.untracked.bool(False)
if useGoodTightCollection:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodTightTracksDirect")
elif useHighPtTrackCollection:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiHighPtTracks")
else:
    process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodMergedTracks")
process.PFJetAnalyzer.jetTag2 = cms.InputTag("akPu5PFpatJets")
process.PFJetAnalyzer.recoJetTag2 = cms.InputTag("akPu5PFJets")
process.PFJetAnalyzer.jetTag3 = cms.InputTag("akPu4PFpatJets")
process.PFJetAnalyzer.recoJetTag3 = cms.InputTag("akPu4PFJets")
process.PFJetAnalyzer.jetTag4 = cms.InputTag("akPu3PFpatJets")
process.PFJetAnalyzer.recoJetTag4 = cms.InputTag("akPu3PFJets")

#Frank's analyzer

process.load("PbPbTrackingTools.HiTrackValidator.hitrackvalidator_cfi")

if useGoodTightCollection:
    process.hitrkvalidator.trklabel=cms.untracked.InputTag("hiGoodTightTracksDirect")
elif useHighPtTrackCollection:
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

from Saved.QM11Ana.Analyzers_cff import trkAnalyzer
process.trkAnalyzer = trkAnalyzer

if useGoodTightCollection:
    process.trkAnalyzer.trackSrc = cms.InputTag("hiGoodTightTracksDirect")
elif useHighPtTrackCollection:
    process.trkAnalyzer.trackSrc = cms.InputTag("hiHighPtTracks")
else:
    process.trkAnalyzer.trackSrc = cms.InputTag("hiGoodTracks")


#Yetkin's analyzer

from CmsHi.JetAnalysis.ak3PFTowerJetsAna_cff import *

process.yetkinsAna = akPu3PFtowerJetsAna.clone()
process.RandomNumberGeneratorService.yetkinsAna = process.RandomNumberGeneratorService.generator.clone(initialSeed = 1)

# track efficiency anlayzer
#process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")
#for tree output
process.TFileService = cms.Service("TFileService",
                                   #fileName=cms.string("JetAnalysisTTrees_hiGoodTracks_v2.root"))
                                   #fileName=cms.string("JetAnalysisTTrees_hiHighPtTracks_v2.root"))
                                   fileName=cms.string("JetAnalysisTTrees_hiGoodTightTracks_v1.root"))



process.load("MNguyen.Configuration.HI_JetSkim_cff")
process.hltJetHI.HLTPaths = ["HLT_HIJet35U"]

if useGoodTightCollection:
    process.jetSkimPath = cms.Path(
        process.jetSkimSequence*
        process.hiTrackReco*
        process.hiGoodTracksSelection*
        process.hiGoodTightTracksDirectSelection*
        process.hiHighPtTrackDirectSelection*
        process.muonRecoPbPb*
        process.HiParticleFlowRecoNoJets*
        process.PFTowers*
        process.runAllJets*
        process.inclusiveJetAnalyzerSequence*
        process.PFJetAnalyzerSequence
        *process.trkAnalyzer
        *process.hiValidatorSequence
        *process.yetkinsAna
        )
elif useHighPtTrackCollection:
    process.jetSkimPath = cms.Path(
        process.jetSkimSequence*
        process.hiTrackReco*
        process.hiHighPtTrackSelection*
        process.muonRecoPbPb*
        process.HiParticleFlowRecoNoJets*
        process.PFTowers*
        process.runAllJets*
        process.inclusiveJetAnalyzerSequence*
        process.PFJetAnalyzerSequence
        *process.trkAnalyzer
        *process.hitrkvalidator
        *process.yetkinsAna
        )

else:
    process.jetSkimPath = cms.Path(
        process.jetSkimSequence*
        process.hiTrackReco*
        process.hiGoodTracksSelection*
        #process.conformalPixelTrackReco *
        #process.hiGoodMergedTracks*
        process.muonRecoPbPb*
        process.HiParticleFlowRecoNoJets*
        process.PFTowers*
        process.runAllJets*
        process.inclusiveJetAnalyzerSequence*
        process.PFJetAnalyzerSequence
        *process.trkAnalyzer
        *process.hitrkvalidator
        *process.yetkinsAna
        )

#process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")
#
#process.output = cms.OutputModule("PoolOutputModule",
#                                  process.jetTrkSkimContent,
#                                  SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('jetSkimPath')
#    ),
#                                  fileName = cms.untracked.string("PAT.root"),
#                                  dataset = cms.untracked.PSet(
#    filterName = cms.untracked.string('jetSkimPath'),
#    dataTier = cms.untracked.string('PAT')
#    )
#                                  )
#
#
## Save some extra stuff
##?
#process.output.outputCommands.extend(["keep *_heavyIon_*_*"])
## triggger
#process.output.outputCommands.extend(["keep *_TriggerResults_*_*"])
##tracks
#process.output.outputCommands.extend(["keep *_hiGoodMergedTracks_*_*"])
## reco jets
#process.output.outputCommands.extend(["keep recoCaloJets_*_*_*"])
#process.output.outputCommands.extend(["keep recoPFJets_*_*_*"])
##particle flow
#process.output.outputCommands.extend(["keep recoPFCandidates_particleFlow_*_*"])
##process.output.outputCommands.extend(["keep recoPFClusters_*_*_*"])
##process.output.outputCommands.extend(["keep recoPFRecHits_*_*_*"])
##fast jet pf stuff
#process.output.outputCommands.extend(["keep doubles_*PF*_*_*"])
#process.output.outputCommands.extend(["keep *_PFTowers_*_*"])
##calorimeter stuff
#process.output.outputCommands.extend(["keep *_towerMaker_*_*"])
#process.output.outputCommands.extend(["keep *_caloTowers_*_*"])
#process.output.outputCommands.extend(["keep *_hcalnoise_*_*"])
##process.output.outputCommands.extend(["keep *_hbhereco_*_*"])
##process.output.outputCommands.extend(["keep *_horeco_*_*"])
##process.output.outputCommands.extend(["keep *_hfreco_*_*"])
##process.output.outputCommands.extend(["keep *_ecalRecHit_*_*"])
#
##JPT
#process.output.outputCommands.extend(["keep *_jptic*_*_*"])
#process.output.outputCommands.extend(["keep *_recoJPT*_*_*"])
#process.output.outputCommands.extend(["keep *_JetPlusTrack*_*_*"])
#
#
#
## Save all RECO!
##process.output.outputCommands.extend(["keep *_*_*_RECO"])
#
#process.out_step = cms.EndPath(process.output)




# Schedule definition
#process.schedule = cms.Schedule(process.jetSkimPath,process.out_step)


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

process.MessageLogger.cerr.FwkReport.reportEvery = 1









