
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## Magnet, geometry, detector condition ( needed for patTuple )
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("RecoHI.HiEgammaAlgos.HiEgamma_cff")
## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.Tracer = cms.Service("Tracer")
# Centrality
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.HeavyIonGlobalParameters = cms.PSet(
    centralitySrc = cms.InputTag("hiCentrality"),
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("")
    )
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR10_P_V12::All')  # for data global run.
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/d101/kimy/recoFiles/HiPhoton15_HIMinBiasHF_skim_99_1_LVg-151027.root'
    # lumisToProcess =  cms.untracked.VLuminosityBlockRange(
    # '150431:1-150431:1000'
    ),
                            inputCommands = cms.untracked.vstring(
    'keep *',
    'drop recoSuperClusters_*_*_*',
    'drop recoPhotons_*_*_*',
    'drop recoPhotonCores_*_*_*')
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ecalRechit_photon15trig.root'),
                                   closeFileFast = cms.untracked.bool(True)
                                   )
########################################################################################

## Standard Heavy Ion PAT Configuration File
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)
from PhysicsTools.PatAlgos.tools.coreTools import *
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# Modification for HI
process.load("CmsHi.PhotonAnalysis.MultiPhotonAnalyzer_cfi")
process.multiPhotonAnalyzer.GenParticleProducer = cms.InputTag("hiGenParticles")
process.multiPhotonAnalyzer.PhotonProducer = cms.InputTag("selectedPatPhotons")
process.multiPhotonAnalyzer.VertexProducer = cms.InputTag("hiSelectedVertex")
process.multiPhotonAnalyzer.OutputFile = cms.string('mpa.root')
process.multiPhotonAnalyzer.doStoreCompCone = cms.untracked.bool(True)
process.multiPhotonAnalyzer.doStoreJets = cms.untracked.bool(False)

# HiGoodMergedTrack
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")    #process.trksel_step  = cms.Path(process.hiGoodTracksSelection)
#process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
#process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
#process.hiGoodMergTrackSequence = cms.Sequence(
#    process.hiGoodTracksSelection*
#    process.conformalPixelTrackReco *
#    process.hiGoodMergedTracks
#    )

# detector responce
process.load("CmsHi.PhotonAnalysis.isoConeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.ecalHistProducer_cfi")
process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")

# spike cleaner for island superclsters
process.cleanPhotons.maxHoverEBarrel = cms.double(100)
photonObj = "cleanPhotons"
process.patPhotons.photonSource = cms.InputTag(photonObj)
process.photonMatch.src = cms.InputTag(photonObj)
process.PhotonIDProd.photonProducer  = cms.string(photonObj)
process.gamIsoDepositTk.src = cms.InputTag(photonObj)
process.gamIsoDepositEcalFromHits.src = cms.InputTag(photonObj)
process.gamIsoDepositHcalFromTowers.src = cms.InputTag(photonObj)
process.gamIsoDepositHcalDepth1FromTowers.src = cms.InputTag(photonObj)
process.gamIsoDepositHcalDepth2FromTowers.src = cms.InputTag(photonObj)
from RecoHI.HiEgammaAlgos.HiCoreTools import *

# random Cone sequence
process.load("RandomConeAna.Configuration.randomConeSequence_cff")
process.multiPhotonAnalyzer.compPhotonProducer = cms.InputTag("compleCleanPhoton")

# turn off MC matching for data
removeMCMatching(process, ['Photons'])
removeMCMatching(process, ['Jets'])
removeMCMatching(process, ['Muons'])
process.patHeavyIonDefaultSequence.remove(process.heavyIon)
process.patHeavyIonDefaultSequence.remove(process.genPartons)
process.patHeavyIonDefaultSequence.remove(process.hiPartons)
process.patHeavyIonDefaultSequence.remove(process.heavyIonCleanedGenJets)
process.patHeavyIonDefaultSequence.remove(process.patJetGenJetMatch)
process.patHeavyIonDefaultSequence.remove(process.photonMatch)
process.patHeavyIonDefaultSequence.remove(process.patJetGenJetMatch)

# trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
# =============== Trigger selection ====================
process.HIminbiasTrig = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = cms.vstring('HLT_HIMinBiasHfOrBSC_Core')
    )

# clean collision selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

# the path! 
process.p = cms.Path(
    process.HIminbiasTrig * #    process.HIphotonTrig * 
    process.collisionEventSelection *
    process.hiGoodTracksSelection *  # process.hiGoodMergTrackSequence *
    process.hiEcalClusteringSequence *
    process.hiPhotonCleaningSequence *
    process.compleCleanPhotonSequence *
    process.patHeavyIonDefaultSequence *
    process.multiPhotonAnalyzer 
    )
