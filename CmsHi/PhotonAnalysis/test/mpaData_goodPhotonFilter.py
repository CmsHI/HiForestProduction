
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
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_9/RelValPyquen_GammaJet_pt20_2760GeV/GEN-SIM-RECO/START39_V7HI-v1/0000/A83BF2A6-F93D-E011-9DFB-002354EF3BDB.root'
    
    ),
                            inputCommands = cms.untracked.vstring(
    'keep *',
    'drop recoSuperClusters_islandSuperClusters_*_*',
    #    'drop recoPhotons_*_*_*',
    'drop recoPhotonCores_*_*_*',
    'drop recoCaloClusters_*_*_*',
    'drop recoGsfElectrons_*_*_*',
    'drop recoGsfElectronCores_*_*_*'
    ),
                            dropDescendantsOfDroppedBranches = cms.untracked.bool( False )
                            )

# USE JASON FILE For local jobs
#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#import FWCore.ParameterSet.Types as CfgTypes
#myLumis = LumiList.LumiList(filename = 'Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt').getCMSSWString().split(',')
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('___TFoutf___.root'),
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
removeMCMatching(process, ['Photons'])
removeMCMatching(process, ['Jets'])
removeMCMatching(process, ['Muons'])
process.patHeavyIonDefaultSequence.remove(process.heavyIon)
process.patHeavyIonDefaultSequence.remove(process.makeHeavyIonJets)
process.patHeavyIonDefaultSequence.remove(process.makeHeavyIonMuons)
process.patHeavyIonDefaultSequence.remove(process.selectedPatJets)
process.patHeavyIonDefaultSequence.remove(process.selectedPatMuons)

# trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
# =============== Trigger selection ====================
process.HIphotonTrig = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = cms.vstring('HLT_HIPhoton15')
    )

# clean collision selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")


# photon stkim!!!!
# photon selection

########### electron reco with island ##########
process.load("RecoEgamma.EgammaElectronProducers.ecalDrivenElectronSeedsParameters_cff")
process.load("RecoHI.HiEgammaAlgos.HiElectronSequence_cff")                 # gsf electrons
process.ecalDrivenElectronSeeds.barrelSuperClusters = cms.InputTag("hiSpikeCleanedSC")
process.ecalDrivenElectronSeeds.endcapSuperClusters = cms.InputTag("correctedIslandEndcapSuperClusters")
# chqnge the hoe to 0 22
# change the Et cut 13 GeV
process.ecalDrivenElectronSeeds.SeedConfiguration.maxHOverEBarrel = cms.double(0.25)
process.ecalDrivenElectronSeeds.SeedConfiguration.maxHOverEEndcaps = cms.double(0.01)
process.ecalDrivenElectronSeeds.SeedConfiguration.SCEtCut  = cms.double(13.0)
process.gsfElectrons.ctfTracks = cms.InputTag("hiGlobalPrimTracks")

########### random number seed
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   multiPhotonAnalyzer = cms.PSet(
    engineName = cms.untracked.string("TRandom3"),
    initialSeed = cms.untracked.uint32(98236)
    )
                                                   )

############ Make HiGoodTrack tight!!
process.hiGoodTracks.src = cms.InputTag("hiPostGlobalPrimTracks")
process.hiGoodTracks.keepAllTracks = cms.bool(False)
process.hiGoodTracks.qualityBit = cms.string('highPurity')
process.hiGoodTracks.min_nhits = cms.uint32(13)
process.hiGoodTracks.chi2n_par = cms.double(0.15)

################ good photon filter #######
process.goodPhotons = cms.EDFilter("PhotonSelector",
                                   src = cms.InputTag("cleanPhotons"),
                                   cut = cms.string("pt > 18 & hadronicOverEm < 0.2 & abs(eta) < 1.45"),
                                   )
process.filterGoodPhotons = cms.EDFilter("PhotonCountFilter",
                                         src = cms.InputTag("goodPhotons"),
                                         minNumber = cms.uint32(1),
                                         maxNumber = cms.uint32(999999),
                                         )

process.barrelPhotonFilter = cms.Sequence ( process.goodPhotons * process.filterGoodPhotons )
 

# the path! 
process.p = cms.Path(
#    process.HIphotonTrig *  All Physics!!
    process.collisionEventSelection *
    process.hiEcalClusteringSequence *
    process.hiPhotonCleaningSequence *
    process.barrelPhotonFilter *
    process.hiPostGlobalPrimTracks * process.hiGoodTracks *
    process.siPixelRecHits*
    process.hiPixel3PrimTracks*
    process.hiPixelTrackSeeds*
    process.hiElectronSequence*
    process.compleCleanPhotonSequence *
    process.patHeavyIonDefaultSequence *
    process.multiPhotonAnalyzer 
    )
