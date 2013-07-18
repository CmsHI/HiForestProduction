
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## Magnet, geometry, detector condition ( needed for patTuple )
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("RecoHI.HiEgammaAlgos.HiEgamma_cff")
## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# Centrality
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.HeavyIonGlobalParameters = cms.PSet(
    centralitySrc = cms.InputTag("hiCentrality"),
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("")
    )
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_39Y_V3::All')  # for data global run.
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    ___inf___
    # lumisToProcess =  cms.untracked.VLuminosityBlockRange(
    # '150431:1-150431:1000'
    ),
                            inputCommands = cms.untracked.vstring(
    'keep *',
    'drop recoSuperClusters_islandSuperClusters_*_*',
    'drop recoPhotons_*_*_*',
    'drop recoPhotonCores_*_*_*',
    'drop recoCaloClusters_*_*_*'
    ),
                            dropDescendantsOfDroppedBranches = cms.untracked.bool( False )
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('___TFoutf___'),
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
process.multiPhotonAnalyzer.doStoreMET = cms.untracked.bool(False)
process.multiPhotonAnalyzer.doStoreJets = cms.untracked.bool(False)
process.multiPhotonAnalyzer.OutputFile = cms.string('___outf___')
process.multiPhotonAnalyzer.doStoreCompCone = cms.untracked.bool(True)

# HiGoodTrack
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")    #process.trksel_step  = cms.Path(process.hiGoodTracksSelection)

# detector responce
process.load("CmsHi.PhotonAnalysis.isoConeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.ecalHistProducer_cfi")
process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")

# spike cleaner for island superclsters
process.cleanPhotons.maxHoverEBarrel = cms.double(1)
process.cleanPhotons.maxHoverEBarrel = cms.double(1)
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

# trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.HIphotontrig = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.HIphotontrig.HLTPaths = ["HLT_HIPhoton15"]
process.HIphotontrig.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

# clean collision selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

############ Changes for data mixed evenets ################
# make hi gen particle
#process.load("PhysicsTools.HepMCCandAlgos.HiGenParticles_cfi")
# remove some parts in pat
process.patHeavyIonDefaultSequence.remove(process.heavyIonCleanedGenJets)
process.patHeavyIonDefaultSequence.remove(process.patJetGenJetMatch)
process.patHeavyIonDefaultSequence.remove(process.patJetGenJetMatch)
# change the genEventScale to calculate the pT hat 
process.multiPhotonAnalyzer.GenEventScale = cms.InputTag("hiSignal")
############################################################

####################### tune ##############
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
########################################################################

# the path! 
process.p = cms.Path(
    # process.HIphotontrig *
    #    process.collisionEventSelection *
 #   process.hiGenParticles * 
    process.hiGoodTracksSelection * # process.highPurityTracks *
    process.hiEcalClusteringSequence *
    process.hiPhotonCleaningSequence *
    process.patHeavyIonDefaultSequence *
    process.compleCleanPhotonSequence *
    process.multiPhotonAnalyzer 
    )
