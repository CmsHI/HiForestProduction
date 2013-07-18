
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
    ),
                            inputCommands = cms.untracked.vstring(
    'keep *',
    #    'drop recoPhotons_*_*_*',
    'drop recoPhotonCores_*_*_*',
    'drop recoGsfElectrons_*_*_*',
    'drop recoGsfElectronCores_*_*_*'
    ),                            
                            dropDescendantsOfDroppedBranches = cms.untracked.bool( False ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string('tfile_1.root'),
#                                   closeFileFast = cms.untracked.bool(True)
#                                   )
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

# HiGoodMergedTrack
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff") 

# detector responce
process.load("CmsHi.PhotonAnalysis.isoConeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.ecalHistProducer_cfi")
process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.evtCounter_cfi")

# spike cleaner for island superclsters
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
process.load("PhysicsTools.HepMCCandAlgos.HiGenParticles_cfi")
# remove some parts in pat
process.patHeavyIonDefaultSequence.remove(process.heavyIonCleanedGenJets)
process.patHeavyIonDefaultSequence.remove(process.patJetGenJetMatch)
process.patHeavyIonDefaultSequence.remove(process.patJetGenJetMatch)
# change the genEventScale to calculate the pT hat 
process.multiPhotonAnalyzer.GenEventScale = cms.InputTag("generator")
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


# random number generator for mpa
process.load('Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff')
process.RandomNumberGeneratorService.multiPhotonAnalyzer = cms.PSet(
        engineName = cms.untracked.string("TRandom3"),
            initialSeed = cms.untracked.uint32(982346)
            )
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()


process.genParticleCounter = cms.EDAnalyzer("GenParticleCounter",
                                            src = cms.untracked.string("hiGenParticles"),
                                            vertexProducer_ = cms.untracked.string("hiSelectedVertex"),
                                            doCentrality =  cms.untracked.bool(True)
                                            )

 ############
process.load('MNguyen.iterTracking.secondStep_triplet_cff')
process.load('MNguyen.iterTracking.thirdStep_pairs_cff')
process.load('MNguyen.iterTracking.MergeTrackCollectionsHI_cff')
process.newCombinedSeeds = cms.EDProducer("SeedCombiner",
                                          seedCollections = cms.VInputTag(cms.InputTag("hiThrdSeedFromPairs"), cms.InputTag("hiScndSeedFromTriplets")),
                                          clusterRemovalInfos = cms.VInputTag(cms.InputTag("hiThrdClusters"), cms.InputTag("hiScndClusters"))
                                          )
process.iterTracking = cms.Sequence(process.hiGlobalPrimTracks * process.hiGoodTightTracks * process.secondStep * process.thirdStep * process.iterTracks * process.newCombinedSeeds)
#hiThrdTrackCandidates*hiThrdGlobalPrimTracks*ppThrdLoose*ppThrdTight*ppThrdTracks*hiThrdGoodTightTracks
process.iterTracking.remove(process.hiThrdTrackCandidates)
process.iterTracking.remove(process.hiThrdGlobalPrimTracks)
process.iterTracking.remove(process.ppThrdLoose)
process.iterTracking.remove(process.ppThrdTight)
process.iterTracking.remove(process.ppThrdTracks)
process.iterTracking.remove(process.hiThrdGoodTightTracks)
 
############ Make HiGoodTrack tight!!  this is the tight track but used for photon anlaysis
process.hiGoodTracks = process.hiGoodTightTracks.clone()
process.hiGoodTracks.src = cms.InputTag("hiPostGlobalPrimTracks")
process.hiGoodTracks.keepAllTracks = cms.bool(False)
process.hiGoodTracks.qualityBit = cms.string('highPurity')
process.hiGoodTracks.min_nhits = cms.uint32(13)
process.hiGoodTracks.chi2n_par = cms.double(0.15)


#### for electron seeds!!!!!###################################
process.hiPhotonCleaningSequence = cms.Sequence(process.hiSpikeCleanedSC*process.ecalDrivenElectronSeeds*process.cleanPhotonCore*process.cleanPhotons)
process.ecalDrivenElectronSeeds.SeedConfiguration.initialSeeds = cms.InputTag("newCombinedSeeds")
process.cleanPhotonCore.pixelSeedProducer = cms.string('ecalDrivenElectronSeeds')
process.hiElectronSequence.remove(process.ecalDrivenElectronSeeds)
#######################################

## pre filter ##
#cms.Sequence(hiSpikeCleanedSC*ecalDrivenElectronSeeds*cleanPhotonCore*cleanPhotons)
process.hiSC1 = process.hiSpikeCleanedSC.clone()
process.cleanPhotonCore1 = process.cleanPhotonCore.clone()
process.cleanPhotons1    = process.cleanPhotons.clone()
################ good photon filter #######
process.goodPhotons = cms.EDFilter("PhotonSelector",
                                   src = cms.InputTag("cleanPhotons1"),
                                   cut = cms.string("pt > 18 & hadronicOverEm < 0.23 & abs(eta) < 1.5"),
                                   )
process.filterGoodPhotons = cms.EDFilter("PhotonCountFilter",
                                         src = cms.InputTag("goodPhotons"),
                                         minNumber = cms.uint32(1),
                                         maxNumber = cms.uint32(999999),
                                         )
process.barrelPhotonFilter = cms.Sequence ( process.hiSC1 * process.cleanPhotonCore1 * process.cleanPhotons1 * process.goodPhotons * process.filterGoodPhotons )



# the path! 
process.p = cms.Path(
    # process.HIphotontrig *
    #    process.collisionEventSelection *
    process.barrelPhotonFilter *
    process.hiGenParticles *
    process.hiPostGlobalPrimTracks * process.hiGoodTracks *
    process.siPixelRecHits*
    #    process.hiPixel3PrimTracks*     process.hiPixelTrackSeeds*
    process.iterTracking *
    process.hiEcalClusteringSequence *
    process.hiPhotonCleaningSequence *
    process.hiElectronSequence*
    process.patHeavyIonDefaultSequence *
    process.compleCleanPhotonSequence *
    process.multiPhotonAnalyzer  #*
    #    process.genParticleCounter
    )
