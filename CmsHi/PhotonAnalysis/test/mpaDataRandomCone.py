
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## Magnet, geometry, detector condition ( needed for patTuple )
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoHI.HiEgammaAlgos.HiEgamma_cff")
## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.Tracer = cms.Service("Tracer")
# Centrality
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff'),
process.HeavyIonGlobalParameters = cms.PSet(
    centralitySrc = cms.InputTag("hiCentrality"),
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("")
    )
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_39X_V6B::All')  # for data global run.
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/hidata/HIRun2010/HIAllPhysics/RECO/SDmaker_3SD_1CS_PDHIAllPhysicsZSv2_SD_PhotonHI-v1/0048/000E18FD-B24C-E011-8820-003048C9C1D4.root'
    ),
                            inputCommands = cms.untracked.vstring(
    'keep *',
    #    'drop recoSuperClusters_*_*_*',
    'drop recoPhotons_*_*_*',
    'drop recoPhotonCores_*_*_*'
    #    'drop recoCaloClusters_*_*_*'
    ),
    dropDescendantsOfDroppedBranches = cms.untracked.bool( False )
    )
                            
# USE JASON FILE For local jobs
import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
myLumis = LumiList.LumiList(filename = 'Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)



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
process.multiPhotonAnalyzer.OutputFile = cms.string('___outf___')
process.multiPhotonAnalyzer.doStoreCompCone = cms.untracked.bool(True)
process.multiPhotonAnalyzer.doStoreJets = cms.untracked.bool(False)

############ Make HiGoodTrack tight!!
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")
process.hiGoodTracks.src = cms.InputTag("hiPostGlobalPrimTracks")
process.hiGoodTracks.keepAllTracks = cms.bool(False)
process.hiGoodTracks.qualityBit = cms.string('highPurity')
process.hiGoodTracks.min_nhits = cms.uint32(13)
process.hiGoodTracks.chi2n_par = cms.double(0.15)



# detector responce
process.load("CmsHi.PhotonAnalysis.isoConeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.ecalHistProducer_cfi")
process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")

# spike cleaner for island superclsters
process.cleanPhotons.maxHoverEBarrel = cms.double(100)
photonObj = "randomPhoton"
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
process.compleCleanSuperCluster.photonProducer= cms.untracked.InputTag(photonObj)


# random number generator for mpa
process.load('Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff')
process.RandomNumberGeneratorService.multiPhotonAnalyzer = cms.PSet(
        engineName = cms.untracked.string("TRandom3"),
                        initialSeed = cms.untracked.uint32(982346)
                        )

process.RandomNumberGeneratorService.randomSuperCluster = cms.PSet(
    engineName = cms.untracked.string("TRandom3"),
    initialSeed = cms.untracked.uint32(982346)
    )

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# turn off MC matching for data

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


# trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
# =============== Trigger selection ====================
process.HIphotonTrig = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = cms.vstring('HLT_HIPhoton15')
    )
process.HIminbiasTrig = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = cms.vstring('HLT_HIMinBiasHfOrBSC_Core')
    )

# clean collision selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

# random cone path! 
process.p = cms.Path(
    process.HIminbiasTrig *
    process.collisionEventSelection *
    #  process.siPixelRecHits*
    process.hiPostGlobalPrimTracks * process.hiGoodTracks *
    process.randomConeSqeunceWithCalo *
    process.patHeavyIonDefaultSequence *
    process.compleCleanPhotonSequence * 
    process.multiPhotonAnalyzer
    )       


process.isoCC1.photons = cms.InputTag(photonObj)
process.isoCC2.photons = cms.InputTag(photonObj)
process.isoCC3.photons = cms.InputTag(photonObj)
process.isoCC4.photons = cms.InputTag(photonObj)
process.isoCC5.photons = cms.InputTag(photonObj)
process.isoCR1.photons = cms.InputTag(photonObj)
process.isoCR2.photons = cms.InputTag(photonObj)
process.isoCR3.photons = cms.InputTag(photonObj)
process.isoCR4.photons = cms.InputTag(photonObj)
process.isoCR5.photons = cms.InputTag(photonObj)

process.isoT11.photons = cms.InputTag(photonObj)
process.isoT12.photons = cms.InputTag(photonObj)
process.isoT13.photons = cms.InputTag(photonObj)
process.isoT14.photons = cms.InputTag(photonObj)
process.isoT21.photons = cms.InputTag(photonObj)
process.isoT22.photons = cms.InputTag(photonObj)
process.isoT23.photons = cms.InputTag(photonObj)
process.isoT24.photons = cms.InputTag(photonObj)
process.isoT31.photons = cms.InputTag(photonObj)
process.isoT32.photons = cms.InputTag(photonObj)
process.isoT33.photons = cms.InputTag(photonObj)
process.isoT34.photons = cms.InputTag(photonObj)
process.isoT41.photons = cms.InputTag(photonObj)
process.isoT42.photons = cms.InputTag(photonObj)
process.isoT43.photons = cms.InputTag(photonObj)
process.isoT44.photons = cms.InputTag(photonObj)
process.isoDR11.photons = cms.InputTag(photonObj)
process.isoDR12.photons = cms.InputTag(photonObj)
process.isoDR13.photons = cms.InputTag(photonObj)
process.isoDR14.photons = cms.InputTag(photonObj)
process.isoDR21.photons = cms.InputTag(photonObj)
process.isoDR22.photons = cms.InputTag(photonObj)
process.isoDR23.photons = cms.InputTag(photonObj)
process.isoDR24.photons = cms.InputTag(photonObj)
process.isoDR31.photons = cms.InputTag(photonObj)
process.isoDR32.photons = cms.InputTag(photonObj)
process.isoDR33.photons = cms.InputTag(photonObj)
process.isoDR34.photons = cms.InputTag(photonObj)
process.isoDR41.photons = cms.InputTag(photonObj)
process.isoDR42.photons = cms.InputTag(photonObj)
process.isoDR43.photons = cms.InputTag(photonObj)
process.isoDR44.photons = cms.InputTag(photonObj)

