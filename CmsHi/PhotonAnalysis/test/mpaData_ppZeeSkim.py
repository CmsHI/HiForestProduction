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
    ___inf___
    # lumisToProcess =  cms.untracked.VLuminosityBlockRange(
    # '150431:1-150431:1000'
    )
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
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)

# Modification for HI
process.load("CmsHi.PhotonAnalysis.MultiPhotonAnalyzer_cfi")
process.multiPhotonAnalyzer.OutputFile = cms.string('___outf1___')
process.multiPhotonAnalyzer.doStoreCompCone = cms.untracked.bool(True)
process.multiPhotonAnalyzer.doStoreJets = cms.untracked.bool(False)

# HiGoodTrack
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")    #process.trksel_step  = cms.Path(process.hiGoodTracksSelection)

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
#process.patHeavyIonDefaultSequence.remove(process.hiEgammaIsolationSequence)

 
# the path! 
process.p = cms.Path(
    #  process.evtCounter *
    #    process.hiGoodTracksSelection * 
    process.hiEcalClusteringSequence * # for pp only
    process.hiPhotonCleaningSequence *
    process.compleCleanPhotonSequence *
    process.patHeavyIonDefaultSequence *
    process.multiPhotonAnalyzer 
    )

# for pp data
process.multiPhotonAnalyzer.VertexProducer = cms.InputTag("offlinePrimaryVertices")
process.multiPhotonAnalyzer.doStoreCentrality = False
process.multiPhotonAnalyzer.TrackProducer     = cms.InputTag('generalTracks')

process.cleanPhotons.primaryVertexProducer = cms.InputTag('offlinePrimaryVerticesWithBS')
process.cleanPhotons.isolationSumsCalculatorSet.trackProducer = cms.InputTag('generalTracks')
process.compleCleanPhoton.primaryVertexProducer = cms.InputTag('offlinePrimaryVerticesWithBS')
process.compleCleanPhoton.isolationSumsCalculatorSet.trackProducer = cms.InputTag('generalTracks')

process.gamIsoDepositTk.ExtractorPSet.inputTrackCollection   =  cms.InputTag("generalTracks")
process.isolationInputParameters.track = cms.InputTag("generalTracks")

process.isoCC1.track = cms.InputTag('generalTracks')
process.isoCC2.track = cms.InputTag('generalTracks')
process.isoCC3.track = cms.InputTag('generalTracks')
process.isoCC4.track = cms.InputTag('generalTracks')
process.isoCC5.track = cms.InputTag('generalTracks')

process.isoCR1.track = cms.InputTag('generalTracks')
process.isoCR2.track = cms.InputTag('generalTracks')
process.isoCR3.track = cms.InputTag('generalTracks')
process.isoCR4.track = cms.InputTag('generalTracks')
process.isoCR5.track = cms.InputTag('generalTracks')

process.isoT11.track = cms.InputTag('generalTracks')
process.isoT21.track = cms.InputTag('generalTracks')
process.isoT31.track = cms.InputTag('generalTracks')
process.isoT41.track = cms.InputTag('generalTracks')
process.isoT12.track = cms.InputTag('generalTracks')
process.isoT22.track = cms.InputTag('generalTracks')
process.isoT32.track = cms.InputTag('generalTracks')
process.isoT42.track = cms.InputTag('generalTracks')
process.isoT13.track = cms.InputTag('generalTracks')
process.isoT23.track = cms.InputTag('generalTracks')
process.isoT33.track = cms.InputTag('generalTracks')
process.isoT43.track = cms.InputTag('generalTracks')
process.isoT14.track = cms.InputTag('generalTracks')
process.isoT24.track = cms.InputTag('generalTracks')
process.isoT34.track = cms.InputTag('generalTracks')
process.isoT44.track = cms.InputTag('generalTracks')

process.isoDR11.track = cms.InputTag('generalTracks')
process.isoDR21.track = cms.InputTag('generalTracks')
process.isoDR31.track = cms.InputTag('generalTracks')
process.isoDR41.track = cms.InputTag('generalTracks')
process.isoDR12.track = cms.InputTag('generalTracks')
process.isoDR22.track = cms.InputTag('generalTracks')
process.isoDR32.track = cms.InputTag('generalTracks')
process.isoDR42.track = cms.InputTag('generalTracks')
process.isoDR13.track = cms.InputTag('generalTracks')
process.isoDR23.track = cms.InputTag('generalTracks')
process.isoDR33.track = cms.InputTag('generalTracks')
process.isoDR43.track = cms.InputTag('generalTracks')
process.isoDR14.track = cms.InputTag('generalTracks')
process.isoDR24.track = cms.InputTag('generalTracks')
process.isoDR34.track = cms.InputTag('generalTracks')
process.isoDR44.track = cms.InputTag('generalTracks')


