
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.Timing = cms.Service("Timing")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    "/store/hidata/HIRun2010/HIAllPhysics/RECO/PromptReco-v3/000/151/935/029A130C-AEF5-DF11-9BC7-003048F1C58C.root"
    ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("RecoHI.HiEgammaAlgos.HiEgamma_cff")

#process.GlobalTag.globaltag = cms.string('START36_V10::All')
process.GlobalTag.globaltag = cms.string('GR10_P_V12::All')  # for data global run.

process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

from PhysicsTools.PatAlgos.tools.coreTools import *

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)



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

process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
process.filter_step = cms.Path(process.collisionEventSelection)

# minbias events
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMinBiasHFOrBSC = process.hltHighLevel.clone()
process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_Core"]


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('hltAna_spikeCleanIsland_5.root'),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

#process.load("RecoHI.HiEgammaAlgos.hiEcalSpikeFilter_cfi")

process.load("CmsHi.PhotonAnalysis.isoConeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.ecalHistProducer_cfi")
process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")

# Clean the supercluster
process.cleanedIslandBarrelSuperClusters  = cms.EDProducer("HiSpikeCleaner",
                                                           recHitProducerBarrel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                                           recHitProducerEndcap = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
                                                           originalSuperClusterProducer = cms.InputTag("correctedIslandBarrelSuperClusters"),
                                                           outputColl  = cms.string( "" ),
                                                           etCut          = cms.double(8),
                                                           TimingCut    = cms.untracked.double(9999999.0),
                                                           swissCutThr    = cms.untracked.double(0.95)
                                                           )

# Ajust it..
process.ecalHistProducer.superClusterBarrel = cms.InputTag("cleanedIslandBarrelSuperClusters")


process.isoConeMap.etCut = 5;  # For hydjet

process.p = cms.Path(
    process.hltMinBiasHFOrBSC *
    process.collisionEventSelection * 
    #  process.isoConeMap  *
    process.cleanedIslandBarrelSuperClusters  * 
    process.ecalHistProducer 
    #process.spikeInspector
    )


#process.e = cms.EndPath(process.out)    
    
    
