
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.Timing = cms.Service("Timing")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_4/RelValTTbar/GEN-SIM-RECO/MC_39Y_V6-v1/0000/06C48E3D-EAF7-DF11-A9F3-002354EF3BDE.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_4/RelValTTbar/GEN-SIM-RECO/MC_39Y_V6-v1/0000/BCC068A2-F8F7-DF11-A156-00304867903E.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_4/RelValTTbar/GEN-SIM-RECO/MC_39Y_V6-v1/0000/CC64E46A-E7F7-DF11-8806-0026189438D8.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_9_4/RelValTTbar/GEN-SIM-RECO/MC_39Y_V6-v1/0000/E894061C-E1F7-DF11-ADCF-0030486792B4.root'
    ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.GlobalTag.globaltag = cms.string('MC_39Y_V6::All')

process.load("Configuration.StandardSequences.MagneticField_cff")


process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('isoRecHitCounter.root'),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.genParticleCounter = cms.EDAnalyzer("GenParticleCounter",
                                            src = cms.untracked.string("genParticles"),
                                            doCentrality =  cms.untracked.bool(False)
                                            )



process.p = cms.Path(
    process.genParticles *
    process.genParticleCounter
    )



    
