import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("edwenger.EvtSelAnalyzer.evtselanalyzer_cfi")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = 'MC_3XY_V26::All'
process.GlobalTag.globaltag = 'START3X_V26A::All'


# message logger
process.MessageLogger.categories = ['EvtSelAnalyzer']
process.MessageLogger.debugModules = ['*']
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    DEBUG = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    EvtSelAnalyzer = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/relval/CMSSW_3_5_7/RelValMinBias/GEN-SIM-RECO/MC_3XY_V26-v1/0013/98E230F5-7F49-DF11-BE20-001A928116FA.root'
    '/store/relval/CMSSW_3_5_7/RelValMinBias/GEN-SIM-RECO/START3X_V26-v1/0012/B8D5AE61-6949-DF11-879B-00261894393F.root'    
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('testMC.root'))

process.evtselanalyzer.isGEN = True
process.p = cms.Path(process.evtselanalyzer)
