import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("edwenger.EvtSelAnalyzer.evtselanalyzer_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/d101/edwenger/data7TeV/v3/trkAnaSkimAOD_1_2.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('test.root'))

process.p = cms.Path(process.evtselanalyzer)
