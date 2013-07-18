import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("edwenger.VertexAnalyzer.vertexanalyzer_cfi")
process.load("edwenger.VertexAnalyzer.SelectedVertex_cfi")

# message logger
process.MessageLogger.categories = ['VertexAnalyzer']
process.MessageLogger.debugModules = ['*']
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    DEBUG = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    VertexAnalyzer = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/d101/edwenger/data7TeV/v3/trkAnaSkimAOD_1_2.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('test.root'))

process.vtx = cms.Path(process.selectedVertex * process.vertexanalyzer)

process.output = cms.OutputModule("PoolOutputModule",
       SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('vtx')),
       fileName = cms.untracked.string('edmtest.root'))

process.out = cms.EndPath(process.output)
