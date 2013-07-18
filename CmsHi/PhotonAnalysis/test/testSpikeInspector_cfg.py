import FWCore.ParameterSet.Config as cms

process = cms.Process("SpikeDist")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = 'START38_V3::All'
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('SpikeDist')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'INPUTFIRSTFILE'
    ),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

INPUTFILEEXTEND

process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('OUTPUTFILENAME')
        )

process.p = cms.Path(process.spikeInspector)
