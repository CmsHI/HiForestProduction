import FWCore.ParameterSet.Config as cms

process = cms.Process("hltAnaMu")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')

process.GlobalTag.globaltag = 'GR_P_V43D::All'
#process.GlobalTag.globaltag = 'GR_R_44_V4::All'
#process.GlobalTag.globaltag = 'STARTHI44_V4::All'
#process.GlobalTag.globaltag = 'START44_V4::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
#_input_
"/store/hidata/HIRun2013/PAMuon/RECO/PromptReco-v1/000/211/607/00000/ECBCA1C0-B974-E211-A3C7-BCAEC5364CFB.root",
    )
)

process.TFileService = cms.Service("TFileService",
#    fileName = cms.string(_output_)
    fileName = cms.string("test.root")
)

# load centrality
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality"),
  pPbRunFlip = cms.untracked.uint32(211313)
  )

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

process.analysis = cms.EDAnalyzer('HLTMuTree',
  muons = cms.InputTag("muons"),
#  vertices = cms.InputTag("hiSelectedVertex"),
  vertices = cms.InputTag("offlinePrimaryVertices"),
  doReco = cms.untracked.bool(True),
  doGen = cms.untracked.bool(False),
  genparticle = cms.InputTag("genParticles"),
  simtrack = cms.InputTag("mergedtruth","MergedTrackTruth"),
)

process.p = cms.Path(process.siPixelRecHits * process.PAcollisionEventSelection * process.pACentrality * process.analysis)
