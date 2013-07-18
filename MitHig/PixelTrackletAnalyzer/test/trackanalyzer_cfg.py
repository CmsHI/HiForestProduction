import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/mnt/hadoop/cms/store/user/davidlw/Hydjet_Bass_MinBias_2760GeV/Summer11-STARTHI44_V4-v1-GEN-SIM-RECO_v1/a3a91e311bdf37073736c6b7978adab1/hiReco_RAW2DIGI_RECO_1938_3_znQ.root'
    )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'STARTHI44_V7::All'

process.demo = cms.EDAnalyzer('TrackAnalyzer',
                          trackPtMin = cms.untracked.double(0.4),
                          simTrackPtMin = cms.untracked.double(0.4),
                          vertexSrc = cms.vstring('hiSelectedVertex'),
                          trackSrc = cms.InputTag('hiSelectedTracks'),
                          pfCandSrc = cms.InputTag('particleFlow'),
                          beamSpotSrc = cms.untracked.InputTag('offlineBeamSpot'),
                          doPFMatching = cms.untracked.bool(False),
                          doSimTrack = cms.untracked.bool(False),
                          useQuality = cms.untracked.bool(False)
                          )

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('testntuple.root')
)        

process.p = cms.Path(process.demo)
