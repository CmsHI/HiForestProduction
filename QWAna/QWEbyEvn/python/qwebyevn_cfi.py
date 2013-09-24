import FWCore.ParameterSet.Config as cms

qwebye = cms.EDAnalyzer('QWEbyEvn'
    ,tracks_ = cms.untracked.InputTag('hiGoodTightTracks')
    ,vtxCollection_ = cms.InputTag('hiSelectedVertex')
    ,centrality_ = cms.InputTag('centralityBin')
    ,fillntuple_ = cms.untracked.bool(False)
)
