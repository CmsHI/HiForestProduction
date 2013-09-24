import FWCore.ParameterSet.Config as cms

spectra = cms.EDAnalyzer('QWSpectra'
    ,tracks_ = cms.untracked.InputTag('hiGoodTightTracks')
    ,vtxCollection_ = cms.InputTag('hiSelectedVertex')
    ,centrality_ = cms.InputTag('centralityBin')
    ,eventplane_ = cms.InputTag('hiEvtPlane', 'recoLevel', 'HiEvtPlane')
)
