import FWCore.ParameterSet.Config as cms

evtselanalyzer = cms.EDAnalyzer('EvtSelAnalyzer',
      gtlabel=cms.untracked.InputTag('gtDigis'),
      towerslabel=cms.untracked.InputTag('towerMaker'),
      hfEThreshold=cms.untracked.double(3.0),
      trackslabel=cms.untracked.InputTag('generalTracks'),
      qualityString=cms.untracked.string('highPurity'),
      triglabel=cms.untracked.InputTag('TriggerResults','','HLT'),
      trignames=cms.untracked.vstring('HLT_L1Tech_BSC_minBias_OR',
                                      'HLT_MinBiasPixel_SingleTrack',
                                      'HLT_L1Jet6U',
                                      'HLT_Jet15U',
                                      'HLT_Jet30U',
                                      'HLT_Jet50U'),
      isGEN=cms.untracked.bool(False),
      includeSelTrk=cms.untracked.bool(False),
      etaMaxSpec=cms.untracked.double(1.0),
      ptMin=cms.untracked.double(0.0),
      vtxlabel=cms.untracked.InputTag("pixel3Vertices"),
      vtxWeight=cms.untracked.bool(False),
      pvtxG_DATA=cms.untracked.vdouble(7.31011e-01,3.04673e+00,6.54385e-02),
      pvtxG_MC=cms.untracked.vdouble(-6.53244e-01,4.47771e+00,4.45992e-02)
)
