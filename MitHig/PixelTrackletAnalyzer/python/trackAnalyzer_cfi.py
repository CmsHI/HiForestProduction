import FWCore.ParameterSet.Config as cms

anaTrack = cms.EDAnalyzer('TrackAnalyzer',
                          trackPtMin = cms.untracked.double(0.4),
                          simTrackPtMin = cms.untracked.double(0.4),
                          vertexSrc = cms.vstring('hiSelectedVertex'),
                          trackSrc = cms.InputTag('hiGoodTightTracks'),
                          pfCandSrc = cms.InputTag('particleFlow'),
			  beamSpotSrc = cms.untracked.InputTag('offlineBeamSpot'),
                          doPFMatching = cms.untracked.bool(False),
                          doSimTrack = cms.untracked.bool(False),
                          doSimVertex = cms.untracked.bool(False),                          
                          useQuality = cms.untracked.bool(False),
                          qualityString = cms.untracked.string('highPurity'),
                          tpFakeSrc = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),
                          tpEffSrc = cms.untracked.InputTag('mergedtruth','MergedTrackTruth')
                          )
