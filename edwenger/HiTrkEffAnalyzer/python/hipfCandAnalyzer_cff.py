import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *
TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

#myTrackAssociatorByHitsRecoDenom = cms.Sequence(TrackAssociatorByHitsRecoDenom)

from edwenger.HiTrkEffAnalyzer.HiTPCuts_cff import *
from edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cfi import *

from edwenger.HiTrkEffAnalyzer.hipfCandAnalyzer_cfi import *
pfCandidateAnalyzer.SimTracks = cms.InputTag("cutsTPForFak")
pfCandidateAnalyzer.ptMin = cms.untracked.double(1.0)
pfCandidateAnalyzer.Tracks = cms.InputTag("hiHighPtTracks")
pfCandidateAnalyzer.isData = cms.untracked.bool(False)
pfCandidateAnalyzer.useQaulityStr = cms.untracked.bool(True)
pfCandidateAnalyzer.qualityString = cms.untracked.string("tight")

pfCandidateAnalyzer_higtight = pfCandidateAnalyzer.clone(Tracks = cms.InputTag("hiGoodTightTracks"),
                                                         qualityString = cms.untracked.string("highPurity")
                                                         )

pfCandidateAnalyzer_higloose = pfCandidateAnalyzer.clone(Tracks = cms.InputTag("hiGoodLooseTracks"),
                                                         qualityString = cms.untracked.string("loose")
                                                         )

hipfCandAnalyzer = cms.Sequence(pfCandidateAnalyzer)
hipfCandAnalyzer_higtight = cms.Sequence(pfCandidateAnalyzer_higtight)
hipfCandAnalyzer_higloose = cms.Sequence(pfCandidateAnalyzer_higloose)


pfCandidateAnalyzer_test = pfCandidateAnalyzer.clone(Tracks = cms.InputTag("hiHighPtCaloTracks"))
hipfCandAnalyzer_test = cms.Sequence(pfCandidateAnalyzer_test)
