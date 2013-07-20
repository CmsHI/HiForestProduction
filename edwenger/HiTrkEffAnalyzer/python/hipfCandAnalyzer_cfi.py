import FWCore.ParameterSet.Config as cms

# PF candidate analyzer
pfCandidateAnalyzer = cms.EDAnalyzer("HiPFCandidateTrackAnalyzer",
               PFCandidates = cms.InputTag("particleFlow",""),
               Vertices = cms.InputTag("hiSelectedVertex"),
               verbose = cms.untracked.bool(False), ## print candidate info
               printBlocks = cms.untracked.bool(False), ## print block/element info
               ptMin = cms.untracked.double(3.0), ## of PF candidate
               SimTracks = cms.InputTag("mergedtruth","MergedTrackTruth"),
               Tracks = cms.InputTag("hiSelectedTracks"),
               hasSimInfo = cms.untracked.bool(False),
               funcCaloComp = cms.string("(x>10) * (1.2*pow(TMath::Abs(x-10),6.0/9)) * (2+1.4/(exp(-10*(x-9))-1))"),                                     
               applyTrkQCs = cms.untracked.bool(False),
               minHits = cms.untracked.double(0),
               maxPtErr = cms.untracked.double(99999.),
               maxD0 = cms.untracked.double(99999.),
               maxDZ = cms.untracked.double(99999.),
               maxD0Norm = cms.untracked.double(99999.),
               maxDZNorm = cms.untracked.double(99999.),
               pixelSeedOnly = cms.untracked.bool(False),
               prodNtuple = cms.untracked.bool(False),
               isData = cms.untracked.bool(True),
               neededCentBins = cms.untracked.vint32(0,1,3,11,19,39),
               useQaulityStr=cms.untracked.bool(False),
               qualityString=cms.untracked.string("highPurity")
)
