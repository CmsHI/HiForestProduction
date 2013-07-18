import FWCore.ParameterSet.Config as cms


# Charged Candidates

allPixelTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiLowPtPixelTracks"), 
    particleType = cms.string('pi+')
)

allInvertedLowPtPixelTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiInvertedLowPtPixelTracks"),
    particleType = cms.string('pi+')
)

allSelectedTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiSelectedTracks"),
    particleType = cms.string('pi+')
)

allMergedNoPtSplitTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedNoPtSplitTracks"),
    particleType = cms.string('pi+')
)

allMergedPtSplit09Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit09Tracks"),
    particleType = cms.string('pi+')
)

allMergedPtSplit12Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit12Tracks"),
    particleType = cms.string('pi+')
)

allMergedPtSplit15Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit15Tracks"),
    particleType = cms.string('pi+')
)

allMergedPtSplit18Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit18Tracks"),
    particleType = cms.string('pi+')
)

hiAllTrackCandidates = cms.Sequence( allSelectedTracks
                                     * allPixelTracks
                                     * allInvertedLowPtPixelTracks
                                     * allMergedNoPtSplitTracks
                                     * allMergedPtSplit09Tracks
                                     * allMergedPtSplit12Tracks
                                     * allMergedPtSplit15Tracks
                                     * allMergedPtSplit18Tracks
)
