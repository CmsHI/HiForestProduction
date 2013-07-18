import FWCore.ParameterSet.Config as cms

# vertex constraint
doConstraint = cms.EDProducer("VertexConstraintProducer")

# refit tracks
from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter
refitTracks = TrackRefitter.clone(
    constraint = 'vertex',
    srcConstr = 'doConstraint'
    )

# sequence
trackRefit = cms.Sequence(doConstraint * refitTracks)
