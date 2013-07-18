import FWCore.ParameterSet.Config as cms

highPurityTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string(
    'quality("highPurity") && (ptError/pt)<0.20 && algo<=7 && numberOfValidHits>=5')
)

from edwenger.TrkEffAnalyzer.selectTracks_cfi import *
selectTracks.src = cms.InputTag("highPurityTracks")
selectTracks.d0_par2 = cms.vdouble(3, 0.0)
selectTracks.dz_par2 = cms.vdouble(3, 0.0)
selectTracks.max_d0 = cms.double(0.2)
selectTracks.max_z0 = cms.double(0.2)
selectTracks.vertices = cms.InputTag("sortedGoodVertices")

preSelectTracks = selectTracks.clone(
    vertices = cms.InputTag("sortedPreGoodVertices")
    )

looseSelectTracks = selectTracks.clone(
    vertices = cms.InputTag("selectedVertex")
    )

allTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("selectTracks"),
    particleType = cms.string('pi+')
)

goodTracks = cms.EDFilter("PtMinCandViewSelector",
    src = cms.InputTag("allTracks"),
    ptMin = cms.double(0.0)
)

highMultFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("allTracks"),
    minNumber = cms.uint32(80)
)

chargedCandidates = cms.Sequence(highPurityTracks *
                                 preSelectTracks*
                                 selectTracks *
                                 looseSelectTracks *
                                 allTracks *
                                 goodTracks)
