import FWCore.ParameterSet.Config as cms

from edwenger.Skims.eventSelection_cff import *
from edwenger.Skims.ExtraVertex_cff import *
from edwenger.Skims.Analyzers_cff import *
from edwenger.Skims.TrackSelection_cff import *
from edwenger.Skims.RootpleProducer_cfi import * 
from edwenger.Skims.TrackRefit_cff import *

extraReco = cms.Sequence(extraVertex *         ## agglomerative pixel vertexing
                         postEvtSelTest *
                         selectedVertex *      ## most-populated (filters)
                         goodVertices *        ## all vertices pass ndof>4, etc.
                         sortedPreGoodVertices *
                         sortedGoodVertices *
                         sortedSumPtOrdVertices * ## vertices sorted by sum(pt*pt)
                         chargedCandidates *   ## selected tracks -> charged candidates
                         postVtxTest *
                         preTrkVtxTest *
                         primaryVertexFilter * ## non-fake, ndof>4, abs(z)<15
                         #trackRefit * 
                         postTrkVtxTest)

