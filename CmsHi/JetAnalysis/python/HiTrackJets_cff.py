import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5TrackJets_cfi import *

from PhysicsTools.RecoAlgos.TrackWithVertexRefSelector_cfi import *
from RecoJets.JetProducers.TracksForJets_cff import *

trackWithVertexRefSelector.src = cms.InputTag("hiTracks")
trackWithVertexRefSelector.ptErrorCut = cms.double(99999)
trackWithVertexRefSelector.d0Max = cms.double(99999)
trackWithVertexRefSelector.numberOfValidHits = cms.uint32(999999)
trackWithVertexRefSelector.useVtx = cms.bool(False)
trackWithVertexRefSelector.quality = cms.string("highPurity")
trackWithVertexRefSelector.vertexTag = cms.InputTag("offlinePrimaryVertices")

ak5TrackJets.doAreaFastjet   = cms.bool(True)
ak5TrackJets.doRhoFastjet    = cms.bool(True)
ak5TrackJets.doFastJetNonUniform = cms.bool(True)
ak5TrackJets.puCenters = cms.vdouble(-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
ak5TrackJets.puWidth = cms.double(0.5)
ak5TrackJets.nExclude = cms.uint32(0)
ak5TrackJets.Rho_EtaMax = cms.double(5.5)
ak5TrackJets.Ghost_EtaMax = cms.double(5.5)

kt5TrackJets = ak5TrackJets.clone(jetAlgorithm = cms.string("Kt"))

ak1TrackJets05 = ak5TrackJets.clone(rParam = cms.double(0.1),puWidth = 0.5)
ak2TrackJets05 = ak5TrackJets.clone(rParam = cms.double(0.2),puWidth = 0.5)
ak3TrackJets05 = ak5TrackJets.clone(rParam = cms.double(0.3),puWidth = 0.5)
ak4TrackJets05 = ak5TrackJets.clone(rParam = cms.double(0.4),puWidth = 0.5)
ak5TrackJets05 = ak5TrackJets.clone(rParam = cms.double(0.5),puWidth = 0.5)
ak6TrackJets05 = ak5TrackJets.clone(rParam = cms.double(0.6),puWidth = 0.5)

ak1TrackJets1 = ak5TrackJets.clone(rParam = cms.double(0.1),puWidth = 1)
ak2TrackJets1 = ak5TrackJets.clone(rParam = cms.double(0.2),puWidth = 1)
ak3TrackJets1 = ak5TrackJets.clone(rParam = cms.double(0.3),puWidth = 1)
ak4TrackJets1 = ak5TrackJets.clone(rParam = cms.double(0.4),puWidth = 1)
ak5TrackJets1 = ak5TrackJets.clone(rParam = cms.double(0.5),puWidth = 1)
ak6TrackJets1 = ak5TrackJets.clone(rParam = cms.double(0.6),puWidth = 1)

ak1TrackJets2 = ak5TrackJets.clone(rParam = cms.double(0.1),puWidth = 2)
ak2TrackJets2 = ak5TrackJets.clone(rParam = cms.double(0.2),puWidth = 2)
ak3TrackJets2 = ak5TrackJets.clone(rParam = cms.double(0.3),puWidth = 2)
ak4TrackJets2 = ak5TrackJets.clone(rParam = cms.double(0.4),puWidth = 2)
ak5TrackJets2 = ak5TrackJets.clone(rParam = cms.double(0.5),puWidth = 2)
ak6TrackJets2 = ak5TrackJets.clone(rParam = cms.double(0.6),puWidth = 2)

kt1TrackJets05 = kt5TrackJets.clone(rParam = cms.double(0.1),puWidth = 0.5)
kt2TrackJets05 = kt5TrackJets.clone(rParam = cms.double(0.2),puWidth = 0.5)
kt3TrackJets05 = kt5TrackJets.clone(rParam = cms.double(0.3),puWidth = 0.5)
kt4TrackJets05 = kt5TrackJets.clone(rParam = cms.double(0.4),puWidth = 0.5)
kt5TrackJets05 = kt5TrackJets.clone(rParam = cms.double(0.5),puWidth = 0.5)
kt6TrackJets05 = kt5TrackJets.clone(rParam = cms.double(0.6),puWidth = 0.5)

kt1TrackJets1 = kt5TrackJets.clone(rParam = cms.double(0.1),puWidth = 1)
kt2TrackJets1 = kt5TrackJets.clone(rParam = cms.double(0.2),puWidth = 1)
kt3TrackJets1 = kt5TrackJets.clone(rParam = cms.double(0.3),puWidth = 1)
kt4TrackJets1 = kt5TrackJets.clone(rParam = cms.double(0.4),puWidth = 1)
kt5TrackJets1 = kt5TrackJets.clone(rParam = cms.double(0.5),puWidth = 1)
kt6TrackJets1 = kt5TrackJets.clone(rParam = cms.double(0.6),puWidth = 1)

kt1TrackJets2 = kt5TrackJets.clone(rParam = cms.double(0.1),puWidth = 2)
kt2TrackJets2 = kt5TrackJets.clone(rParam = cms.double(0.2),puWidth = 2)
kt3TrackJets2 = kt5TrackJets.clone(rParam = cms.double(0.3),puWidth = 2)
kt4TrackJets2 = kt5TrackJets.clone(rParam = cms.double(0.4),puWidth = 2)
kt5TrackJets2 = kt5TrackJets.clone(rParam = cms.double(0.5),puWidth = 2)
kt6TrackJets2 = kt5TrackJets.clone(rParam = cms.double(0.6),puWidth = 2)


recoTrackJets = cms.Sequence(trackWithVertexRefSelector*trackRefsForJets
                             *
                             ak1TrackJets05+ak2TrackJets05+ak3TrackJets05+ak4TrackJets05+ak5TrackJets05+ak6TrackJets05
                             +
                             ak1TrackJets1+ak2TrackJets1+ak3TrackJets1+ak4TrackJets1+ak5TrackJets1+ak6TrackJets1
                             +
                             ak1TrackJets2+ak2TrackJets2+ak3TrackJets2+ak4TrackJets2+ak5TrackJets2+ak6TrackJets2
                             +
                             kt1TrackJets05+kt2TrackJets05+kt3TrackJets05+kt4TrackJets05+kt5TrackJets05+kt6TrackJets05
                             +
                             kt1TrackJets1+kt2TrackJets1+kt3TrackJets1+kt4TrackJets1+kt5TrackJets1+kt6TrackJets1
                             +
                             kt1TrackJets2+kt2TrackJets2+kt3TrackJets2+kt4TrackJets2+kt5TrackJets2+kt6TrackJets2
                             )











