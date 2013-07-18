import FWCore.ParameterSet.Config as cms

# clone event selection analyzer to run after successive steps
from edwenger.EvtSelAnalyzer.evtselanalyzer_cfi import *
preTrgAna = evtselanalyzer.clone()
postTrgAna = evtselanalyzer.clone()
postEvtSelAna = evtselanalyzer.clone()
postVtxAna = evtselanalyzer.clone(includeSelTrk=cms.untracked.bool(True),
                                  vtxlabel=cms.untracked.InputTag("selectedVertex"))
postTrkVtxAna = evtselanalyzer.clone(includeSelTrk=cms.untracked.bool(True),
                                     vtxlabel=cms.untracked.InputTag("sortedGoodVertices"))

# clone vertex analyzer to run with various inputs
from edwenger.VertexAnalyzer.vertexanalyzer_cfi import *
preVtxSel = vertexanalyzer.clone(
    vtxlabel = cms.untracked.InputTag("pixel3Vertices"))
postVtxSel = vertexanalyzer.clone(
    vtxlabel = cms.untracked.InputTag("selectedVertex"))
preTrkVtxSel = vertexanalyzer.clone(
    vtxlabel = cms.untracked.InputTag("offlinePrimaryVertices"))
postTrkVtxSel = vertexanalyzer.clone(
    vtxlabel = cms.untracked.InputTag("sortedGoodVertices"))

# this is for studing trk-jet-vtx correlation
postExtraTrkVtxSel = vertexanalyzer.clone(
    vtxlabel = cms.untracked.InputTag("sortedPreGoodVertices"),
    trklabel = cms.untracked.InputTag("preSelectTracks"),
    jetlabel=cms.untracked.InputTag("cleanedPatJets"),
    jetTrkVerticesCorr=cms.untracked.bool(True)) 

postExtraTrkVtxSel_jet60 = postExtraTrkVtxSel.clone(minJetEt=cms.untracked.double(60))
postExtraTrkVtxSel_jet120 = postExtraTrkVtxSel.clone(minJetEt=cms.untracked.double(120))

# track and jet spectra analyzer
from edwenger.TrackSpectraAnalyzer.trackspectraanalyzer_cfi import *
trackAna.src = cms.untracked.InputTag("selectTracks")
trackAna.jsrc = cms.untracked.InputTag("cleanedPatJets")
trackAna.vsrc = cms.untracked.InputTag("sortedGoodVertices")

from edwenger.Skims.customise_cfi import getEvtSelEff
preTrackAna = trackAna.clone(isGEN=cms.untracked.bool(True),
                             pureGENmode=cms.untracked.bool(True),
                             #pthatCut = cms.untracked.double(15.0), #needed for 0-15 dijet sample
                             nsdOnly=cms.untracked.bool(True))

refitTrackAna = trackAna.clone(src=cms.untracked.InputTag("refitTracks"))

trackAna_STD = trackAna.clone(src = cms.untracked.InputTag("selectTracks"),
                              applyEvtEffCorr=cms.untracked.bool(True),
                              evtEffCorrType=cms.untracked.int32(0),  # STD
                              evtMultCut=cms.untracked.int32(3),
                              evtSelEffv = getEvtSelEff('STD_NSD_TrkVtx'))

looseTrackAna = trackAna.clone(src=cms.untracked.InputTag("looseSelectTracks"),
                               vsrc=cms.untracked.InputTag("selectedVertex"))

looseTrackAna_STD = looseTrackAna.clone(applyEvtEffCorr=cms.untracked.bool(True),
                                        evtEffCorrType=cms.untracked.int32(0), # STD
                                        evtMultCut=cms.untracked.int32(0),
                                        #evtMultCut=cms.untracked.int32(3),
				        evtSelEffv = getEvtSelEff('STD_NSD_PixVtx'))

                                        
# tracking efficiency analyzer
from edwenger.TrkEffAnalyzer.trkEffAnalyzer_cff import *
trkEffAnalyzer.tracks = cms.untracked.InputTag("selectTracks")
trkEffAnalyzer.vertices = cms.untracked.InputTag("sortedGoodVertices")
trkEffAnalyzer.fillNtuples = cms.bool(False)
trkEffAnalyzer.constPtBins = cms.bool(False)

loosetrkEffAnalyzer = trkEffAnalyzer.clone(tracks=cms.untracked.InputTag("looseSelectTracks"),
                                           vertices = cms.untracked.InputTag('selectedVertex',''),
                                           )

trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("selectTracks")

# sequences
preTrgTest = cms.Sequence(preTrgAna)
postTrgTest = cms.Sequence(postTrgAna)
postEvtSelTest = cms.Sequence(postEvtSelAna * preVtxSel)
postVtxTest = cms.Sequence(postVtxAna * postVtxSel)
preTrkVtxTest = cms.Sequence(preTrkVtxSel)
postTrkVtxTest = cms.Sequence(postTrkVtxAna * postTrkVtxSel)
