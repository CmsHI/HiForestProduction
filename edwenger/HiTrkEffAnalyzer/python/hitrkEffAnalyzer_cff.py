import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *
TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

from edwenger.HiTrkEffAnalyzer.HiTPCuts_cff import *
from edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cfi import *

## full tracking
hitrkEffAnalyzer.tracks = cms.untracked.InputTag('hiGoodTightTracks')
hitrkEffAnalyzer.label_tp_effic = cms.untracked.InputTag("cutsTPForEff")
hitrkEffAnalyzer.label_tp_fake = cms.untracked.InputTag("cutsTPForFak")
hitrkEffAnalyzer.hasSimInfo = cms.untracked.bool(True) # without this no sim track info
hitrkEffAnalyzer.doAssociation = cms.untracked.bool(True)
hitrkEffAnalyzer.pixelMultMode =  cms.untracked.bool(False)
hitrkEffAnalyzer.useJetEtMode = cms.int32(2)
hitrkEffAnalyzer.useQaulityStr=cms.untracked.bool(False)
hitrkEffAnalyzer.qualityString = cms.untracked.string("highPurity")

hitrkEffAna = cms.Sequence(cutsTPForFak*
                           cutsTPForEff*
                           hitrkEffAnalyzer)

higloosetrkEffAnalyzer = hitrkEffAnalyzer.clone(tracks = cms.untracked.InputTag('hiGoodLooseTracks'),
                           fillNtuples = cms.bool(False),
                           qualityString = cms.untracked.string("loose")                                                
                         )

hihightrkEffAnalyzer = hitrkEffAnalyzer.clone(tracks = cms.untracked.InputTag('hiHighPtTracks'),
                      fillNtuples = cms.bool(False),
                      qualityString = cms.untracked.string("tight")                                              
                      )

higloosetrkEffAna = cms.Sequence(cutsTPForFak*
                                 cutsTPForEff*
                                 higloosetrkEffAnalyzer)

hihightrkEffAna = cms.Sequence(cutsTPForFak*
                               cutsTPForEff*
                               hihightrkEffAnalyzer)

hicalotrkEffAnalyzer = hitrkEffAnalyzer.clone(tracks = cms.untracked.InputTag('hiCaloCompTracks'),
                                              fillNtuples = cms.bool(False)
                                              )

hicalotrkEffAna = cms.Sequence(cutsTPForFak*
                               cutsTPForEff*
                               hicalotrkEffAnalyzer)

## pixel tracking
hipxltrkEffAnalyzer = hitrkEffAnalyzer.clone(tracks = cms.untracked.InputTag('hiConformalPixelTracks'), # hiPixel3PrimTracks
                                             label_tp_effic = cms.untracked.InputTag('cutsTPForEffPxl'),
                                             label_tp_fake = cms.untracked.InputTag('cutsTPForFakPxl'),
                                             fillNtuples = cms.bool(False) #True
                                             )

hipxltrkEffAna = cms.Sequence(cutsTPForFakPxl*
                              cutsTPForEffPxl*
                              hipxltrkEffAnalyzer)
