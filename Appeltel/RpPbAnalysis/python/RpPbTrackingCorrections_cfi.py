import FWCore.ParameterSet.Config as cms

from SimTracker.TrackHistory.TrackClassifier_cff import *

pPbTrkCorr  = cms.EDAnalyzer('RpPbTrackingCorrections',
    trackClassifier,
    trackSrc = cms.InputTag("generalTracks"),
    vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
    tpEffSrc = cms.InputTag('mergedtruth','MergedTrackTruth'),
    tpFakSrc = cms.InputTag('mergedtruth','MergedTrackTruth'),
    associatorMap = cms.InputTag('tpRecoAssocGeneralTracks'),
    jetSrc = cms.InputTag('akPu3PFpatJets'),
    ptBins = cms.vdouble(
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.05, 1.1, 1.15, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0,
        20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 
        180.0, 300.0, 500.0
    ),
    etaBins = cms.vdouble( 
        -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
        0.4, 0.8, 1.2, 1.6, 2.0, 2.4
    ),
    vtxWeightParameters = cms.vdouble( 4.49636e-02, 1.36629e-01, 5.30010e+00,
                                       2.50170e-02, 4.59123e-01, 9.64888e+00 ),
    doVtxReweighting = cms.bool(False),
    occBins = cms.vdouble( 0., 100.),
    occByCentrality = cms.bool(True),
    occByNPixelTrk = cms.bool(False),
    occByJetEt = cms.bool(False),
    occByLeadingJetEt = cms.bool(False),
    jetEtaMax = cms.double(2.0),
    jetRadius = cms.double(0.3),
    selectSpecies = cms.bool(False),
    pdgIdList = cms.vint32(211,-211),
    selectHeavyFlavorDecays = cms.bool(False),
    selectNonHeavyFlavorDecays = cms.bool(False),
    applyJetCuts = cms.bool(False),
    applyVertexZCut = cms.bool(False),
    vertexZMax = cms.double(15.),
    invertJetCuts = cms.bool(False),
    jetEtMin = cms.double(80.),
    jetEtMax = cms.double(999.),
    applyTrackCuts = cms.bool(True),
    qualityString = cms.string("highPurity"),
    dxyErrMax = cms.double(3.0),
    dzErrMax = cms.double(3.0),
    ptErrMax = cms.double(0.1),
    fillGenHistos = cms.bool(False),
    fillTrkPerfHistos = cms.bool(False),
    doMomRes = cms.bool(False)
)
