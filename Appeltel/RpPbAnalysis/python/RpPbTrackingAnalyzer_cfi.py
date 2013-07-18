import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

trkAna  = cms.EDAnalyzer('RpPbTrackingAnalyzer',
   trackSrc = cms.InputTag("generalTracks"),
   vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
   jetSrc = cms.InputTag('akPu3PFpatJets'),
   vertexZMax = cms.double(99.0),
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
    occBins = cms.vdouble( 0., 100.),
    occByCentrality = cms.bool(True),
    occByNPixelTrk = cms.bool(False),
    occByNcoll = cms.bool(False),
    occByLeadingJetEt = cms.bool(False),
    jetEtaMax = cms.double(2.0),
    applyCuts = cms.bool(False),
    qualityString = cms.string("highPurity"),
    dxyErrMax = cms.double(3.0),
    dzErrMax = cms.double(3.0),
    ptErrMax = cms.double(0.1),
    doMC = cms.bool(False),
    genSrc = cms.InputTag("genParticles"),
    doMCbyTP = cms.bool(False),
    tpSrc = cms.InputTag('mergedtruth','MergedTrackTruth'),
    doTrigEffCorrection = cms.bool(False),
    trigEffByMult = cms.vdouble(0.0,1.0),
    trigContaminationByMult = cms.vdouble(0.0,0.0),
    zeroMultFraction = cms.double(0.0),
    vtxWeightParameters = cms.vdouble( 4.49636e-02, 1.36629e-01, 5.30010e+00,
                                       2.50170e-02, 4.59123e-01, 9.64888e+00 ),
    doVtxReweighting = cms.bool(False)
)
