import FWCore.ParameterSet.Config as cms

import Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cfi

trkAna_quality = Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cfi.trkAna.clone(
        applyCuts = cms.bool(True)
)

trkAna_pixel = Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cfi.trkAna.clone(
   trackSrc = cms.InputTag("pixelTracks")
)

trkAna_noBS = Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cfi.trkAna.clone(
   vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


trkAna_HIN12017 = Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cfi.trkAna.clone(
   applyCuts = cms.bool(True),
   qualityString = cms.string('highPurity'),
   dxyErrMax = cms.double(3.0),
   dzErrMax = cms.double(3.0),
   ptErrMax = cms.double(0.1),
   vertexSrc = cms.InputTag("offlinePrimaryVertices"),
   ptBins = cms.vdouble(
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4,
        7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6,
        130.0
   ),
   etaBins = cms.vdouble( 
        -2.46, -2.36, -2.26, -2.16, -2.06, -1.96, -1.86, -1.76, -1.66, -1.56,
        -1.46, -1.36, -1.26, -1.16, -1.06, -0.96, -0.86, -0.76, -0.66, -0.56,
        -0.46, -0.36, -0.26, -0.16, -0.06, 0.04, 0.14, 0.24, 0.34, 0.44, 0.54, 
        0.64, 0.74, 0.84, 0.94, 1.04, 1.14, 1.24, 1.34, 1.44, 1.54, 1.64, 1.74,
        1.84, 1.94, 2.04, 2.14, 2.24, 2.34, 2.44
   ),
   occByLeadingJetEt = cms.bool(True),
   jetEtaMax = cms.double(2.5),
   occBins = cms.vdouble( 
        0., 20., 40., 60., 80., 100., 120., 140., 160., 
        180., 200., 250., 300., 500.
   ),
   vertexZMax = cms.double(15.),
   doTrigEffCorrection = cms.bool(True),
   trigEffByMult = cms.vdouble( 
        0.0, 0.0448119, 0.513363, 0.77141, 0.883644, 0.926159, 0.948583, 
        0.960119, 0.966817, 0.972921, 0.976898, 0.980945, 0.983142, 0.986, 0.988247, 
        0.989673, 0.991689, 0.993215, 0.993808, 0.994555, 0.996014, 0.996399, 0.996832, 
        0.996901, 0.997528, 0.997738, 0.998142, 0.998412, 0.99869, 0.998518
   ),
   trigContaminationByMult = cms.vdouble(
        0.0, 0.0150659, 0.0134671, 0.00927223, 0.00814763, 0.00572825, 0.00545534, 
        0.00500624, 0.00392851, 0.00290192, 0.0027658, 0.00256728, 0.00172565, 0.00151442, 0.00191661, 
        0.0013661, 0.000786178, 0.000747828, 0.000770812, 0.000412993, 0.000491937, 0.000288172, 0.000299312, 
        0.00021031, 0.000264633, 0.000192915, 0.000148894, 0.000151408, 0.000180004, 0.000156213
   ),
   zeroMultFraction = cms.double( 0.00142393)                        
)

trkAna_HIN12017_MC = trkAna_HIN12017.clone(
    doMC = cms.bool(True),
    genSrc = cms.InputTag("genParticles")
)

trkAna_HIN12017_TP = trkAna_HIN12017.clone(
    doMC = cms.bool(True),
    genSrc = cms.InputTag("genParticles"),
    doMCbyTP = cms.bool(True)
)

trkAna_HIN12017_Pbp = trkAna_HIN12017.clone(
   etaBins = cms.vdouble(
        -2.44, -2.34, -2.24, -2.14, -2.04, -1.94, -1.84, -1.74, -1.64, -1.54,
        -1.44, -1.34, -1.24, -1.14, -1.04, -0.94, -0.84, -0.74, -0.64, -0.54,
        -0.44, -0.34, -0.24, -0.14, -0.04, 
        0.06, 0.16, 0.26, 0.36, 0.46, 0.56,
        0.66, 0.76, 0.86, 0.96, 1.06, 1.16, 1.26, 1.36, 1.46, 1.56, 1.66, 1.76,
        1.86, 1.96, 2.06, 2.16, 2.26, 2.36, 2.46
   )
   # Note: May need to use different trigger efficiency
)
