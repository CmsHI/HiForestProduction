import FWCore.ParameterSet.Config as cms

selectTracks = cms.EDProducer("HackedAnalyticalTrackSelector",

    src = cms.InputTag("generalTracks"),
    keepAllTracks = cms.bool(False), ## if set to true tracks failing this filter are kept in the output
    beamspot = cms.InputTag("offlineBeamSpot"),

    vertices = cms.InputTag("offlinePrimaryVertices"),
    vtxNumber = cms.int32(0),
    vtxTracks = cms.uint32(0), ## at least 3 tracks
    vtxChi2Prob = cms.double(0.0), ## at least 1% chi2nprobability (if it has a chi2)

    copyTrajectories = cms.untracked.bool(False),
    copyExtras = cms.untracked.bool(False), ## set to false on AOD
    qualityBit = cms.string(''), ## set to '' or comment out if you don't want to set the bit

    # parameters for adapted optimal cuts on chi2 and primary vertex compatibility
    chi2n_par = cms.double(99999.), # no longer in use
    chi2n_par_pt = cms.vdouble(0.522, 0.432, 0.024, -2.05, 2.3), 
    chi2n_par_pt2 = cms.vdouble(0.522, 0.45, 0.0294, -38.5, 1.775),
    res_par = cms.vdouble(99999., 99999.),
    d0_par1 = cms.vdouble(99999., 0.0),
    dz_par1 = cms.vdouble(99999., 0.0),
    d0_par2 = cms.vdouble(3, 0.0),
    dz_par2 = cms.vdouble(3, 0.0),
    # Boolean indicating if adapted primary vertex compatibility cuts are to be applied.
    applyAdaptedPVCuts = cms.bool(True),

    # Impact parameter absolute cuts.
    max_d0 = cms.double(0.2),
    max_z0 = cms.double(0.2),

    # Cuts on numbers of layers with hits/3D hits/lost hits. 
    minNumberLayers = cms.uint32(0),
    minNumber3DLayers = cms.uint32(0),
    maxNumberLostLayers = cms.uint32(99999)

)
