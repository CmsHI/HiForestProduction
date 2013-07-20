import FWCore.ParameterSet.Config as cms

# Track efficiency analyzer
hitrkEffAnalyzer = cms.EDAnalyzer('HiTrkEffAnalyzer',
    tracks = cms.untracked.InputTag('hiSelectedTracks'),
    jets = cms.untracked.InputTag('icPu5patJets'),
    label_tp_effic = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),
    label_tp_fake = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),                 
    associatormap = cms.untracked.InputTag('trackingParticleRecoTrackAsssociation'),
    vertices = cms.untracked.InputTag('hiSelectedVertex',''),
    beamspot = cms.untracked.InputTag('offlineBeamSpot'),                         
    fillHistograms = cms.bool(True),
    fillNtuples = cms.bool(False),
    ptBinScheme = cms.int32(1), # 1 = original, 2 = constant d(1/pt) binning, 3 = coarse, other then this, constant
    lowPtMode = cms.bool(False),
    doAssociation = cms.untracked.bool(True),
    hasSimInfo = cms.untracked.bool(False),
    pixelMultMode =  cms.untracked.bool(False),
    neededCentBins = cms.untracked.vint32(0,1,3,11,19,39), # 0~5% (0,1), 5~10% (2,3), 10~30% (4,11) etc..
    useJetEtMode = cms.int32(2), # 0 = no jet, 1 = leading jet, 2 = closest jets
    trkPtMin = cms.double(-1),
    trkAcceptedJet = cms.untracked.bool(True), # limit jets to |eta|<2.0 
    useSubLeadingJet = cms.untracked.bool(False),
    jetTrkOnly = cms.untracked.bool(False),
    fiducialCut = cms.untracked.bool(True), # fiducial cut for FPIX  
    useQaulityStr=cms.untracked.bool(False),
    qualityString=cms.untracked.string("highPurity"),
    usePxlPair = cms.untracked.bool(True) # include pixel pair in acceptance definition!
)
