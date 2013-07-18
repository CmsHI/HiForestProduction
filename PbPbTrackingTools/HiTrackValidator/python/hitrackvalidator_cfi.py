import FWCore.ParameterSet.Config as cms

hitrkvalidator = cms.EDAnalyzer('HiTrackValidator',
      bslabel=cms.untracked.InputTag("offlineBeamSpot"),
      vtxlabel=cms.untracked.InputTag("hiSelectedVertex"),
      trklabel=cms.untracked.InputTag("hiSelectedTracks"),
      jetlabel=cms.untracked.InputTag("icPu5patJets"),
      simtrklabel=cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
      towerlabel=cms.untracked.InputTag("towerMaker"),
      associatorMap=cms.untracked.InputTag('trackingParticleRecoTrackAsssociation'),
      etaMax=cms.untracked.double(1.0),
      jetEtMin=cms.untracked.double(-1.0), # cut on leading jet ET of event (if -1.0, no cut)
      hasSimInfo=cms.untracked.bool(False),
      selectFake=cms.untracked.bool(False),
      hasCaloMat=cms.untracked.bool(True), # check track-calo matching
      towerPtMin = cms.untracked.double(5.0), # min tower pt
      useQaulityStr=cms.untracked.bool(False),
      qualityString=cms.untracked.string("highPurity"),
      fiducialCut = cms.untracked.bool(False), # fiducial cut for FPIX
      # root syntax, pt dependent calo-compatibility cut
      funcCaloComp = cms.string("0.2"),
      # root syntax, pt dependent deltaR matching cut
      funcDeltaRTowerMatch = cms.string("0.087/(1.0+0.1*exp(-0.28*(x-20.)))"),
      neededCentBins = cms.untracked.vint32(0,1,3,11,19,39)
)
