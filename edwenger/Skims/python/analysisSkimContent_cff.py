import FWCore.ParameterSet.Config as cms

analysisSkimContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *',
                                           
      # event
      'keep *_offlinePrimaryVertices_*_*',
      'keep *_pixelVertices_*_*',
      'keep *_pixel3Vertices_*_*',                                     
      'keep *_offlineBeamSpot_*_*',
      'keep *_TriggerResults_*_HLT',
      'keep *_TriggerResults_*_REDIGI',
      'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*RECO',

      # mc gen info
      'keep *_*GenJet*_*_*',
      'keep reco*_genParticles_*_*',
      'keep GenEventInfoProduct_*_*_*',
                                           
      # jet
      'keep *_towerMaker_*_*',
      'keep patJets_selectedPatJets*_*_*',
                                           
      # tracks
      'keep recoTracks_generalTracks_*_*RECO',  
      'keep recoTracks_refitTracks_*_*',
      'keep recoRecoChargedCandidates_allTracks_*_*',
      'keep recoTracks_selectFake*_*_*',
      'keep recoTracks_selectReal*_*_*',
                                           
      # sim track matching
      'keep *_trackingParticleRecoTrackAsssociation_*_*',
      'keep TrackingParticles_mergedtruth_MergedTrackTruth_*'
                                           
      )
    )
