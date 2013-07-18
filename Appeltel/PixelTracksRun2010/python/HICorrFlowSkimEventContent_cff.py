import FWCore.ParameterSet.Config as cms


HICorrFlowEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( 'drop *',
        'keep recoTracks_hiGoodTracks_*_*',
        'keep recoTracks_hiQualityMergedTracks_*_*',
        'keep recoTracks_hiSelectedTracks_*_*', 
        'keep recoTracks_hiGoodPixelTracks_*_*',
        'keep recoTracks_hiConformalPixelTracks_*_*',
        'keep recoTracks_hiGoodMergedTracks_*_*',
        'keep recoTracks_hiLowPtPixelTracks_*_*',
        'keep recoTracks_hiLowPtWideCutPixelTracks_*_*',
        'keep recoTracks_hiMergedNoPtSplitTracks_*_*',
        'keep *_hiSelectedVertex_*_*',
        'keep *_hiCentrality_*_*',
        'keep *_hiEvtPlane_*_*',
        'keep recoCaloJets_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*',
        'keep CaloTowersSorted_towerMaker_*_*',
        'keep recoCaloClusters_*_*_*',
        'keep *_hiGenParticles_*_*',
        'keep *_generator_*_*',
        'keep recoGenJets_*_*_*'
    )
)
 
