import FWCore.ParameterSet.Config as cms

from Appeltel.PixelTracksRun2010.hiTrackListMerger_cfi import *

hiMergedNoPtSplitTracks = hiTrackListMerger.clone(
    TrackProducer1 = 'hiSelectedTracks',
    TrackProducer2 = 'hiLowPtPixelTracks',
    promoteTrackQuality = False,
    allowFirstHitShare = False
)

hiMergedPtSplit09Tracks = hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 0.9,
    trackCollection1MinPt = 0.9
)

hiMergedPtSplit12Tracks = hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 1.2,
    trackCollection1MinPt = 1.2
)

hiMergedPtSplit15Tracks = hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 1.5,
    trackCollection1MinPt = 1.5
)

hiMergedPtSplit18Tracks = hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 1.8,
    trackCollection1MinPt = 1.8
)

hiQualityMergedTracks = hiTrackListMerger.clone(
    TrackProducer1 = 'hiGoodTracks',
    TrackProducer2 = 'hiGoodPixelTracks',
    promoteTrackQuality = False,
    allowFirstHitShare = False
)

hiGoodMergedTracks = hiTrackListMerger.clone(
    TrackProducer1 = 'hiGoodTracks',
    TrackProducer2 = 'hiConformalPixelTracks',
    promoteTrackQuality = False,
    trackCollection2MaxPt = 1.8,
    trackCollection1MinPt = 1.5,
    allowFirstHitShare = False
)



multipleMergedTracks = cms.Sequence ( hiMergedNoPtSplitTracks
                                      * hiMergedPtSplit09Tracks
                                      * hiMergedPtSplit12Tracks
                                      * hiMergedPtSplit15Tracks
                                      * hiMergedPtSplit18Tracks
                                      * hiQualityMergedTracks
)
