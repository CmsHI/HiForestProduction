import FWCore.ParameterSet.Config as cms

process = cms.Process('CORRFLOWSKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
process.load('Appeltel.PixelTracksRun2010.HiTrackCandidates_cff')
process.load('Appeltel.PixelTracksRun2010.HICorrFlowSkimEventContent_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)
process.options = cms.untracked.PSet(
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'file:/scratch/appelte1/pixelMC/Hydjet_Bass_MinBias_2760GeV_1.root' #T3_US_Vanderbilt
#'file:/net/hidsk0001/d00/scratch/appelte1/Hydjet_Bass_MinBias_2760GeV_1.root' #T2_US_MIT
)
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.HICorrFlowEventContent.outputCommands,
    fileName = cms.untracked.string( 'file:association.root' )
)

process.output.outputCommands += [ 'keep *_tpRecoAssocLowPtPixelTracks_*_*' ]
process.output.outputCommands += [ 'keep *_tpRecoAssocMerged12Tracks_*_*' ]
process.output.outputCommands += [ 'keep *_mergedtruth_*_*' ]
process.output.outputCommands += [ 'keep recoTracks_hiMergedPtSplit12Tracks_*_*' ]
process.output.outputCommands += [ 'keep recoTrackExtras_hiMergedPtSplit12Tracks_*_*' ]
process.output.outputCommands += [ 'keep recoTrackExtras_hiLowPtPixelTracks_*_*' ]
process.output.outputCommands += [ 'keep TrackingRecHitsOwned_hiLowPtPixelTracks_*_*' ]
process.output.outputCommands += [ 'keep TrackingRecHitsOwned_hiMergedPtSplit12Tracks_*_*' ]
process.output.outputCommands += [ 'keep *_offlineBeamSpot_*_*' ]
process.output.outputCommands += [ 'keep *_siPixelClusters_*_*' ]
process.output.outputCommands += [ 'keep *_siStripClusters_*_*' ]

# Other statements
process.GlobalTag.globaltag = 'START39_V4HI::All'


process.hiMergedPtSplit12Tracks.copyExtras = True

#####################################################
#
# TrackAssociation
#
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")


process.tpRecoAssocMerged12Tracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocMerged12Tracks.label_tr = cms.InputTag("hiMergedPtSplit12Tracks")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

process.Timing = cms.Service("Timing")



# Path and EndPath definitions

process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
process.lowptpixel_step = cms.Path(process.lowPtPixelTrackReco)
process.merge_step = cms.Path(process.multipleMergedTracks)
process.candidates_step = cms.Path( process.hiAllTrackCandidates)
process.association_step = cms.Path( process.tpRecoAssocMerged12Tracks)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.reconstruction_step,
    process.lowptpixel_step,
    process.merge_step,
    process.candidates_step,
    process.association_step,
    process.endjob_step,
    process.out_step
)



