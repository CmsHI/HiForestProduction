import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKINGDEBUG')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
)

process.source = cms.Source("PoolSource",
   fileNames =
cms.untracked.vstring('/store/hidata/HIRun2010/HIAllPhysics/RECO/SDmaker_3SD_1CS_PDHIAllPhysicsZSv2_SD_JetHI-v1/0055/A2915765-B04F-E011-BB14-003048FEADDC.root'
)
)

process.GlobalTag.globaltag = 'GR_R_39X_V6B::All'

process.rechit_step = cms.Path(process.siPixelRecHits)
process.lowptpixel_step = cms.Path(process.hiConformalPixelTracks)

process.schedule = cms.Schedule(
   process.rechit_step,
   process.lowptpixel_step
)
