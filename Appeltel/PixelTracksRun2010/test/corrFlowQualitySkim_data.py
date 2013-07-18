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

process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")

process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
process.load('Appeltel.PixelTracksRun2010.HiTrackCandidates_cff')
process.load('Appeltel.PixelTracksRun2010.HICorrFlowSkimEventContent_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet(
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2A14D410-44EE-DF11-949B-003048F118DE.root'
  )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.HICorrFlowEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('filter_step')
                ),
    fileName = cms.untracked.string( 'file:skim_quality_test.root')
)

#Trigger Selection

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMinBiasHF = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMinBiasHF.HLTPaths = ["HLT_HIMinBiasHfOrBSC_Core"]

# Other statements
process.GlobalTag.globaltag = 'GR10_P_V12::All'


# Path and EndPath definitions

process.reco_step = cms.Path( process.hltMinBiasHF * 
                   process.collisionEventSelection *
                   process.hiGoodTracksSelection *
                   process.goodPixelTrackReco *
                   process.hiQualityMergedTracks 
)

process.filter_step = cms.Path( process.hltMinBiasHF * process.collisionEventSelection )
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(
    process.filter_step,
    process.reco_step,
    process.endjob_step,
    process.out_step
)

