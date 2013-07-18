import FWCore.ParameterSet.Config as cms

process = cms.Process('HLTAna')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('HLTrigger.Configuration.HLT_HIon_data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Analyzers
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.UseTFileService			= cms.untracked.bool(True)
process.hltbitorg = process.hltbitanalysis.clone(
  l1GtReadoutRecord = 'hltGtDigis::'+process.name_() # if running from reco use: "gtDigis"
  ) 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Other statements
process.GlobalTag.globaltag = 'GR_R_44_V2::All'
fname='openhlt.root'

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
      #'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/150/590/02D89852-5DEC-DF11-9E11-001D09F282F5.root'
			'/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0001/D8E2C34A-D2AE-E011-9E1E-001A6478AC14.root',
      #'/store/data/Run2011A/MinimumBias/RAW/v1/000/173/692/3C024415-BECC-E011-96B3-001D09F244BB.root'
      )
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('l1EmulatorFromRaw nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('l1EmulatorFromRaw_RAW2DIGI_L1.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI-RECO')
    )
)
process.FEVTDEBUGHLToutput.outputCommands = ['keep *']


# Additional output definition
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(fname)                                  
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess*process.hltbitorg)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])
