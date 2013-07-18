# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: l1EmulatorFromRaw -s RAW2DIGI,L1,HLT:HIon -n 100 --conditions auto:hltonline --datatier DIGI-RECO --eventcontent FEVTDEBUGHLT --data --filein /store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/150/590/02D89852-5DEC-DF11-9E11-001D09F282F5.root --customise=L1Trigger/Configuration/customise_l1EmulatorFromRaw --custom_conditions=L1GtTriggerMenu_L1Menu_Collisions2011_v5,L1GtTriggerMenuRcd,frontier://FrontierProd/CMS_COND_31X_L1T --processName=L1EmulRawHLT --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('L1EmulRawHLT')

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/hidata/HIRun2010/HICorePhysics/RAW-RECO/v3/000/151/019/745B4D47-F3EE-DF11-92F8-0025B3E06424.root')
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    annotation = cms.untracked.string('l1EmulatorFromRaw nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('l1EmulatorFromRaw_RAW2DIGI_L1_HLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI-RECO')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'GR_R_44_V5::All'
#process.GlobalTag.toGet = cms.VPSet()
#process.GlobalTag.toGet.append(cms.PSet(tag=cms.string("L1GtTriggerMenu_L1Menu_Collisions2011_v5_mc"),record=cms.string("L1GtTriggerMenuRcd"),connect=cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T"),))

# Analyzers
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.UseTFileService			= cms.untracked.bool(True)

process.hltbitnew = process.hltbitanalysis.clone(
  l1GtObjectMapRecord	= cms.InputTag("hltL1GtObjectMap","",process.name_()),
  l1GtReadoutRecord		= cms.InputTag("simGtDigis","",process.name_()),
  hltresults           = cms.InputTag("TriggerResults","",process.name_()),
  HLTProcessName       = cms.string(process.name_())
  )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("openhlt_data_l1emul.root")
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess*process.hltbitnew)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1simulation_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customise_l1EmulatorFromRaw
from L1Trigger.Configuration.customise_l1EmulatorFromRaw import customise 

#call to customisation function customise imported from L1Trigger.Configuration.customise_l1EmulatorFromRaw
process = customise(process)

# set l1 menu from xml
import L1Trigger.Configuration.L1Trigger_custom
process = L1Trigger.Configuration.L1Trigger_custom.customiseL1Menu( process )

# customize the HLT to use the emulated results
import HLTrigger.Configuration.customizeHLTforL1Emulator
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimGmtGctGtDigis( process )

# End of customisation functions
