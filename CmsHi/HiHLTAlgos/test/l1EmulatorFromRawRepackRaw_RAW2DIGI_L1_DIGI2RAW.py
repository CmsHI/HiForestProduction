# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: l1EmulatorFromRawRepackRaw -s RAW2DIGI,L1,DIGI2RAW -n 10 --conditions auto:hltonline --datatier DIGI-RAW --eventcontent FEVTDEBUGHLT --data --filein file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/0EDFD1BD-C9AE-E011-963E-003048CB82AC.root --customise=L1Trigger/Configuration/customise_l1EmulatorFromRaw --custom_conditions=L1GtTriggerMenu_L1Menu_Collisions2011_v5_mc,L1GtTriggerMenuRcd,frontier://FrontierProd/CMS_COND_31X_L1T --processName=L1EmulRawRepackRaw --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('L1EmulRawRepackRaw')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/0EDFD1BD-C9AE-E011-963E-003048CB82AC.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('l1EmulatorFromRawRepackRaw nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    #outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = process.RAWEventContent.outputCommands,
    fileName = cms.untracked.string('l1EmulatorFromRawRepackRaw_RAW2DIGI_L1_DIGI2RAW.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI-RAW')
    )
)
process.FEVTDEBUGHLToutput.outputCommands.extend(['keep FEDRawDataCollection_*_*_*'])
process.FEVTDEBUGHLToutput.outputCommands.extend(['drop FEDRawDataCollection_source_*_*'])

# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'GR_R_44_V5::All'
process.GlobalTag.globaltag = 'START44_V6::All'
process.GlobalTag.toGet = cms.VPSet()
process.GlobalTag.toGet.append(cms.PSet(tag=cms.string("L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2011_v0_mc"),record=cms.string("L1GtTriggerMenuRcd"),connect=cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T"),))

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

# customisation of the process
process.gctDigiToRaw.rctInputLabel = 'simRctDigis'
process.gctDigiToRaw.gctInputLabel = 'simGctDigis'
process.l1GtPack.DaqGtInputTag = 'simGtDigis'
process.l1GtPack.MuGmtInputTag = 'simGmtDigis'
process.l1GtEvmPack.EvmGtInputTag = 'simGtDigis'
from EventFilter.RawDataCollector.rawDataCollectorByLabel_cfi import rawDataCollector
process.sourcenew = rawDataCollector.clone(
    verbose = cms.untracked.int32(0),
    RawCollectionList = cms.VInputTag( cms.InputTag('gctDigiToRaw'),
                                       cms.InputTag('l1GtPack'),
                                       cms.InputTag('l1GtEvmPack'),
                                       cms.InputTag('source')
                                       #cms.InputTag('rawDataCollector')
                                       )
    )
process.DigiToRaw = cms.Sequence(#csctfpacker*dttfpacker*
  process.gctDigiToRaw*process.l1GtPack*process.l1GtEvmPack*
  #*siPixelRawData*SiStripDigiToRaw*ecalPacker*esDigiToRaw*hcalRawData*cscpacker*dtpacker*rpcpacker*castorRawData*
  process.sourcenew
  )

# Automatic addition of the customisation function from L1Trigger.Configuration.customise_l1EmulatorFromRaw
from L1Trigger.Configuration.customise_l1EmulatorFromRaw import customise 

#call to customisation function customise imported from L1Trigger.Configuration.customise_l1EmulatorFromRaw
process = customise(process)

# set l1 menu from xml
#import L1Trigger.Configuration.L1Trigger_custom
#process = L1Trigger.Configuration.L1Trigger_custom.customiseL1Menu( process )

# End of customisation functions
