# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: l1EmulatorFromRaw -s RAW2DIGI,L1,HLT:HIon -n 100 --conditions auto:hltonline --datatier DIGI-RECO --eventcontent FEVTDEBUGHLT --data --filein /store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/150/590/02D89852-5DEC-DF11-9E11-001D09F282F5.root --customise=L1Trigger/Configuration/customise_l1EmulatorFromRaw --custom_conditions=L1GtTriggerMenu_L1Menu_Collisions2011_v5,L1GtTriggerMenuRcd,frontier://FrontierProd/CMS_COND_31X_L1T --processName=L1EmulRawHLT --no_exec
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
#process.load('HLTrigger.Configuration.HLT_HIon_Jet_data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
    #'file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/0EDFD1BD-C9AE-E011-963E-003048CB82AC.root',
    #'file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/B03AE47D-C2AE-E011-A408-003048CB96F8.root',
    #'file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/F8A193D3-C2AE-E011-8D37-003048CF667C.root',
    #'file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/AA944F78-C2AE-E011-8D7D-003048CB87A6.root',
    #'file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/50D1C4BF-C3AE-E011-931A-003048CF6578.root'
    #'file:l1EmulatorFromRawRepackRaw_RAW2DIGI_L1_DIGI2RAW_officialL1Menu_mctag_100.root'
    'file:/net/hisrv0001/home/davidlw/scratch1/HLTStudies/CMSSW_4_4_1_L1Repack/src/l1EmulatorFromRawRepackRaw_RAW2DIGI_L1_DIGI2RAW.root'
    )
)
newSource = True
doL1Emul = False
doL1Xml = False # the new HIon cff contains the customization to use xml

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.4 $'),
    annotation = cms.untracked.string('l1EmulatorFromRaw nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('HLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI-RECO')
    )
)

# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'GR_R_44_V5::All'
#process.GlobalTag.globaltag = 'GR_R_44_V5::All'
process.GlobalTag.globaltag = 'START44_V6::All'
process.GlobalTag.toGet = cms.VPSet()
process.GlobalTag.toGet.append(cms.PSet(tag=cms.string("L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2011_v0_mc"),record=cms.string("L1GtTriggerMenuRcd"),connect=cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T"),))

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.endjob_step = cms.EndPath(process.endOfProcess)
process.endjob_step = cms.EndPath()

# Schedule definition
process.schedule = cms.Schedule()
if doL1Emul:
  process.schedule.extend([process.raw2digi_step,process.L1simulation_step])
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])

# customisation of the process.

if doL1Emul:
  # Automatic addition of the customisation function from L1Trigger.Configuration.customise_l1EmulatorFromRaw
  from L1Trigger.Configuration.customise_l1EmulatorFromRaw import customise 
  #call to customisation function customise imported from L1Trigger.Configuration.customise_l1EmulatorFromRaw
  process = customise(process)
  import HLTrigger.Configuration.customizeHLTforL1Emulator
  # customize the HLT to use the emulated results
  process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
  process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimGmtGctGtDigis( process )

# set l1 menu from xml
if doL1Xml:
  import L1Trigger.Configuration.L1Trigger_custom
  process = L1Trigger.Configuration.L1Trigger_custom.customiseL1Menu( process )

# report
process.HLTAnalyzerEndpath.remove(process.hltL1GtTrigReport)
# L1 reports
from L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi import l1GtTrigReport
process.hltL1GtTrigReport = l1GtTrigReport.clone(L1GtRecordInputTag = cms.InputTag( 'hltGtDigis' ))
if doL1Emul:
  process.hltL1GtTrigReport.L1GtRecordInputTag = cms.InputTag( 'simGtDigis' )
process.MessageLogger.categories.append('L1GtTrigReport')
process.endjob_step*=process.hltL1GtTrigReport

if (newSource):
    from FWCore.ParameterSet import Mixins
    for module in process.__dict__.itervalues():
        if isinstance(module, Mixins._Parameterizable):
            for parameter in module.__dict__.itervalues():
                if isinstance(parameter, cms.InputTag):
                    if parameter.moduleLabel == 'source':
                        parameter.moduleLabel = 'sourcenew'
                        
# End of customisation functions
