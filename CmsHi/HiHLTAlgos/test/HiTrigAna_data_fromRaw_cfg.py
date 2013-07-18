# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: l1EmulatorFromRaw -s RAW2DIGI,L1,HLT:HIon -n 100 --conditions auto:hltonline --datatier DIGI-RECO --eventcontent FEVTDEBUGHLT --data --filein /store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/150/590/02D89852-5DEC-DF11-9E11-001D09F282F5.root --customise=L1Trigger/Configuration/customise_l1EmulatorFromRaw --custom_conditions=L1GtTriggerMenu_L1Menu_Collisions2011_v5,L1GtTriggerMenuRcd,frontier://FrontierProd/CMS_COND_31X_L1T --processName=L1EmulRawHLT --no_exec
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process('HiTrigAna')

# setup 'analysis'  options
#options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
#options.outputFile = 'output.root'
#options.inputFiles = 'rfio:/castor/cern.ch/cms/store/data/ComissioningHI/HIExpressPhysics/RAW/v1/000/180/892/F6B10997-2B08-E111-A133-BCAEC532970F.root'
#options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
#options.parseArguments()

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
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/data/ComissioningHI/HIHighPt/RAW/v1/000/180/891/9AC450E9-3508-E111-8E63-003048CF94A6.root')
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Other statements
process.GlobalTag.globaltag = 'GR_E_V22::All'

# Analyzers
# openhlt
process.load('HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi')
process.hltbitanalysis.UseTFileService			= cms.untracked.bool(True)

process.hltana = process.hltbitanalysis.clone(
   l1GtReadoutRecord		= cms.InputTag("hltGtDigis"),
   l1GctHFBitCounts     = cms.InputTag("hltGctDigis"),
   l1GctHFRingSums      = cms.InputTag("hltGctDigis"),
   l1extramu            = cms.string('hltL1extraParticles'),
   l1extramc            = cms.string('hltL1extraParticles'),
   l1GtObjectMapRecord  = cms.InputTag("hltL1GtObjectMap"),
   hltresults           = cms.InputTag("TriggerResults","","HLT"),
   HLTProcessName       = cms.string("HLT")
  )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("openhlt_data.root")
)

# Path and EndPath definitions
process.recoextra_step = cms.Path(process.HLTL1UnpackerSequence)
process.ana_step = cms.Path(process.hltana)

# Schedule definition
process.schedule = cms.Schedule(process.recoextra_step,process.ana_step)

repacked=True
if (repacked):
    from FWCore.ParameterSet import Mixins
    for module in process.__dict__.itervalues():
        if isinstance(module, Mixins._Parameterizable):
            for parameter in module.__dict__.itervalues():
                if isinstance(parameter, cms.InputTag):
                    if parameter.moduleLabel == 'source':
                        parameter.moduleLabel = 'rawDataRepacker'
