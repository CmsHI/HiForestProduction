# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: l1EmulatorFromRaw -s RAW2DIGI,L1,HLT:HIon -n 100 --conditions auto:hltonline --datatier DIGI-RECO --eventcontent FEVTDEBUGHLT --data --filein /store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/150/590/02D89852-5DEC-DF11-9E11-001D09F282F5.root --customise=L1Trigger/Configuration/customise_l1EmulatorFromRaw --custom_conditions=L1GtTriggerMenu_L1Menu_Collisions2011_v5,L1GtTriggerMenuRcd,frontier://FrontierProd/CMS_COND_31X_L1T --processName=L1EmulRawHLT --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RawSkim')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
    '/store/user/davidlw/HICorePhysics/l1EmulatorFromRawRepackRaw_441_v2/ddf8720368adbbb267958624404a73f3/l1EmulatorFromRawRepackRaw_RAW2DIGI_L1_DIGI2RAW_80_1_7ST.root'
    )
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.output = cms.OutputModule("PoolOutputModule",
           fileName = cms.untracked.string('l1EmulatorFromRawRepackRaw_RawSkim.root'),
           outputCommands = cms.untracked.vstring(
   'drop *_*_*_*',
   'keep *_sourcenew_*_*',
   'keep L1GlobalTriggerObjectMapRecord_*_*_*',
   'keep *_TriggerResults_*_*',
   'keep *_hltTriggerSummaryAOD_*_*'
   )
)

process.out_step = cms.EndPath(process.output)
