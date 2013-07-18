
import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("SKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/EventContent/EventContent_cff')


# =============== input file setting =====================

# =============== 2.36 TeV MC Sample =====================

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/cms/store/caf/user/frankma/HR10Exp3/r150308HFSkim/skim_RECO_9_1_2w4.root")
                            )


# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# =============== Trigger selection ====================
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.HLTPaths = ["HLT_HIMinBiasHF"]

# =============== Output ================================
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
    fileName = cms.untracked.string("skim_output.root")
)

process.eventFilter_step = cms.Path(process.hltHighLevel)
process.output_step = cms.EndPath(process.output)




