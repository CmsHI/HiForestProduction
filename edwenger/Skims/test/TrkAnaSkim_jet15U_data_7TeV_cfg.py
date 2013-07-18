import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 7 TeV Collision Data =====================

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Commissioning10/MinimumBias/RECO/May6thPDSkim2_SD_JetMETTau-v1/0137/5869167B-D25D-DF11-9BE5-00248C0BE013.root'),
    '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/035/96AC9E3D-B865-DF11-90B2-0030487CD180.root',
    '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/035/48B7C943-D965-DF11-A486-0030487CD700.root'),
    lumisToProcess = cms.untracked.VLuminosityBlockRange(
    '136035:1-136035:max'))

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.GlobalTag.globaltag = 'GR_R_35X_V7A::All' #357
process.GlobalTag.globaltag ='GR_R_35X_V8B::All' #358p3

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/TrkAnaSkim_jet15U_data_7TeV_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('trkhists.root')
                                   )
process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
    ignoreTotal=cms.untracked.int32(0),
    oncePerEventMode = cms.untracked.bool(False)
    )


# =============== Import Sequences =====================

process.load("edwenger.Skims.EventFilter_cff")
process.load("edwenger.Skims.ExtraReco_cff")
process.load("edwenger.Skims.Analysis_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Summer09_7TeV_ReReco332") # get the 7 TeV jet corrections

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All']) # turn off MC matching for data

from edwenger.Skims.customise_cfi import *
#process = enableAOD(process)
process = enableHLTJet(process)

# =============== Final Paths =====================

process.eventFilter_step = cms.Path(process.eventFilter)
process.extraReco_step   = cms.Path(process.eventFilter * process.extraReco)
process.ana_step         = cms.Path(process.eventFilter * process.analysisSeq)

# =============== Output ================================

process.load("edwenger.Skims.analysisSkimContent_cff")
process.output = cms.OutputModule("PoolOutputModule",
    process.analysisSkimContent,
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD'),
      filterName = cms.untracked.string('TrkAnaFilter')),
    fileName = cms.untracked.string('trkAnaSkimAOD.root')
    )

process.output_step = cms.EndPath(process.output)

# =============== Schedule =====================

process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.extraReco_step,
    process.ana_step,
    process.output_step
    )
