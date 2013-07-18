import FWCore.ParameterSet.Config as cms
process = cms.Process("TRKANA")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 7 TeV Collision Data =====================

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/data/Commissioning10/MinimumBias/RECO/Apr20ReReco-v1/0164/D237203A-DC4C-DF11-BBF5-0018F3D095FC.root'))
      '/store/user/edwenger/MinimumBias/MB-C10-PR9-TRKANASKIM-v3/ae98f896d123ace1f592d26e790fa90c/trkAnaSkimAOD_100_1.root'))
      #'/store/user/edwenger/MinimumBias/MB-C10-A20RR-TRKANASKIM-v2/86d28cd0599312fbc0b38fb077d9e1fc/trkAnaSkimAOD_10_1.root'))
      #'/store/user/edwenger/MinimumBias/MB-C10-PR9-MBskim-v0/41db2a25372977708462def4fd0d0a92/trkAnaSkimAOD_89_1.root'))
      #'/store/user/edwenger/MinimumBias/MB-C10-PR9-JMTskim-v0/e8f20a07d94a41822e93583a611264a6/trkAnaSkimAOD_97_1.root'))
      #'/store/user/frankma/MinimumBias/MB-C10-M6RR-MBHfLooseskim-v0/731e11dd8899c883730fd3f1e4a373e3/trkAnaSkimAOD_282_2_hxv.root'))
      # jm pr4: run > 141950
      #'/store/user/frankma/JetMET/JM-R10A-PR4-Jet50Uskim-v0/122e6996b8ad1b72a263641eb6b3d37b/trkAnaSkimAOD_9_1_eQU.root'))

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR10_P_V7::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/prod/condor/TrkAnaSkim_data_7TeV_cfg.py,v $'),
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
#switchJECSet( process, "Summer09_7TeV_ReReco332") # if not explicit then takes from Relesae default

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All']) # turn off MC matching for data

from edwenger.Skims.customise_cfi import *
#enableEitherHFEvtSel(process) # Replace HF coinc with the looser any HF hit evt selection
#enableHLTJet(process,"HLT_Jet50U")
updateEvtSelEff(process.trackAna_STD,"STD_NSD_TrkVtx")
updateEvtSelEff(process.looseTrackAna_STD,"STD_NSD_PixVtx")
enableDJetAna(process,"Data","FF")
#setMaxNumberVertices(process,1)
#runWithsortedSumPtOrdVertices(process)
enableAOD(process)

# =============== Final Paths =====================

process.eventFilter_step = cms.Path(process.eventFilter)
process.extraReco_step   = cms.Path(process.extraReco)
process.ana_step         = cms.Path(process.analysisSeq)

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
    process.ana_step#,
    #process.output_step
    )
