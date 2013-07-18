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
      # may6rr
      #'/store/data/Commissioning10/MinimumBias/RAW-RECO/May6thPDSkim_GOODCOLL-v1/0003/60245B69-BF5C-DF11-A0D5-0018F3D09634.root',
      #'/store/data/Commissioning10/MinimumBias/RAW-RECO/May6thPDSkim_GOODCOLL-v1/0003/5ECC499A-815C-DF11-9B24-0018F3D09688.root'
      # jun14rr, Run 134725
      #'/store/data/Commissioning10/MinimumBias/RECO/SD_JetMETTau-Jun14thSkim_v1/0149/826CD0DA-7D84-DF11-BA7D-0018F3D09648.root',
      #'/store/data/Commissioning10/MinimumBias/RECO/SD_JetMETTau-Jun14thSkim_v1/0129/DEB2CC2F-4380-DF11-A8DE-003048678F02.root'
      # jmt pr4: run >137436 (CMSSW_3_6_1_patch4) 
      '/store/data/Run2010A/JetMETTau/RECO/v4/000/140/059/28ADC2FF-708E-DF11-B715-001D09F2527B.root'
      ),
    )

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_R_36X_V12A::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.10 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/TrkAnaSkim_data_900GeV_cfg.py,v $'),
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
#switchJECSet(process, "Summer09_7TeV_ReReco332") # get the 900 GeV jet corrections, for the moment use default 7TeV JEC

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All']) # turn off MC matching for data

from edwenger.Skims.customise_cfi import *
updateEvtSelEff(process.trackAna_STD,"STD_NSD_900GeV_TrkVtx")
updateEvtSelEff(process.looseTrackAna_STD,"STD_NSD_900GeV_PixVtx")
#enableHLTJet(process)
enableDJetAna(process,"Data","LIGHT") # anaModes: "MC","Data", outLevels: "LIGHT","FF","FULL"
run900GeVmode(process) # set trigger names and appropriate binning for 900 GeV
#process = enableAOD(process)

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
