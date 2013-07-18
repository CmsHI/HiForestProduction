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
      #'/store/data/Commissioning10/MinimumBias/RECO/Apr20ReReco-v1/0164/D237203A-DC4C-DF11-BBF5-0018F3D095FC.root'))
      #'/store/user/edwenger/MinimumBias/MB-C10-A20RR-TRKANASKIM-v4/38a0275e6a2a62980cf9b1868edbd493/trkAnaSkimAOD_97_1.root'))
      # jun14rr: run 132440 - 137436
      #'/store/data/Run2010A/JetMETTau/RECO/Jun14thReReco_v2/0000/5C4D249B-4179-DF11-8CE3-0025B3E06396.root',
      #'/store/data/Run2010A/JetMETTau/RECO/Jun14thReReco_v2/0000/56D1F17A-3979-DF11-A340-003048D476DA.root',
      #'/store/data/Run2010A/JetMETTau/RECO/Jun14thReReco_v2/0000/4A5A42A3-4379-DF11-A264-003048636198.root'
      # jmt pr4: run >137436 (CMSSW_3_6_1_patch4) 
      #'/store/data/Run2010A/JetMETTau/RECO/v4/000/140/059/28ADC2FF-708E-DF11-B715-001D09F2527B.root',
      #'/store/data/Run2010A/JetMETTau/RECO/v4/000/140/059/24F181DC-8A8E-DF11-BCFE-003048F1183E.root',
      #'/store/data/Run2010A/JetMETTau/RECO/v4/000/140/059/22228061-A18E-DF11-9FAE-003048F1183E.root'
      # jm pr4: run > 141950
      '/store/data/Run2010A/JetMET/RECO/v4/000/144/114/F4E53E34-31B4-DF11-9A2C-0030487CD6D8.root',
      '/store/data/Run2010A/JetMET/RECO/v4/000/144/114/C8E66A3C-28B4-DF11-BA17-001D09F24763.root'
      # jm pr4 skim
      #'file:/home/sungho/sctch101/data/spectra/7TeV/test/JM-R10A-PR4-Jet50Uskim-v0-trkAnaSkimAOD_530_3_kOu.root'
      ))
    # =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR10_P_V7::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.23 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/TrkAnaSkim_data_7TeV_cfg.py,v $'),
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
enableHLTJet(process,"HLT_Jet50U")
updateEvtSelEff(process.trackAna_STD,"STD_NSD_TrkVtx")
updateEvtSelEff(process.looseTrackAna_STD,"STD_NSD_PixVtx")
enableDJetAna(process,"Data")
setMaxNumberVertices(process,1)
runWithsortedSumPtOrdVertices(process)
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
    #process.output_step
    )
