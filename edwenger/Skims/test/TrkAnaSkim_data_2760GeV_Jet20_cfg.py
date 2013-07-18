import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 2760 GeV Collision Data =====================

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/home/sungho/sctch101/data/spectra/2760GeV/AA4A07CD-775A-E011-990E-001617C3B6E2.root'
      ))

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V17::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/TrkAnaSkim_data_2760GeV_Jet20_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('trkhists.root')
                                   )
#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#    ignoreTotal=cms.untracked.int32(0),
#    oncePerEventMode = cms.untracked.bool(False)
#    )


# =============== Import Sequences =====================

process.load("edwenger.Skims.EventFilter_cff")
process.load("edwenger.Skims.ExtraReco_cff")
process.load("edwenger.Skims.Analysis_cff")

from edwenger.Skims.customise_cfi import *
#enableEitherHFEvtSel(process) # Replace HF coinc with the looser any HF hit evt selection
enableHLTJet(process,"HLT_Jet20_v1")
updateEvtSelEff(process.trackAna_STD,"STD_NSD_TrkVtx")
updateEvtSelEff(process.looseTrackAna_STD,"STD_NSD_PixVtx")
#enableDJetAna(process,"Data")
removeDijetAna(process) # remove dijetAna_seq 
setMaxNumberVertices(process,1)
runWithsortedSumPtOrdVertices(process)
run2760GeVmode(process)
#process = enableAOD(process)

# =============== Final Paths =====================

process.eventFilter_step = cms.Path(process.eventFilter)
process.extraReco_step   = cms.Path(process.eventFilter * process.extraReco)
process.ana_step         = cms.Path(process.eventFilter * process.analysisSeq)

# =============== Output ================================

#process.load("edwenger.Skims.analysisSkimContent_cff")
#process.output = cms.OutputModule("PoolOutputModule",
#    process.analysisSkimContent,
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
#    dataset = cms.untracked.PSet(
#      dataTier = cms.untracked.string('AOD'),
#      filterName = cms.untracked.string('TrkAnaFilter')),
#    fileName = cms.untracked.string('trkAnaSkimAOD.root')
#    )

#process.output_step = cms.EndPath(process.output)

# =============== Schedule =====================

process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.extraReco_step,
    process.ana_step,
    #process.output_step
    )
