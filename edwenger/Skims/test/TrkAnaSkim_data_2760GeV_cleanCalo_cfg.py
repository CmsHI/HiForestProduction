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
      #'file:/home/sungho/sctch101/data/spectra/2760GeV/trkAnaSkimAOD.root'
      ))

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V17::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/TrkAnaSkim_data_2760GeV_cleanCalo_cfg.py,v $'),
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
#enableHLTJet(process,"HLT_Jet20_v1")
updateEvtSelEff(process.trackAna_STD,"STD_NSD_TrkVtx")
updateEvtSelEff(process.looseTrackAna_STD,"STD_NSD_PixVtx")
#enableDJetAna(process,"Data")
#setMaxNumberVertices(process,1)
#runWithsortedSumPtOrdVertices(process)
run2760GeVmode(process)
#process = enableAOD(process)
enableOpenHlt(process,process.analysisSeq)

# === jet ana ===
process.load('MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_cff')
process.ak5JetAnalyzer = process.inclusiveJetAnalyzer.clone(
	jetTag = 'cleanedPatJets', vtxTag=cms.untracked.InputTag('sortedGoodVertices'),
	isMC=False,useCentrality=False)
process.analysisSeq*=process.ak5JetAnalyzer

process.load('Saved.CaloTools.calotowerscleaner_cfi')
process.cleanedTowers.excludeEtaMin = cms.double(-1.5)
process.cleanedTowers.excludeEtaMax = cms.double(0)
process.cleanedTowers.excludePhiMin = cms.double(0.2)
process.cleanedTowers.excludePhiMax = cms.double(0.8)

process.load('RecoJets.Configuration.RecoJets_cff')
process.ak5CaloJets = process.ak5CaloJets.clone(src = cms.InputTag("cleanedTowers"))
#process.patJets.addJetID = cms.bool(False) # uncomment if running on skim
process.patJets.embedPFCandidates = cms.bool(False)
process.patJets.embedGenJetMatch = cms.bool(False)
process.patJets.addAssociatedTracks = cms.bool(False)
process.patJets.trackAssociationSource = ""

process.load('RecoJets.Configuration.JetIDProducers_cff')

# =============== Final Paths =====================

process.eventFilter_step = cms.Path(process.eventFilter)
process.cleanReco_step = cms.Path(process.cleanedTowers * process.ak5CaloJets * process.ak5JetID)
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
process.output.outputCommands = ['keep *']

process.output_step = cms.EndPath(process.output)

# =============== Schedule =====================

process.schedule = cms.Schedule(
    process.eventFilter_step,
		process.cleanReco_step,
    process.extraReco_step,
    process.ana_step,
    process.output_step
    )
