import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 2.76 TeV MC Sample =====================

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
    '/store/mc/Spring10/Hydjet_Quenched_MinBias_2760GeV/GEN-SIM-RECO/MC_36Y_V7A-v1/0005/E82E5989-BE70-DF11-9952-001EC9AA9FCC.root'))

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'MC_37Y_V4::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/HiTrkAnaSkim_data_2TeV_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('trkhistsMC.root')
                                   )

# =============== Import Sequences =====================
process.load("edwenger.Skims.HiAnalysis_cff")
process.load("edwenger.Skims.HiExtraReco_cff")

from edwenger.Skims.hicustomise_cfi import *
#process = enableSIM(process)    # activate isGEN in analyzers

# =============== Pat jet in HI ========================
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from Saved.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Spring10") # Spring10 is the new default

# =============== Final Paths =====================
process.prod        = cms.Path(process.hiextraReco + process.makeHeavyIonJets)
process.ana_step    = cms.Path(process.hiAnalysisSeq)

#process = disableLowPt(process) # disable low pt pixel
process = disableEffAna(process) # disable eff

from Saved.JulyExercise.customise_cfi import *
disableMC(process)

# =============== Output ================================
#process.load("edwenger.Skims.analysisSkimContent_cff")
#process.output = cms.OutputModule("PoolOutputModule",
#    process.analysisSkimContent,
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
#    dataset = cms.untracked.PSet(
#      dataTier = cms.untracked.string('AODSIM'),
#      filterName = cms.untracked.string('TrkAnaFilter')),
#    fileName = cms.untracked.string('trkAnaSkimAODSIM.root')
#    )

#process.output_step = cms.EndPath(process.output)



# =============== Schedule =====================

process.schedule = cms.Schedule(
    process.prod,
    process.ana_step,
    #process.output_step
    )
