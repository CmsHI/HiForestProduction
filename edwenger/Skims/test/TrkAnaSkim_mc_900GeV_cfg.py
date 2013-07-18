import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("ANASKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# ============= pre-setting ============================
# setup 'standard'  options
options = VarParsing.VarParsing ('standard')

# my own variable
options.register('inputType',
                 "MinBias900GeV",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Input file type - MB or QCDPtX")

# get and parse the command line arguments
options.parseArguments()


# =============== 7 TeV MC Sample =====================

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
    'file:/home/sungho/sctch101/mc/spectra/900GeV/test/QCDPtX_Recodebug_491_1_1BL.root'))

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/TrkAnaSkim_mc_900GeV_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('trkhistsMC.root')
                                   )

# =============== Import Sequences =====================

process.load("edwenger.Skims.EventFilter_cff")
process.load("edwenger.Skims.ExtraReco_cff")
process.load("edwenger.Skims.Analysis_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *

from edwenger.Skims.customise_cfi import *
enableSIM(process)    # activate isGEN in analyzers
setGlobTagAndRedigi(process,options.inputType) # this sets glob. tag and redigi name
updateEvtSelEff(process.trackAna_STD,"STD_NSD_900GeV_TrkVtx")
updateEvtSelEff(process.looseTrackAna_STD,"STD_NSD_900GeV_PixVtx")
enableDJetAna(process,"MC","LIGHT") # anaModes: "MC","Data", outLevels: "LIGHT","FF","FULL"
run900GeVmode(process) # set trigger names and appropriate binning for 900 GeV

# =============== Final Paths =====================

process.eventFilter_step = cms.Path(process.eventFilter)
process.extraReco_step   = cms.Path(process.eventFilter * process.extraReco)
process.ana_step         = cms.Path(process.eventFilter * process.analysisSeq)

# Re-make genjets for 35x reprocessing of 31x MC

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

process.gen_step = cms.Path(process.genJetParticles *
                            process.ak5GenJets *
                            process.ak7GenJets)

# =============== Output ================================
process.load("edwenger.Skims.analysisSkimContent_cff")
process.output = cms.OutputModule("PoolOutputModule",
    process.analysisSkimContent,
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AODSIM'),
      filterName = cms.untracked.string('TrkAnaFilter')),
    fileName = cms.untracked.string('trkAnaSkimAODSIM.root')
    )

process.output_step = cms.EndPath(process.output)

# =============== Schedule =====================

process.schedule = cms.Schedule(
    process.gen_step,
    process.eventFilter_step,
    process.extraReco_step,
    process.ana_step,
    #process.output_step
    )
