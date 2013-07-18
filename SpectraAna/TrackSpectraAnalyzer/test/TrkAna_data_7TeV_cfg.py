import FWCore.ParameterSet.Config as cms
process = cms.Process("ANA")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 7 TeV Collisions =====================

#process.load("edwenger.Skims.fileNames7TeV_cff")
#process.source = cms.Source("PoolSource",process.PromptReco7TeV)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/user/s/sungho/cmssw/CMSSW_356_SpectraPPV4/src/edwenger/Skims/python/trkAnaSkimAODSIM.root'),
)



# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_R_35X_V6::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/ASYoon/SpectraAna/TrackSpectraAnalyzer/test/TrkAna_data_7TeV_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)


# =============== Extra Reco Steps =====================
#process.load("edwenger.Skims.BeamSpot7TeV_cff")     # custom beamspot db source
process.load("edwenger.Skims.ChargedCandidates_cff") # make charged candidates from selected tracks
process.load("edwenger.Skims.ExtraVertex_cff")       # agglomerative pixel vertexing
process.load("edwenger.Skims.TrackRefit_cff")        # refit constrained to primary vertex

#================ Track Ana ==========================
process.load("SpectraAna.TrackSpectraAnalyzer.trackspectraanalyzer_cfi")
process.trackAna.isGEN = cms.untracked.bool(False)
process.trackAna.etaMax = cms.untracked.double(2.5)
process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('trackAnaDATA.root')
                                   )



# =============== Final Paths =====================
process.ana = cms.Path(process.chargedCandidates *
                       process.trackAna)

