import FWCore.ParameterSet.Config as cms
process = cms.Process("ANA")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 7 TeV Collisions =====================

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/user/s/sungho/cmssw/CMSSW_356_SpectraPPV4/src/edwenger/Skims/python/trkAnaSkimAODSIM.root'),
    #'rfio:/castor/cern.ch/user/s/sungho/spectra/7TeV/mc/mtv_reco_ana/v4_10M/trkAnaSkimAODSIM_98_1.root'),
    #'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root',
    #'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FED8673E-F53D-DF11-9E58-0026189437EB.root'),
)

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.GlobalTag.globaltag = 'START3X_V26A::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/ASYoon/SpectraAna/TrackSpectraAnalyzer/test/TrkAna_mc_7TeV_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)

# =============== Extra Reco Steps =====================
process.load("edwenger.Skims.ExtraVertex_cff")       # agglomerative pixel vertexing
process.load("edwenger.Skims.ChargedCandidates_cff") # make charged candidates from selected tracks

#================ Track Association =====================
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
#process.trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("selectTracks")
#process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

#================ Fake and Real track producer ==========
#process.load("PPTrackingTools.FakeAndRealTrackSelector.selectFakeAndReal_cff")


#================= MTV for eff/fake rate ================
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("Validation.RecoTrack.cuts_cff")
process.load("CmsHi.TrackAnalysis.findableSimTracks_cfi")
process.load("edwenger.Skims.MTVConfiguration_cff")

# for eff
process.primaryChgSimTracks.minHit = cms.int32(0)

# some of the cuts are already applied in chargedCandidates
process.cutsRecoTracks.quality = ['highPurity']
process.cutsRecoTracks.source = cms.InputTag('selectTracks')
process.cutsRecoTracks.minHit   = cms.int32(8)
process.cutsRecoTracks.minRapidity  = cms.double(-2.5)
process.cutsRecoTracks.maxRapidity  = cms.double(2.5)
process.cutsRecoTracks.ptMin = cms.double(0.1)
process.cutsRecoTracks.tip = cms.double(3.5)
process.cutsRecoTracks.lip = cms.double(30)

# for fake
process.cutsTPFake.minRapidity = cms.double(-2.5)
process.cutsTPFake.maxRapidity = cms.double(2.5)
process.cutsTPFake.tip = cms.double(3.5)
process.cutsTPFake.lip = cms.double(30)
process.cutsTPFake.ptMin = cms.double(0.1)
process.cutsTPFake.minHit   = cms.int32(0)

# mtv config
process.multiTrackValidator.outputFile = 'mtv.root'
process.multiTrackValidator.label_tp_effic = cms.InputTag("primaryChgSimTracks")
process.multiTrackValidator.label_tp_fake  = cms.InputTag("cutsTPFake")
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.label = ['cutsRecoTracks']
process.multiTrackValidator.UseAssociators = cms.bool(False)

process.mtv = cms.Sequence(
    process.primaryChgSimTracks *
    process.cutsTPFake  *
    process.cutsRecoTracks *
    process.multiTrackValidator
)

#================ Track Ana ==========================
process.load("SpectraAna.TrackSpectraAnalyzer.trackspectraanalyzer_cfi")
process.trackAna.etaMax = cms.untracked.double(2.5)
#process.trackAna.doJet = cms.untracked.bool(False)
process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('trackAnaMC.root')
                                   )

# =============== Final Paths =====================
process.trkAna = cms.Path(process.chargedCandidates *
                          #process.extraVertex *
                          #process.mtv *
                          process.trackAna)
