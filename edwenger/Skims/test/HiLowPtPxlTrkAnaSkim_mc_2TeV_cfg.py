import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# =============== 2.36 TeV MC Sample =====================

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
    '/store/mc/Spring10/Hydjet_Quenched_MinBias_2760GeV/GEN-SIM-RECO/MC_36Y_V7A-v1/0005/E82E5989-BE70-DF11-9952-001EC9AA9FCC.root'))

# =============== Other Statements =====================

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'MC_37Y_V4::All'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/edwenger/Skims/test/HiLowPtPxlTrkAnaSkim_mc_2TeV_cfg.py,v $'),
    annotation = cms.untracked.string('BPTX_AND + BSC_OR + !BSCHALO')
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('trkhistsMC.root')
                                   )

# =============== Import Sequences =====================
process.load("edwenger.Skims.HiAnalysis_cff")

from edwenger.Skims.hicustomise_cfi import *
process = enableSIM(process)    # activate isGEN in analyzers

# =============== Pat jet in HI ========================
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from Saved.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Spring10") # Spring10 is the new default

from Saved.JulyExercise.customise_cfi import *
#disableMC(process)

# =============== Low-Pt Pixel Tracking ==============

process.load("RecoLocalTracker/SiPixelRecHits/SiPixelRecHits_cfi")

process.load("RecoHI/HiTracking/LowPtTracking_PbPb_cff")

process.hiLowPtPixelTracks = process.hiPixel3PrimTracks.clone()
process.hiLowPtPixelTracks.RegionFactoryPSet.RegionPSet.ptMin = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.ptMin = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.nSigmaTipMaxTolerance = cms.double(4.0)
process.hiLowPtPixelTracks.FilterPSet.tipMax = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.nSigmaLipMaxTolerance = cms.double(4.0)
process.hiLowPtPixelTracks.FilterPSet.lipmax = cms.double(0.2)
process.hiLowPtPixelTracks.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000

process.trackAna.src = cms.untracked.InputTag("hiLowPtPixelTracks")
process.trackAna.src_evtCorr = cms.untracked.InputTag("hiLowPtPixelTracks")
process.trkEffAnalyzer.tracks = cms.untracked.InputTag("hiLowPtPixelTracks")
process.trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("hiLowPtPixelTracks")

# =============== Final Paths w/ Pixel Tracks =====================
process.extraReco_step   = cms.Path(
    process.siPixelRecHits *
    process.hiLowPtPixelTracks *
    process.makeHeavyIonJets
)
process.ana_step         = cms.Path(process.hiAnalysisSeq)


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
    process.extraReco_step,
    process.ana_step,
    #process.output_step
    )

