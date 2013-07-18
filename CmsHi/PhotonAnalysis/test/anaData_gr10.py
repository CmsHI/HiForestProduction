
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.Timing = cms.Service("Timing")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

## Source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_100_1_Gnr.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_101_1_byR.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_102_1_NvE.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_103_1_QZM.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_104_1_G8d.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_105_1_AO7.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_106_1_lJ7.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_107_1_g6R.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_108_1_4Pf.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_109_1_ntv.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_10_1_QDZ.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_159_1_JtW.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_15_1_WdS.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_160_1_8WU.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_161_1_Zqg.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_162_1_TZk.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_163_1_kV3.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_164_1_qeE.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_165_1_TgB.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_166_1_j9h.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_167_1_BA6.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_168_1_7jL.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_169_1_Ugi.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_178_1_wLy.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_179_1_Veh.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_17_1_Mkl.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_180_1_Hz4.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_181_1_fgN.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_182_1_shE.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_183_1_tE0.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_184_1_euC.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_185_1_zjH.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_186_1_5pL.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_187_1_J1X.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_188_1_eH2.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_189_1_zDB.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_18_1_3nO.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_190_1_Qhy.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_191_1_7mY.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_192_1_kqi.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_193_1_wmG.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_194_1_dti.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_195_1_K8N.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_196_1_IG8.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_197_1_obm.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_198_1_KJk.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_199_1_SRJ.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_19_1_Ltf.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_230_1_B53.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_231_1_rZd.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_232_1_zKE.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_233_1_B4L.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_234_1_Yi9.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_235_1_hn3.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_236_1_dqo.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_237_1_voF.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_238_1_ARV.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_239_1_sLp.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_23_1_aIL.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_240_1_VgL.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_241_1_Npd.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_242_1_yeZ.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_243_1_Tx4.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_244_1_OMP.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_245_1_nEv.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_246_1_2eq.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_288_1_OwP.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_289_1_hTS.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_28_1_3G1.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_290_1_snJ.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_291_1_PmO.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_292_1_Nzk.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_293_1_t4N.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_294_1_zU9.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_295_1_Jhn.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_296_1_dlb.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_297_1_Zm0.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_298_1_snR.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_299_1_Z9J.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_29_1_4vv.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_2_1_lLc.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_300_1_qrB.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_301_1_a4N.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_302_1_CVZ.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_303_1_ctT.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_304_1_jZO.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_305_1_D3O.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_306_1_Lk3.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_307_1_cIH.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_308_1_X5q.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_309_1_ZYS.root',
'/store/caf/user/kimy/bscMinbiasSkim_run150476AllPhysicsReHlt_v003/Hydjet_Quenched_MinBias_2760GeV_cfi_py_REHLT_RECO_30_1_SP9.root'

# lumisToProcess =  cms.untracked.VLuminosityBlockRange(
# '150431:1-150431:1000'
)
)



## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("RecoHI.HiEgammaAlgos.HiEgamma_cff")


# Centrality
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.HeavyIonGlobalParameters = cms.PSet(
       centralitySrc = cms.InputTag("hiCentrality"),
          centralityVariable = cms.string("HFhits"),
          nonDefaultGlauberModel = cms.string("")
          )
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR10_P_V12::All')  # for data global run.
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)
from PhysicsTools.PatAlgos.tools.coreTools import *
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# Modification for HI
process.load("CmsHi.PhotonAnalysis.MultiPhotonAnalyzer_cfi")
process.multiPhotonAnalyzer.GenParticleProducer = cms.InputTag("hiGenParticles")
process.multiPhotonAnalyzer.PhotonProducer = cms.InputTag("selectedPatPhotons")
process.multiPhotonAnalyzer.VertexProducer = cms.InputTag("hiSelectedVertex")
process.multiPhotonAnalyzer.doStoreMET = cms.untracked.bool(False)
process.multiPhotonAnalyzer.doStoreJets = cms.untracked.bool(False)
process.multiPhotonAnalyzer.OutputFile = cms.string('mpa_hfSkim_icludingSpikes.root')
#process.multiPhotonAnalyzer.isMC_      = cms.untracked.bool(False)
#process.singlePhotonAnalyzer.isMC_      = cms.untracked.bool(False)

#  Not embedding sample..
process.multiPhotonAnalyzer.GenEventScale = cms.InputTag("generator")   

#
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('___outf2___.root'),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


# Comp Cone Analysis
process.multiPhotonAnalyzer.doStoreCompCone = cms.untracked.bool(False)
process.photons.isolationSumsCalculatorSet.trackProducer = process.multiPhotonAnalyzer.TrackProducer

 # spike cleaner
process.patPhotons.photonSource = cms.InputTag("cleanPhotons")
process.photonMatch.src = cms.InputTag("cleanPhotons")
process.PhotonIDProd.photonProducer  = cms.string('cleanPhotons')

# detector responce
process.load("CmsHi.PhotonAnalysis.isoConeInspector_cfi")
process.load("CmsHi.PhotonAnalysis.ecalHistProducer_cfi")
process.load("CmsHi.PhotonAnalysis.SpikeInspector_cfi")

# high purity tracks
process.highPurityTracks = cms.EDFilter("TrackSelector",
                                        src = cms.InputTag("hiSelectedTracks"),
                                        cut = cms.string('quality("highPurity")')
                                        )

# get the pixel tree
process.pixelVertexFromClusters = cms.EDProducer('PixelVertexProducerClusters')
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("AMPT_2760GeV"),
    centralitySrc = cms.InputTag("hiCentrality")
    )


process.anaPixelHit = cms.EDAnalyzer('PixelHitAnalyzer',
                                     vertexSrc = cms.vstring('hiSelectedVertex'),
                                     trackSrc = cms.InputTag('hiSelectedTracks'),
                                     doTracking = cms.untracked.bool(False),
                                     L1gtReadout = cms.InputTag('gtDigis'),
                                     hltTrgResults = cms.untracked.string("TriggerResults::HLT"),
                                     doMC        = cms.untracked.bool(False),
                                     doPixel     = cms.untracked.bool(False)
                                     )



# for data
process.patHeavyIonDefaultSequence.remove(process.heavyIon)
process.patHeavyIonDefaultSequence.remove(process.genPartons)
process.patHeavyIonDefaultSequence.remove(process.hiPartons)

process.patHeavyIonDefalutSequencePhoton = cms.Sequence(process.isoCC1+process.isoCC2+process.isoCC3+process.isoCC4+process.isoCC5+process.isoCR1+process.isoCR2+process.isoCR3+process.isoCR4+process.isoCR5+process.isoT11+process.isoT12+process.isoT13+process.isoT14+process.isoT21+process.isoT22+process.isoT23+process.isoT24+process.isoT31+process.isoT32+process.isoT33+process.isoT34+process.isoT41+process.isoT42+process.isoT43+process.isoT44+process.isoDR11+process.isoDR12+process.isoDR13+process.isoDR14+process.isoDR21+process.isoDR22+process.isoDR23+process.isoDR24+process.isoDR31+process.isoDR32+process.isoDR33+process.isoDR34+process.isoDR41+process.isoDR42+process.isoDR43+process.isoDR44+process.PhotonIDProd+process.interestingEcalDetIdEB*process.interestingEcalDetIdEE*process.reducedEcalRecHitsEB*process.reducedEcalRecHitsEE*process.selectDigi*process.gamIsoDepositTk+process.gamIsoDepositEcalFromHits+process.gamIsoDepositHcalFromTowers+process.gamIsoDepositHcalDepth1FromTowers+process.gamIsoDepositHcalDepth2FromTowers*process.gamIsoFromDepsTk*process.gamIsoFromDepsEcalFromHits*process.gamIsoFromDepsHcalFromTowers*process.gamIsoFromDepsHcalDepth1FromTowers*process.gamIsoFromDepsHcalDepth2FromTowers*process.gamIsoDepositTk*process.gamIsoFromDepsTk+process.gamIsoDepositEcalFromHits*process.gamIsoFromDepsEcalFromHits+process.gamIsoDepositHcalFromTowers*process.gamIsoFromDepsHcalFromTowers*process.patPhotons*process.selectedPatPhotons)
process.interestingEcalDetIdEB.basicClustersLabel = cms.InputTag("islandBasicClusters","islandBarrelBasicClusters")
process.interestingEcalDetIdEE.basicClustersLabel = cms.InputTag("islandBasicClusters","islandEndcapBasicClusters")

# no cut at this moment
process.hiSpikeCleanedSC.TimingCut = cms.untracked.double(4.0)
process.hiSpikeCleanedSC.swissCutThr = cms.untracked.double(0.9)



import HLTrigger.HLTfilters.hltHighLevel_cfi
process.HIMinBiasBSC = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.HIMinBiasBSC.HLTPaths = ["HLT_HIMinBiasBSC"]
process.HIMinBiasBSC.TriggerResultsTag = cms.InputTag("TriggerResults","","REHLT")
process.HIphotontrig = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.HIphotontrig.HLTPaths = ["HLT_HIPhoton15"]
process.HIphotontrig.TriggerResultsTag = cms.InputTag("TriggerResults","","REHLT")

process.p = cms.Path(
    process.HIMinBiasBSC *
    process.HIphotontrig *
    process.highPurityTracks *
    #  process.siPixelRecHits*process.anaPixelHit*
    process.hiPhotonCleaningSequence*
    process.photonIDSequence*
    # process.patHeavyIonDefaultSequence *
    process.patHeavyIonDefalutSequencePhoton*
    process.multiPhotonAnalyzer *
    #    process.isoConeMap  *
    #    process.ecalHistProducer *
    process.spikeInspector
    )

#process.e = cms.EndPath(process.out)    

