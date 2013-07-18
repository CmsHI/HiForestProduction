import FWCore.ParameterSet.Config as cms

process = cms.Process('GENPARTANA')


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#from CmsHi.Analysis2010.CommonFunctions_cff import *
#overrideCentrality(process)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('particleSpecies.root')
)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring( 
#'/store/user/appeltel/mcrecoMIT/mnt/hadoop/cms/store/user/kjung/Pythia/RECO/VBatch3/pPbRecoShift_f1041_e.root'
#'/store/user/appeltel/mcrecoMIT/mnt/hadoop/cms/store/user/kjung/Hydjet/RECO/MIT2/pPbRecoShift_f331_e.root'
#'file:/scratch/appelte1/Hijing2_pPbRecoShift_f101_e.root'
#'file:/scratch/appelte1/hijing_pPbRecoShift_f991_e.root'
#'file:hydjetDefault_cfi.root'
#'file:/scratch/appelte1/HIJING2/27/pPbRecoShift_f181_e.root'
#'file:/scratch/appelte1/RECO_TEST_HIJING2_pPbRecoShift_f41_e.root'
'file:/scratch/appelte1/Hijing2_sameEventTestReco.root'
  )
)


process.GlobalTag.globaltag = 'START53_V14::All'

process.genPart = cms.EDAnalyzer('RpPbMCAnalyzer',
   genSrc = cms.InputTag("genParticles"),
   vertexSrc = cms.InputTag("offlinePrimaryVertices"),
   printDebug = cms.bool(True),
   etaMax = cms.double(2.5),
   rapidityShift = cms.double(0.0),
   vzMax = cms.double(999.0)
)

process.p = cms.Path( process.genPart )
