import FWCore.ParameterSet.VarParsing as VarParsing

import FWCore.ParameterSet.Config as cms

process = cms.Process('hiForestAna2011')

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)


hiTrackQuality = "highPurity"              # iterative tracks
#hiTrackQuality = "highPuritySetWithPV"    # calo-matched tracks

doElectrons = False
doRegitForBjets = False


#####################################################################################

process.load("CmsHi.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V2 for pPb",
                                          "PF : generalTracks (pp) with highPurity",
                                          "EP Flattening OFF",
                                          "Electrons OFF",
                                          "Preshower OFF"
                                          )

#####################################################################################
# Input source
#####################################################################################
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
'file:/d101/icali/pARun/Run201355/Data.00210355.0001.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0002.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0003.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0004.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0005.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0006.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0007.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0008.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0009.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0010.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0011.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0012.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0013.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0014.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0015.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0016.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0017.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0018.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0019.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0020.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0021.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0022.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0023.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0024.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0025.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0026.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0027.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0028.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0029.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0030.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0031.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0032.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0033.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0034.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0035.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0036.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0037.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0038.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0039.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0040.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0041.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0042.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0043.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0044.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0045.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0046.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0047.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0048.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0050.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0052.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0057.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0062.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0064.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0065.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0069.A.storageManager.00.0000.dat',
'file:/d101/icali/pARun/Run201355/Data.00210355.0074.A.storageManager.00.0000.dat',
    ))    

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1))


#####################################################################################
# Load some general stuff
#####################################################################################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('Configuration.StandardSequences.RawToDigi_cff')

process.GlobalTag.globaltag = 'GR_E_V33::All'


# load centrality

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality")
)
        
#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("pPb_openHLT.root"))

process.load('L1Trigger.Configuration.L1Extra_cff')
process.load('CmsHi.HiHLTAlgos.hltanalysis_cff')

process.hltanalysis.dummyBranches = cms.untracked.vstring()

process.hltAna = cms.EndPath(process.hltanalysis)
process.reco_extra=cms.Path(process.RawToDigi*process.L1Extra)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.reco_extra,
    process.hltAna
    )
