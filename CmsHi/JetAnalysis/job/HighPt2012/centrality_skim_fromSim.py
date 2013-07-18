# Auto generated configuration file
# using: 
# Revision: 1.341.2.2 
# Source: /local/reps/CMSSW.admin/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: skim --conditions STARTHI44_V7::All --scenario HeavyIons -n 10 --eventcontent FEVTDEBUGHLT -s DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --inputCommands keep *,drop *_simEcalPreshowerDigis_*_* --filein file:step1.root --fileout file:skim.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIMSKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:/mnt/hadoop/cms/store/mc/Summer12/HidjetQuenchedMinBias/GEN-SIM/STARTHI44_V12-v5/0001/309FAE25-29CD-E111-BD31-001517E73D34.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('skim nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Reconstruction Sequences
# Global + High-Level Reco Sequence
process.globalRecoPbPb = cms.Sequence(
                              process.hiEcalClusters
                              * process.hiRecoJets
                              * process.hiCentrality
                              )


process.hiCentrality.producePixelhits = cms.bool(False)
process.hiCentrality.produceTracks = cms.bool(False)
process.hiCentrality.producePixelTracks = cms.bool(False)

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
	nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
	centralitySrc = cms.InputTag("hiCentrality")
	)
process.hiCentrality.pixelBarrelOnly = False


#####################################################################################
# Event Filter
#####################################################################################
process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
process.centralityFilter.selectedBins = range(0,12)
process.superFilterSequence = cms.Sequence(process.centralityFilter)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('reconstruction_step')),
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('skim.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-HLTDEBUG')
    )
)
process.FEVTDEBUGHLToutput.outputCommands = ['drop *','keep *_*_*_SIM','keep *_hiCentrality_*_*']

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'STARTHI44_V7::All'

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons*process.centralityFilter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
