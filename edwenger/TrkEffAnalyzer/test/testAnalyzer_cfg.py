# Auto generated configuration file
# using: 
# Revision: 1.168.2.1 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 -s RAW2DIGI,RECO --datatier GEN-SIM-RECO --eventcontent RECODEBUG --conditions auto:mc --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

# message logger
process.MessageLogger.categories = ['TrkEffAnalyzer']
process.MessageLogger.debugModules = ['*']
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    DEBUG = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    TrkEffAnalyzer = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    )
)

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.6 $'),
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.options = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True)
)

# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
    '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root',
    '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FED8673E-F53D-DF11-9E58-0026189437EB.root'),
)

# Track selection
process.highPurityTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('quality("highPurity")')
)

# Track efficiency analyzer
process.load("edwenger.TrkEffAnalyzer.trkEffAnalyzer_cff")
process.trkEffAnalyzer.tracks = cms.untracked.InputTag('highPurityTracks')
process.trackingParticleRecoTrackAsssociation.label_tr = cms.untracked.InputTag('highPurityTracks')
process.trkEffAnalyzer.doAssociation = cms.untracked.bool(True) # association done inside analyzer (not from map)
process.trkEffAnalyzer.hasSimInfo = cms.untracked.bool(True) # for MC matching

# Output definition
process.RECODEBUGEventContent.outputCommands.extend(
    ['keep *_trackingParticleRecoTrackAsssociation_*_*'] )
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECODEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('output.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('hists.root')
                                   )

# Other statements
process.GlobalTag.globaltag = 'MC_3XY_V27::All' #357
#process.GlobalTag.globaltag = 'MC_36Y_V6::All'   #360

# Path and EndPath definitions
process.ana_step = cms.Path(process.highPurityTracks *
                            #process.trackingParticleRecoTrackAsssociation *
                            process.trkEffAnalyzer)
#process.out_step = cms.EndPath(process.output)

