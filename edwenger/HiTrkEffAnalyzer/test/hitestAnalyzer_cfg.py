import FWCore.ParameterSet.Config as cms

process = cms.Process('HiTEST')

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
process.MessageLogger.categories = ['HiTrkEffAnalyzer']
process.MessageLogger.debugModules = ['*']
process.MessageLogger.cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
            DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
                ),
            INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
                ),
        HiTrkEffAnalyzer = cms.untracked.PSet(
          limit = cms.untracked.int32(-1)
          )
)

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.3 $'),
            annotation = cms.untracked.string('step2 nevts:1'),
            name = cms.untracked.string('PyReleaseValidation')
        )

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(10)
        )
process.options = cms.untracked.PSet(
        #wantSummary = cms.untracked.bool(True)
    )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'dcache:/pnfs/cmsaf.mit.edu/hibat/cms/users/jazzitup/relval/CMSSW_3_8_4/RECO/dijet30/gen-sim-reco-embedded_1.root')
    'file:/net/hisrv0001/home/y_alive/cmssw/CMSSW_384p2_HiTracking/src/RecoHI/HiTracking/test/test_out_numEvent5.root')
)

# Track selection 
#process.load("edwenger.HiTrkEffAnalyzer.HiTPCuts_cff")

# Track efficiency analyzer 
process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")

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
process.GlobalTag.globaltag = 'MC_38Y_V8::All' #38X

# Path and EndPath definitions
process.ana_step  = cms.Path(process.hitrkEffAna)
                             #*process.hipxltrkEffAna)

#process.out_step = cms.EndPath(process.output)


