import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    annotation = cms.untracked.string('test nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    #RAW
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0046/F0315553-BCD3-DF11-AFC7-00261894387B.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/FE0B79FD-6FD3-DF11-A3A7-003048678F8E.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/B43EC98F-A4D3-DF11-BC5E-00304867920C.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/9049424A-8AD3-DF11-894C-002618943910.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/7CAD6267-70D3-DF11-A92C-002618FDA279.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/7280F344-84D3-DF11-A2AC-002618943856.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/5AFD0B99-76D3-DF11-87CA-001A92810AE2.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/3AEF04F5-6DD3-DF11-B3CA-0018F3D095EC.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/325ED9FE-6ED3-DF11-8020-002618943910.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/2020B1F1-6CD3-DF11-8F19-0018F3D09616.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/168D9A63-7CD3-DF11-A784-00261894387B.root',
    '/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0044/1066E12C-94D3-DF11-AF97-002618943910.root'

    #RECO
    #'/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-RECO/MC_39Y_V2-v1/0045/E0991CAE-B5D3-DF11-A19E-0026189438E4.root',
    #'/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-RECO/MC_39Y_V2-v1/0044/CA67894B-8AD3-DF11-A223-00261894390B.root',
    #'/store/relval/CMSSW_3_9_0_pre7/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-RECO/MC_39Y_V2-v1/0044/7E7C82FE-6FD3-DF11-83AB-0018F3D09608.root'
    )
)

# Additional output definition

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('trkhists.root')
                                   )

# Other statements
process.GlobalTag.globaltag = 'MC_39Y_V2::All'

#process.hiGlobalPrimTracks.Fitter = 'KFFittingSmoother' # without outlier rejection

#process.trkana = cms.EDAnalyzer("OnTrackClusterAnalyzer")
process.trkana = cms.EDAnalyzer("OnTrackClusterAnalyzer",
                                src = cms.untracked.InputTag("hiSelectedTracks"),
                                setQualityBit = cms.untracked.bool(True)
                                )

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
process.ana_step = cms.Path(process.trkana)

# Schedule definition
#process.schedule = cms.Schedule(process.ana_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.ana_step)
