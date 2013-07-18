import FWCore.ParameterSet.Config as cms

process = cms.Process("Vertexing")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("UserCode.FerencSiklerVertexing.NewVertexProducer_cfi")

###############################################################################
# Categories and modules
process.CategoriesAndModules = cms.PSet(
    categories = cms.untracked.vstring('NewVertices'),
    debugModules = cms.untracked.vstring('*')
)

###############################################################################
# Message logger
process.MessageLogger = cms.Service("MessageLogger",
    process.CategoriesAndModules,
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    destinations = cms.untracked.vstring('cerr')
)

###############################################################################
# Source
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
        '/store/relval/CMSSW_3_1_0_pre6/RelValMinBias/GEN-SIM-RECO/IDEAL_31X_v1/0002/6A330844-1733-DE11-974A-001617E30F50.root',
       '/store/relval/CMSSW_3_1_0_pre6/RelValMinBias/GEN-SIM-RECO/IDEAL_31X_v1/0002/1E1124BF-8032-DE11-AAA7-000423D944FC.root')
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/F64ABC78-1733-DE11-866C-000423D98E6C.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/9A2FD94C-2033-DE11-B10A-001617DF785A.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/98E8674A-FB32-DE11-BD4A-000423D99896.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/8823A234-FC32-DE11-A251-001617E30F58.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/6E122CEB-FA32-DE11-9055-000423D99264.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/50D4BADB-FA32-DE11-BA01-000423D98DC4.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/406012BC-FB32-DE11-AD0F-000423D94990.root',
#        '/store/relval/CMSSW_3_1_0_pre6/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0002/12B74DD2-F932-DE11-921C-0016177CA7A0.root')
)

process.GlobalTag.globaltag = 'IDEAL_31X::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

###############################################################################
# Vertex finder
#process.newVertices.Method = cms.int32(0)

###############################################################################
# Path
process.vtxs = cms.Path(process.newVertices)

# Schedule
process.schedule = cms.Schedule(process.vtxs)
