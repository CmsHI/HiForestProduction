import FWCore.ParameterSet.Config as cms
process = cms.Process("ANALYSIS")

# Services
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_35X_V7A::All' #357
process.GlobalTag.globaltag ='GR_R_35X_V8B::All' #358p3

process.MessageLogger.cerr.FwkReport.reportEvery=100
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))


#Input source
process.source = cms.Source (
    "PoolSource",    
    fileNames = cms.untracked.vstring(

    # MinBias data (run 132440)
    #'/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0116/A2EC56ED-3A5C-DF11-AE04-002618943832.root',
    #'/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0110/966DE744-A959-DF11-949F-001A92971B90.root',
    #'/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0110/6CCAB9F0-AC59-DF11-A2F7-002618943948.root'

    # MinBias data (run 133874)
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0116/56853CE3-3A5C-DF11-9CF3-001A92811748.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0108/F27F337E-8359-DF11-8FD9-00261894394A.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FEB3A8ED-5E59-DF11-AF61-00261894394A.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FE48B9A6-6959-DF11-85A0-003048679012.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FE44168B-7259-DF11-85CE-003048D15CC0.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FC7EC61B-4F59-DF11-9832-001A92971B88.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FC5BAF07-8E59-DF11-92A2-0018F3D09600.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FC23349E-4B59-DF11-8E76-00248C0BE012.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FAF871D4-4F59-DF11-9D3B-0018F3D095F0.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FAF348CC-6259-DF11-A4A2-0026189438DF.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FA9D6FD8-4E59-DF11-B544-0018F3D09696.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FA4A5B04-8F59-DF11-8C1C-00261894389D.root',
    '/store/data/Commissioning10/MinimumBias/RECO/May6thReReco-v1/0005/FA352828-6859-DF11-8961-001A92810AC8.root',
    
    ),
    #lumisToProcess = cms.untracked.VLuminosityBlockRange('132440:157-132440:378'),
    lumisToProcess = cms.untracked.VLuminosityBlockRange('133874:166-133874:297','133874:299-133874:721','133874:724-133874:814','133874:817-133874:864','133875:1-133875:20','133875:22-133875:37','133876:1-133876:315','133877:1-133877:77','133877:82-133877:104','133877:113-133877:231','133877:236-133877:294','133877:297-133877:437','133877:439-133877:622','133877:624-133877:853','133877:857-133877:1472','133877:1474-133877:1640','133877:1643-133877:1931','133881:1-133881:71','133881:74-133881:223','133881:225-133881:551','133885:1-133885:132','133885:134-133885:728','133927:1-133927:44','133928:1-133928:645'),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Event selection
process.load("edwenger.Skims.EventFilter_cff")
process.load("edwenger.Skims.ExtraReco_cff")
process.eventFilter *= process.primaryVertexFilter

# Filter on PF candidate pt
process.clones = cms.EDFilter("PtMinCandViewCloneSelector",
    src=cms.InputTag("particleFlow",""),
    ptMin=cms.double(5.0),
    filter=cms.bool(True)
)
process.eventFilter *= process.clones


# track association
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *
TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

# PF candidate analyzer
process.pfCandidateAnalyzer = cms.EDAnalyzer("PFCandidateAnalyzer",
    PFCandidates = cms.InputTag("particleFlow",""),
    Vertices = cms.InputTag("sortedGoodVertices"),                                         
    verbose = cms.untracked.bool(False), ## print candidate info
    printBlocks = cms.untracked.bool(False), ## print block/element info
    ptMin = cms.untracked.double(5.0), ## of PF candidate
    SimTracks = cms.InputTag("mergedtruth","MergedTrackTruth"),
    Tracks = cms.InputTag("generalTracks"),
    hasSimInfo = cms.untracked.bool(False),
    minHits = cms.untracked.double(5),
    maxPtErr = cms.untracked.double(0.05),
    maxD0 = cms.untracked.double(0.2),                                         
    maxDZ = cms.untracked.double(0.2),
    maxD0Norm = cms.untracked.double(3.0),
    maxDZNorm = cms.untracked.double(3.0),
    pixelSeedOnly = cms.untracked.bool(True)                                         
)

# edm output
process.aod = cms.OutputModule("PoolOutputModule",
    outputCommands=cms.untracked.vstring('drop *',
                                         'keep recoPFClusters_*CAL_*_*',
                                         'keep recoPFBlocks_*_*_*', 
                                         'keep recoPFCandidates_*__*',
                                         'keep recoTracks_generalTracks_*_*'),
    fileName = cms.untracked.string('aod.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filter'))
)

# ntuple output
process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("pftuple.root"))

# Paths
process.filter = cms.Path(process.eventFilter)
process.p = cms.Path(process.eventFilter
                     *process.goodVertices
                     *process.sortedGoodVertices
                     *process.pfCandidateAnalyzer)
#process.outpath = cms.EndPath(process.aod)

# Schedule
process.schedule = cms.Schedule(process.filter,
                                process.p,
                                #process.outpath
                                )


def customiseMC(process):
    process.eventFilter.remove(process.physDeclFilter) # always false in MC
    process.eventFilter.remove(process.bptxAnd)        # always false in MC
    process.eventFilter.remove(process.bscOrBptxOr)
    process.eventFilter.remove(process.hltMinBias)
    process.TFileService.fileName="pftupleMC.root"     # new ntuple name
    process.source.fileNames= [                        # QCD Pt80 files
        '/store/mc/Spring10/QCD_Pt80/GEN-SIM-RECO/START3X_V26_S09-v1/0011/FC7F52EF-2F45-DF11-869F-E41F13181D30.root',
        '/store/mc/Spring10/QCD_Pt80/GEN-SIM-RECO/START3X_V26_S09-v1/0011/FC5048AD-3245-DF11-ABFD-E41F131817E8.root',
        '/store/mc/Spring10/QCD_Pt80/GEN-SIM-RECO/START3X_V26_S09-v1/0011/F8D06B88-3045-DF11-ACE5-00215E21D6DE.root',
        '/store/mc/Spring10/QCD_Pt80/GEN-SIM-RECO/START3X_V26_S09-v1/0011/F2D274EA-2D45-DF11-8E2B-E41F13181A74.root',
        '/store/mc/Spring10/QCD_Pt80/GEN-SIM-RECO/START3X_V26_S09-v1/0011/F23800E9-2E45-DF11-ACCF-E41F1318162C.root'
        ]
    
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('1:1-999999:max')
    
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.pfCandidateAnalyzer.hasSimInfo=cms.untracked.bool(True)
    
    return process


def customiseJet15U(process):
    process.minBiasBscFilter.remove(process.bscOrBptxOr)
    process.minBiasBscFilter.remove(process.hltMinBias)
    process.hltJets.HLTPaths = ['HLT_Jet15U']
    process.minBiasBscFilter *= process.hltJets
    
    process.source.fileNames= [
        # JetMetTau PD
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/C6338C7C-7466-DF11-B48C-0030487A322E.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/B6E9C8D6-6E66-DF11-B497-001D09F23A3E.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/AAA8A40E-9466-DF11-A2B2-001D09F2AF96.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/7062F157-7966-DF11-AD70-000423D98E54.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/4A6D37F5-7E66-DF11-948F-001D09F244BB.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/400BD94D-8C66-DF11-8C50-000423D95220.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/3AFD28CD-6766-DF11-837B-000423D33970.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/18E30CD7-5B66-DF11-A361-001D09F29538.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/122E50F4-7E66-DF11-809E-001D09F24DDA.root',
        '/store/data/Run2010A/JetMETTau/RECO/v1/000/136/066/009A51E5-6966-DF11-92AC-0030487CD812.root'
        ]
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('136066:181-136066:297','136066:299-136066:319','136066:321-136066:321','136066:323-136066:336','136066:338-136066:348','136066:350-136066:366','136066:368-136066:529','136066:532-136066:544','136066:546-136066:595','136066:597-136066:719','136066:721-136066:1004','136066:1006-136066:1184','136080:249-136080:250','136080:252-136080:256','136082:1-136082:173','136082:175-136082:422','136082:477-136082:504','136082:506-136082:506')
    
    
    return process

def customiseLooseCuts(process):
    process.pfCandidateAnalyzer.minHits=0
    process.pfCandidateAnalyzer.maxPtErr=1.0
    process.pfCandidateAnalyzer.maxD0=10.0
    process.pfCandidateAnalyzer.maxDZ=10.0
    process.pfCandidateAnalyzer.maxD0Norm=100.0
    process.pfCandidateAnalyzer.maxDZNorm=100.0
    process.pfCandidateAnalyzer.pixelSeedOnly=False

    return process


#process = customiseMC(process)
process = customiseJet15U(process)
process = customiseLooseCuts(process)
