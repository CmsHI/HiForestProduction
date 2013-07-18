import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')


process.load("HeavyIonsAnalysis.Configuration.Sources.Data2010_Run150476_Express_HLT_MinBias_cff")

#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/cms/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/152/640/8604B150-C9FA-DF11-B748-0030487C912E.root")
#                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR10_P_V12::All'  #Data

process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


# --- Pat Stuff --
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

# redo tracker rechits
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoHI.Configuration.Reconstruction_HI_cff")
#process.load("RecoHI.HiTracking.LowPtTracking_PbPb_cff")


process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)

# Load Heavy Ion "Good" Track Selection 
process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")



# Define Jet Algo parameters
#process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')
#JPT jets
process.load("RecoJets.Configuration.RecoJPTJetsHIC_cff")

#Now the correct L1 
#process.load("RecoHI.HiJetAlgos.HiL1Corrector_cff")

process.patJets.embedCaloTowers = cms.bool(False)


#make pat jets for JPT jets
process.jpticPu5corr = process.patJetCorrFactors.clone(jetSource =
                                                       cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                       corrLevels =
                                                       cms.PSet(L2Relative = cms.string("L2Relative_AK5JPT"),
                                                                L3Absolute = cms.string("L3Absolute_AK5JPT"),
                                                                L5Flavor = cms.string("none")))
process.jpticPu5patJets = process.patJets.clone(jetSource  =
                                                cms.InputTag("JetPlusTrackZSPCorJetIconePu5"),
                                                jetCorrFactorsSource =
                                                cms.VInputTag(cms.InputTag("jpticPu5corr")))


process.JPTiterativeConePu5JetTracksAssociatorAtVertex.tracks = 'hiGoodMergedTracks'
process.JPTiterativeConePu5JetTracksAssociatorAtCaloFace.tracks = 'hiGoodMergedTracks'
process.JetPlusTrackZSPCorJetIconePu5.tracks = 'hiGoodMergedTracks'



process.runAllJets = cms.Sequence(process.recoJPTJetsHIC*process.jpticPu5corr*process.jpticPu5patJets) 





process.load("MNguyen.Configuration.HI_JetSkim_cff")
process.hltJetHI.HLTPaths = ["HLT_HIJet50U"]
process.hiEcalRecHitSpikeFilter.minEt = 20.0



process.jetSkimPath = cms.Path(
    process.rechits*
    process.jetSkimSequence*
    process.hiGoodTracksSelection*
    process.runAllJets
    )

process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")




process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContent,
                                  SelectEvents = cms.untracked.PSet(
                                      SelectEvents = cms.vstring('jetSkimPath')
                                    ),
                                  fileName = cms.untracked.string("RECOPAT.root"),
                                  dataset = cms.untracked.PSet(
                                      filterName = cms.untracked.string('jetSkimPath'),
                                  dataTier = cms.untracked.string('PAT')
                                  )
                                )





# Save some extra stuff
#?
process.output.outputCommands.extend(["keep *_heavyIon_*_*"])
# triggger
process.output.outputCommands.extend(["keep *_TriggerResults_*_*"])
#tracks
#process.output.outputCommands.extend(["keep *_hiSelectedTracks_*_*"])
process.output.outputCommands.extend(["keep *_hiGoodTracks_*_*"])
# reco jets
process.output.outputCommands.extend(["keep recoCaloJets_*_*_*"])
#calorimeter stuff
process.output.outputCommands.extend(["keep *_towerMaker_*_*"])
process.output.outputCommands.extend(["keep *_caloTowers_*_*"])
process.output.outputCommands.extend(["keep *_hcalnoise_*_*"])
process.output.outputCommands.extend(["keep *_hbhereco_*_*"])
process.output.outputCommands.extend(["keep *_horeco_*_*"])
process.output.outputCommands.extend(["keep *_hfreco_*_*"])
process.output.outputCommands.extend(["keep *_ecalRecHit_*_*"])
process.output.outputCommands.extend(["keep *_jptic*_*_*"])
process.output.outputCommands.extend(["keep *_recoJPT*_*_*"])
process.output.outputCommands.extend(["keep *_JetPlusTrack*_*_*"])




# Save all RECO!
#process.output.outputCommands.extend(["keep *_*_*_RECO"])

process.out_step = cms.EndPath(process.output)


# Schedule definition
process.schedule = cms.Schedule(process.jetSkimPath,process.out_step)


# And the logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring('Unknown', 
        'ProductNotFound', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ModuleFailure', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError', 
        'FatalRootError', 
        'NotFound')
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100









