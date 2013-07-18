import FWCore.ParameterSet.Config as cms

process = cms.Process('HIJETS')



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'file:/tmp/mnguyen/SD_pp_JetHI.root'
    )

                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#load some general stuff
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR_P_V17::All'#GR_P_V14::All'

process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


#JEC from DB

process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource("PoolDBESSource",
                                                      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
            ),
                                                      timetype = cms.string('runnumber'),
                                                      toGet = cms.VPSet(
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                              tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_IC5Calo"),
                              label = cms.untracked.string("IC5Calo")
                              ),
           cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                 tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_AK3PF"),
                                 label = cms.untracked.string("AK3PFCand")
                                 ),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                  tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_AK4PF"),
                                  label = cms.untracked.string("AK4PFCand")
                                  ),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                  tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_AK5PF"),
                                  label = cms.untracked.string("AK5PFCand")
                                  ),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                 tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_PFTower_AK3PF"),
                                  label = cms.untracked.string("AK3PF")
                                  ),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                  tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_PFTower_AK4PF"),
                                  label = cms.untracked.string("AK4PF")
                                  ),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                  tag = cms.string("JetCorrectorParametersCollection_HIReco_413patch2_PFTower_AK5PF"),
                                  label = cms.untracked.string("AK5PF")
                                  ),
            ),
                           connect = cms.string("sqlite_file:JEC_HIReco_PFTower_413patch2_2011_v2.db"),
                           )


process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')



# --- Pat Stuff --
process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
#configureHeavyIons(process)

process.patJets.jetSource  = cms.InputTag("iterativeConePu5CaloJets")


process.patJets.addBTagInfo         = False
process.patJets.addTagInfos         = False
process.patJets.addDiscriminators   = False
process.patJets.addAssociatedTracks = False
process.patJets.addJetCharge        = False
process.patJets.addJetID            = False
process.patJets.getJetMCFlavour     = False
process.patJets.addGenPartonMatch   = False
process.patJets.addGenJetMatch      = False
process.patJets.embedGenJetMatch    = False
process.patJets.embedGenPartonMatch   = False




# Define Jet Algo parameters
process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')
process.load('RecoHI.HiJetAlgos.HiRecoPFJets_cff')
process.load("RecoJets.Configuration.RecoJPTJetsHIC_cff")

#Now the correct L1 
process.load("RecoHI.HiJetAlgos.HiL1Corrector_cff")


process.ic5sub = process.heavyIonL1SubtractedJets.clone()
process.ic5sub.src = 'ic5CaloJets'
process.ic5sub.rhoTag = 'kt4CaloJets'

process.patJets.embedCaloTowers = cms.bool(False)

# set up vertex matching, does it work for embedded events?
process.iterativeConePu5CaloJets.doPVCorrection = cms.bool(True)
process.iterativeConePu5CaloJets.srcPVs = 'hiSelectedVertex'
process.iterativeConePu5CaloJets.doPUOffsetCorr = False


process.ak3PFJets = process.ak5PFJets.clone()
process.ak3PFJets.rParam       = cms.double(0.3)


process.ak4PFJets = process.ak5PFJets.clone()
process.ak4PFJets.rParam       = cms.double(0.4)

    
process.iterativeConePu5CaloJets.puPtMin = cms.double(10.0)
process.ak5PFJets.puPtMin = cms.double(25.0)
process.ak4PFJets.puPtMin = cms.double(20.0)
process.ak3PFJets.puPtMin = cms.double(15.0)


# PF using a grid of pseudo-towers to emulate Russian-style subtraction

process.load("RecoHI.HiJetAlgos.ParticleTowerProducer_cff")
process.PFTowers = process.particleTowerProducer.clone()



process.icPu5corr = process.patJetCorrFactors.clone(src = cms.InputTag("iterativeConePu5CaloJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('IC5Calo')
                                                  )
process.icPu5patJets = process.patJets.clone(jetSource  = cms.InputTag("iterativeConePu5CaloJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("icPu5corr")))
process.icPu5patSequence = cms.Sequence(process.iterativeConePu5CaloJets*process.icPu5corr*process.icPu5patJets)




### AK5, PF  ###
process.ak5PFJets = process.ak5PFJets.clone()
process.ak5PFJets.src = 'PFTowers'
process.ak5PFJets.jetType = 'BasicJet'

process.ak5PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak5PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
)
process.ak5PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak5PFJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFcorr")))
process.ak5PFpatSequence = cms.Sequence(process.ak5PFJets*process.ak5PFcorr*process.ak5PFpatJets)



### AK4, PF  ###
process.ak4PFJets = process.ak4PFJets.clone()
process.ak4PFJets.src = 'PFTowers'
process.ak4PFJets.jetType = 'BasicJet'

process.ak4PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak4PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK4PF')
)
process.ak4PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak4PFJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4PFcorr")))
process.ak4PFpatSequence = cms.Sequence(process.ak4PFJets*process.ak4PFcorr*process.ak4PFpatJets)



### AK3, PF  ###
process.ak3PFJets = process.ak3PFJets.clone()
process.ak3PFJets.src = 'PFTowers'
process.ak3PFJets.jetType = 'BasicJet'

process.ak3PFcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak3PFJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK3PF')
)
process.ak3PFpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak3PFJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3PFcorr")))
process.ak3PFpatSequence = cms.Sequence(process.ak3PFJets*process.ak3PFcorr*process.ak3PFpatJets)



# Regular PF, no towers

### AK5, PF  ###

process.ak5PFCandJets = process.ak5PFJets.clone()
process.ak5PFCandJets.src = 'particleFlow'

process.ak5PFCandcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak5PFCandJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK5PF')
                                                    )
process.ak5PFCandpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak5PFCandJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFCandcorr")))
process.ak5PFCandpatSequence = cms.Sequence(process.ak5PFCandJets*process.ak5PFCandcorr*process.ak5PFCandpatJets)

### AK4, PF  ###

process.ak4PFCandJets = process.ak4PFJets.clone()
process.ak4PFCandJets.src = 'particleFlow'

process.ak4PFCandcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak4PFCandJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK4PF')
                                                    )
process.ak4PFCandpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak4PFCandJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak4PFCandcorr")))
process.ak4PFCandpatSequence = cms.Sequence(process.ak4PFCandJets*process.ak4PFCandcorr*process.ak4PFCandpatJets)

### AK3, PF  ###

process.ak3PFCandJets = process.ak3PFJets.clone()
process.ak3PFCandJets.src = 'particleFlow'

process.ak3PFCandcorr = process.patJetCorrFactors.clone(src = cms.InputTag("ak3PFCandJets"),
                                                    levels = cms.vstring('L2Relative','L3Absolute'),
                                                    payload = cms.string('AK3PF')
                                                    )
process.ak3PFCandpatJets = process.patJets.clone(jetSource  = cms.InputTag("ak3PFCandJets"),
                                             jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak3PFCandcorr")))
process.ak3PFCandpatSequence = cms.Sequence(process.ak3PFCandJets*process.ak3PFCandcorr*process.ak3PFCandpatJets)





process.runAllJets = cms.Sequence(
    process.icPu5patSequence +
    process.ak5PFpatSequence +
    process.ak4PFpatSequence +
    process.ak3PFpatSequence +
    process.ak5PFCandpatSequence +
    process.ak3PFCandpatSequence +
    process.ak4PFCandpatSequence
    )

# jet analysis trees
# trees need centrality
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


process.load("MNguyen.InclusiveJetAnalyzer.inclusiveJetAnalyzer_cff")
process.inclusiveJetAnalyzer.isMC = cms.untracked.bool(False)
process.inclusiveJetAnalyzer.useCentrality = cms.untracked.bool(False)
process.icPu5JetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.icPu5JetAnalyzer.jetTag = 'icPu5patJets'
process.ak5PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak5PFJetAnalyzer.jetTag = 'ak5PFpatJets'
process.ak4PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak4PFJetAnalyzer.jetTag = 'ak4PFpatJets'
process.ak3PFJetAnalyzer = process.inclusiveJetAnalyzer.clone()
process.ak3PFJetAnalyzer.jetTag = 'ak3PFpatJets'
process.ak3PFCandJetAnalyzer = process.ak3PFJetAnalyzer.clone()
process.ak4PFCandJetAnalyzer = process.ak4PFJetAnalyzer.clone()
process.ak5PFCandJetAnalyzer = process.ak5PFJetAnalyzer.clone()
process.ak3PFCandJetAnalyzer.jetTag = cms.InputTag("ak3PFCandpatJets")
process.ak4PFCandJetAnalyzer.jetTag = cms.InputTag("ak4PFCandpatJets")
process.ak5PFCandJetAnalyzer.jetTag = cms.InputTag("ak5PFCandpatJets")




process.inclusiveJetAnalyzerSequence = cms.Sequence(                  process.icPu5JetAnalyzer
                                                                      *process.ak5PFJetAnalyzer
                                                                      *process.ak4PFJetAnalyzer
                                                                      *process.ak3PFJetAnalyzer
                                                                      *process.ak3PFCandJetAnalyzer
                                                                      *process.ak4PFCandJetAnalyzer
                                                                      *process.ak5PFCandJetAnalyzer
                                                                      )

process.load("MNguyen.InclusiveJetAnalyzer.PFJetAnalyzer_cff")
process.PFJetAnalyzer.isMC = cms.untracked.bool(False)
process.PFJetAnalyzer.trackTag = cms.InputTag("hiGoodTracks")
process.PFJetAnalyzer.jetTag2 = cms.InputTag("ak5PFpatJets")
process.PFJetAnalyzer.recoJetTag2 = cms.InputTag("ak5PFJets")
process.PFJetAnalyzer.jetTag3 = cms.InputTag("ak4PFpatJets")
process.PFJetAnalyzer.recoJetTag3 = cms.InputTag("ak4PFJets")
process.PFJetAnalyzer.jetTag4 = cms.InputTag("ak3PFpatJets")
process.PFJetAnalyzer.recoJetTag4 = cms.InputTag("ak3PFJets")

#Frank's analyzer
from Saved.QM11Ana.Analyzers_cff import trkAnalyzer
process.trkAnalyzer = trkAnalyzer
process.trkAnalyzer.trackSrc = cms.InputTag("hiGoodTracks")


process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("JetAnalysisTTrees_hiGoodTracks_v1.root"))



#filters
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )

#process.load("MNguyen.Configuration.HI_JetSkim_cff")
#process.hltJetHI.HLTPaths = ["HLT_HIJet35U"]

process.jetSkimPath = cms.Path(
    process.scrapingVeto*
    process.primaryVertexFilter*
    process.HBHENoiseFilter*
    process.PFTowers*
    process.runAllJets*
    process.inclusiveJetAnalyzerSequence*
    process.PFJetAnalyzerSequence
    #*process.trkAnalyzer
    )




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

process.MessageLogger.cerr.FwkReport.reportEvery = 1









