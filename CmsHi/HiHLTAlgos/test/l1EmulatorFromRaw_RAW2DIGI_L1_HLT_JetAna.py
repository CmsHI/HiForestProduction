# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: l1EmulatorFromRaw -s RAW2DIGI,L1,HLT:HIon -n 100 --conditions auto:hltonline --datatier DIGI-RECO --eventcontent FEVTDEBUGHLT --data --filein /store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/150/590/02D89852-5DEC-DF11-9E11-001D09F282F5.root --customise=L1Trigger/Configuration/customise_l1EmulatorFromRaw --custom_conditions=L1GtTriggerMenu_L1Menu_Collisions2011_v5,L1GtTriggerMenuRcd,frontier://FrontierProd/CMS_COND_31X_L1T --processName=L1EmulRawHLT --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('L1EmulRawHLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('HLTrigger.Configuration.HLT_HIon_JetMET_data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:/mnt/hadoop/cms/store/results/heavy-ions/HICorePhysics/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/HICorePhysics/USER/StoreResults-HICorePhysics_Skim_MinimumBias_RAW-a606dc809a29a92e17749e5652319ad0-SD_MBHI/0000/0EDFD1BD-C9AE-E011-963E-003048CB82AC.root')
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.6 $'),
    annotation = cms.untracked.string('l1EmulatorFromRaw nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('l1EmulatorFromRaw_RAW2DIGI_L1_HLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DIGI-RECO')
    )
)

# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'GR_R_44_V5::All'
process.GlobalTag.globaltag = 'START44_V6::All'
# customise globaltag
process.GlobalTag.toGet = cms.VPSet()
process.GlobalTag.toGet.append(cms.PSet(tag=cms.string("L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2011_v0_mc"),record=cms.string("L1GtTriggerMenuRcd"),connect=cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T"),))

# Analyzers
# L1 Ntuple
process.load("UserCode.L1TriggerDPG.l1NtupleProducer_cfi")
process.l1NtupleProducer.gctTauJetsSource     = cms.InputTag("none","")
process.l1NtupleProducer.gctCentralJetsSource = cms.InputTag("none","")
process.l1NtupleProducer.gctNonIsoEmSource    = cms.InputTag("none","")
process.l1NtupleProducer.gctForwardJetsSource = cms.InputTag("none","")
process.l1NtupleProducer.gctIsoEmSource       = cms.InputTag("none","")
process.l1NtupleProducer.gctEnergySumsSource  = cms.InputTag("none","")
process.l1NtupleProducer.gctTauJetsSource     = cms.InputTag("none","")
process.l1NtupleProducer.rctSource            = cms.InputTag("none")
process.l1NtupleProducer.dttfSource           = cms.InputTag("none")
process.l1NtupleProducer.csctfTrkSource       = cms.InputTag("none")
process.l1NtupleProducer.csctfLCTSource       = cms.InputTag("none")
process.l1NtupleProducer.csctfStatusSource    = cms.InputTag("none")
process.l1NtupleProducer.csctfDTStubsSource   = cms.InputTag("none")
process.l1NtupleProducer.gtSource             = cms.InputTag("simGtDigis") # put gtDigis for data
process.l1NtupleProducer.gmtSource            = cms.InputTag("simGmtDigis") # put gtDigis (yes, not gmtDigis) for data

# openhlt
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.UseTFileService			= cms.untracked.bool(True)

process.hltbitnew = process.hltbitanalysis.clone(
  l1GtObjectMapRecord	= cms.InputTag("hltL1GtObjectMap","",process.name_()),
  l1GtReadoutRecord		= cms.InputTag("simGtDigis","",process.name_()),
  hltresults           = cms.InputTag("TriggerResults","",process.name_()),
  HLTProcessName       = cms.string(process.name_())
  )

# Jet Analyzers
process.load("PhysicsTools.PatAlgos.patHeavyIonSequences_cff")
process.patJets.jetSource  = cms.InputTag("hltHICaloJetCorrected")
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
process.patJets.embedGenPartonMatch = False
process.patJets.embedCaloTowers	    = False
process.patJets.addJetCorrFactors		= False

from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *
process.icPu5JetAnalyzer = inclusiveJetAnalyzer.clone(
	L1gtReadout = 'hltGtDigis::'+process.name_(),
	hltTrgResults = cms.untracked.string('TriggerResults::'+process.name_()),
	jetTag = 'patJets',
	hltTrgNames = ['HLT_HIJet55_v1','HLT_HIJet65_v1','HLT_HIJet80_v1','HLT_HIJet95_v1'],
	useCentrality = False,
	useJEC = cms.untracked.bool(False),
	useVtx = cms.untracked.bool(False),
	isMC = False
	)

process.load('CmsHi.JetAnalysis.rechitanalyzer_cfi')
process.rechitanalyzer.doEcal = cms.untracked.bool(False)
process.rechitanalyzer.doHcal = cms.untracked.bool(False)
process.rechitanalyzer.hasVtx = cms.untracked.bool(False)
process.rechitanalyzer.hcalHFRecHitSrc = "hltHfreco"
process.rechitanalyzer.towersSrc = "hltTowerMakerForAll"
process.rechitanalyzer.JetSrc = "hltHICaloJetCorrected"
process.rechitanalyzer.TowerTreePtMin = cms.untracked.double(4)

process.load("MitHig.PixelTrackletAnalyzer.METAnalyzer_cff")
process.anaMET.METSrc = "hltMet"

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("openhlt_data_l1emul_mcgt.root")
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.ana_step = cms.EndPath(process.HLTHIRecoJetSequenceIC5Corrected*process.patJets*process.hltMet)
process.endjob_step = cms.EndPath(process.endOfProcess*process.l1NtupleProducer*process.hltbitnew*process.icPu5JetAnalyzer*process.rechitanalyzer*process.anaMET)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1simulation_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.ana_step])
process.schedule.extend([process.endjob_step])

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customise_l1EmulatorFromRaw
from L1Trigger.Configuration.customise_l1EmulatorFromRaw import customise 

#call to customisation function customise imported from L1Trigger.Configuration.customise_l1EmulatorFromRaw
process = customise(process)

# set l1 menu from xml
#import L1Trigger.Configuration.L1Trigger_custom
#process = L1Trigger.Configuration.L1Trigger_custom.customiseL1Menu( process )

# customize the HLT to use the emulated results
import HLTrigger.Configuration.customizeHLTforL1Emulator
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimGmtGctGtDigis( process )

# End of customisation functions
