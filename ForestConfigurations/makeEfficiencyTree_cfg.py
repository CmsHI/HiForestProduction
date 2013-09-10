#flake8: noqa
'''

Generate trees for measuring and comparing L1 and UCT efficiencies with
respect to RECO objects.

Usage:

    ./makeEfficiencyTree_cfg.py

Optional arguments:

    inputFiles=myFile.root outputFile=outputFile.root maxEvents=-1

Authors: L. Dodd, N. Woods, I. Ojalvo, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

# Get command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# Set useful defaults
#options.inputFiles = 'file:step2_RAW2DIGI_L1Reco_RECO.root'
#options.inputFiles = 'file:/mnt/hadoop/cms/store/user/icali/HIHighPt/Skim_HLT_HIJet55_538HIpa2_RECO/a37cd2798e1f092419699a12291f1745/step2_RAW2DIGI_L1Reco_RECO_458_2_oQA.root'
options.outputFile = "uct_efficiency_HIsub.root"

options.register(
    'eicIsolationThreshold',
    3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "EIC Isolation threshold")
options.register(
    'ecalCalib',
    'CALIB_V4',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'Can be CALIB_V1, CALIB_V3, or CALIB_V4')
options.register(
    'eicCardHcalOnly',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'If 1, turn off the ECAL for the stage1 EGTau path.')
options.register(
    'isMC',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.')


options.parseArguments()

process = cms.Process("L1UCTEfficiency")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR_P_V43D::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)

# Load emulation and RECO sequences
if not options.isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")

# Determine which calibration to use
from L1Trigger.UCT2015.emulation_cfi import \
        eg_calib_v1, eg_calib_v3, eg_calib_v4

calib_map = {
    'CALIB_V1': eg_calib_v1,
    'CALIB_V3': eg_calib_v3,
    'CALIB_V4': eg_calib_v4
}

ecal_calibration = calib_map[options.ecalCalib]
process.RCTConfigProducers.eGammaECalScaleFactors = ecal_calibration
process.RCTConfigProducers.jetMETECalScaleFactors = ecal_calibration
process.UCT2015EClusterProducer.ecalCalibration = ecal_calibration

if options.eicCardHcalOnly:
    print "Disabling ECAL in Stage1 EGTau path"
    # Set all scale factors to 0.
    process.RCTConfigProducers.eGammaECalScaleFactors = [
        0.0 for _ in eg_calib_v1]

# Common branches to add to the ntuple
common_ntuple_branches = cms.PSet(
    index = cms.string("index"), # Index of reco object in the event
    nRecoObjects = cms.string("nTotalObjects"), # Number of reco objects in the event
    nPVs = cms.string("nPVs"), # number of reco'ed vertices in the event

    # Run, lumi, event number
    run = cms.string("id.run"),
    lumi = cms.string("id.luminosityBlock"),
    evt = cms.string("id.event"),

    recoPt = cms.string("reco.pt"),
    recoEta = cms.string("reco.eta"),
    recoPhi = cms.string("reco.phi"),

    # Whether there exists a L1/UCT object matched to reco
    l1Match = cms.string("l1Match"),
    l1gMatch = cms.string("l1gMatch"),

    l1Pt = cms.string("? l1Match ? l1.pt : 0"),
    l1Eta = cms.string("? l1Match ? l1.eta : 0"),
    l1Phi = cms.string("? l1Match ? l1.phi : 0"),
    l1Type = cms.string("? l1Match ? l1.type() : -1"),


    # TODO add L1extra eta/phi indices

    l1DPhi = cms.string("? l1Match ? deltaPhi(l1.phi, reco.phi) : -1"),
    l1DR = cms.string("? l1Match ? deltaR(l1.eta, l1.phi, reco.eta, reco.phi) : -1"),

    l1gPt = cms.string("? l1gMatch ? l1g.pt : 0"),
    l1gEta = cms.string("? l1gMatch ? l1g.eta : 0"),
    l1gPhi = cms.string("? l1gMatch ? l1g.phi : 0"),

    # For tuning isolation and PU subtraction
    l1gPU = cms.string("? l1gMatch ? l1g.getFloat('puLevel', -4) : -2"),
    l1gPUUIC = cms.string("? l1gMatch ? l1g.getFloat('puLevelUIC', -4) : -2"),
    l1gRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedRegionEt', -4) : -2"),

    l1gEtaCode = cms.vstring("? l1gMatch ? l1g.getInt('rgnEta') : 0", "I"),
    l1gPhiCode = cms.vstring("? l1gMatch ? l1g.getInt('rgnPhi') : 0", "I"),

    l1gDPhi = cms.string("? l1gMatch ? deltaPhi(l1g.phi, reco.phi) : -1"),
    l1gDEta = cms.string("? l1gMatch ? l1g.eta - reco.eta : -10"),
    l1gDR = cms.string("? l1gMatch ? deltaR(l1g.eta, l1g.phi, reco.eta, reco.phi) : -1"),
)

process.jetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    #recoSrc = cms.VInputTag("recoJets"),
    #recoSrc = cms.VInputTag( cms.InputTag("l1extraParticles", "Central") ),
    recoSrc = cms.VInputTag("iterativeConePu5CaloJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "JetUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    ),
    useVertex = cms.bool(False)
)

# process.corrjetEfficiency = cms.EDAnalyzer(
#     "EfficiencyTree",
#     recoSrc = cms.VInputTag("iterativeConePu5CaloJets:"),
#     l1Src = cms.VInputTag(
#         # Combine central jets + tau + forward jets
#         cms.InputTag("l1extraParticles", "Central"),
#         cms.InputTag("l1extraParticles", "Tau"),
#         cms.InputTag("l1extraParticles", "Forward"),
#     ),
#     l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "CorrJetUnpacked")),
#     l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
#     # Max DR for RECO-trigger matching
#     maxDR = cms.double(0.5),
#     # Ntuple configuration
#     ntuple = cms.PSet(
#         common_ntuple_branches,
#     ),
#     pvSrc = cms.InputTag("hiSelectedVertex")
# )

process.p1 = cms.Path(
#    process.emulationSequence
    process.uctEmulatorStep *
    process.jetEfficiency 
    #process.corrjetEfficiency
)

# Use the HI background subtraction. 
process.UCT2015Producer.useUICrho = cms.bool(False)
process.UCT2015Producer.useHI = cms.bool(True)

# Make a version of UCT without PU corrections.
# process.UCT2015ProducerNoPU = process.UCT2015Producer.clone(
#     puCorrect = False
# )

process.schedule = cms.Schedule(
    process.p1,
)

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

eic = options.eicIsolationThreshold
print "Setting EIC threshold to %i" % eic
process.RCTConfigProducers.eicIsolationThreshold = eic
