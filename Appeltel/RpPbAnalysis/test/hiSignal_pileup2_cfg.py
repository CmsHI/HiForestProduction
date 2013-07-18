import FWCore.ParameterSet.Config as cms

process = cms.Process('HISIGNAL')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.HiEventMixing_cff')
process.load('Configuration.Geometry.GeometryExtendedReco_cff')
process.load('Configuration.Geometry.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.HiGenCommon.VtxSmearedRealisticPPbBoost8TeVCollision_cff')
process.load('SimGeneral.MixingModule.himixGEN_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('SimGeneral.MixingModule.himixSIMExtended_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('SimGeneral.MixingModule.himixDIGI_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_PIon_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.RandomNumberGeneratorService.hiSignal.initialSeed = ivars.randomNumber

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    # input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:../epos_fix_1_1_8RD.root'),
    inputCommands = cms.untracked.vstring('drop *', 
        'keep *_generator_*_*', 
        'keep *_g4SimHits_*_*'),
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
)



process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.6 $'),
    annotation = cms.untracked.string('HIJING Minimum-Bias p-Pb 2 event pileup'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/Appeltel/RpPbAnalysis/test/hiSignal_pileup2_cfg.py,v $')
)

# Additional output definition
#import random
#process.RandomNumberGeneratorService.hiSignal = cms.PSet(
#        initialSeed = cms.untracked.uint32(random.randrange(1,123456789)),
#        engineName = cms.untracked.string('HepJamesRandom')
#)

#different seed
## ----- WARNING ----- Ensure your PYTHONPATH has the same information as sys.path! Else this will fail!
# Without these three lines, the sim vertex will not be smeared correctly
#from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
#randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
#randSvc.populate()

# Output definition

process.RAWDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('HIJING_pileup2.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAWDEBUG')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'STARTHI53_V17::All', '')

process. hiSignal= cms.EDFilter("HijingGeneratorFilter",
    frame = cms.string('CMS     '),
    targ = cms.string('P       '),
    izp = cms.int32(82),
    bMin = cms.double(0),
    izt = cms.int32(1),
    proj = cms.string('A       '),
    comEnergy = cms.double(5020.0),
    iat = cms.int32(1),
    bMax = cms.double(15),
    iap = cms.int32(208),
    rotateEventPlane = cms.bool(True)
)

process.ProductionFilterSequence = cms.Sequence(process.hiSignal)

process.VtxSmearedCommon.src = cms.InputTag("hiSignal")
process.VtxSmeared.src = cms.InputTag("hiSignal")

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen_hi)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWDEBUGoutput_step = cms.EndPath(process.RAWDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWDEBUGoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
        getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)


process.hiGenParticles.srcVector = ["hiSignal","generator"]

# End of customisation functions

