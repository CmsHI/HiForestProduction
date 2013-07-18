import FWCore.ParameterSet.Config as cms

from GeneratorInterface.Pythia6Interface.pythiaDefault_cff import *
from Configuration.Generator.PyquenDefaultSettings_cff import *

hiSignal = cms.EDFilter("Pythia6GeneratorFilter",
                        maxEventsToPrint = cms.untracked.int32(5),
                        pythiaPylistVerbosity = cms.untracked.int32(1),
                        filterEfficiency = cms.untracked.double(1.0),
                        pythiaHepMCVerbosity = cms.untracked.bool(False),
                        comEnergy = cms.double(2760.0),
                        PythiaParameters = cms.PSet(
    pythiaUESettingsBlock,
    processParameters = cms.vstring(
    'MSEL=11             !Z0/gamma* production, ISUB=1',
    'MSTP(43)=2          !Only Z0',
    'MDME(174,1)=0       !Z decay into d dbar',
    'MDME(175,1)=0       !Z decay into u ubar',
    'MDME(176,1)=0       !Z decay into s sbar',
    'MDME(177,1)=0       !Z decay into c cbar',
    'MDME(178,1)=0       !Z decay into b bbar',
    'MDME(179,1)=0       !Z decay into t tbar',
    'MDME(182,1)=1       !Z decay into e- e+',
    'MDME(183,1)=0       !Z decay into nu_e nu_ebar',
    'MDME(184,1)=0       !Z decay into mu- mu+',
    'MDME(185,1)=0       !Z decay into nu_mu nu_mubar',
    'MDME(186,1)=0       !Z decay into tau- tau+',
    'MDME(187,1)=0       !Z decay into nu_tau nu_taubar',
    'CKIN(4) = -1.       !(D=-1 GeV) upper lim pT_hat, if < 0 innactive',
    'CKIN(7) = -3        !(D=-10) lower limit rapidity',
    'CKIN(8) = 3         !(D=10) upper limit rapidity'
    ),
    parameterSets = cms.vstring('pythiaUESettings',
                                'processParameters')
    )
                        )
