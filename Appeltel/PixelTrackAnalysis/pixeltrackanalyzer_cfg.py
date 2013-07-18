import FWCore.ParameterSet.Config as cms

process = cms.Process("PixelTrackAna")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.GlobalTag.globaltag = 'GR10_P_V12::All'

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    #nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_72_1_NF5.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_71_1_u4A.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_70_1_0qM.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_6_1_tKf.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_69_1_WyT.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_68_1_LnF.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_67_1_E99.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_66_1_1IY.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_65_1_7jB.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_64_1_Kvb.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_63_1_UPS.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_62_1_fYv.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_61_1_S6N.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_60_1_TkW.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_5_1_hZU.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_59_1_C6y.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_58_1_RLF.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_57_1_k8H.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_56_1_tSs.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_55_1_hQ3.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_54_1_n0q.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_53_1_OaQ.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_52_1_SOW.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_51_1_StU.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_50_1_8ZM.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_4_1_DIW.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_49_1_usY.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_48_1_DWd.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_47_1_ibn.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_46_1_Rzf.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_45_1_iaF.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_44_1_zY4.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_43_1_Fs0.root',
        '/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_42_1_dHG.root'

    )
)


process.pixelTrackAna = cms.EDAnalyzer('PixelTrackAnalyzer',
    pixelSrc = cms.InputTag('hiLowPtPixelTracks'),
    vertexSrc = cms.InputTag('hiSelectedVertex'),
    genSrc = cms.InputTag('hiGenParticles'),
    doGenerator = cms.bool(False),
    etaCut = cms.double(1.0),
    ptMin = cms.double(0.25)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('pixelTracks.root')
)


process.p = cms.Path(process.pixelTrackAna)
