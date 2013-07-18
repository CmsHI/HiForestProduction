import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')

ivars.files = 'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_7_0_pre2/RelValPyquen_DiJet_pt80to120_2760GeV/GEN-SIM-RECO/MC_37Y_V1-v1/0018/1420F252-8453-DF11-9903-002618943882.root'
ivars.output = 'hi_pat_80to120.root'
ivars.maxEvents = 3

ivars.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(ivars.files)
                            )

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(ivars.maxEvents)
                    )

process.load('PhysicsTools.PatAlgos.patHeavyIonSequences_cff')
from PhysicsTools.PatAlgos.tools.heavyIonTools import *
configureHeavyIons(process)

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Summer09_7TeV_ReReco332")

process.load("HeavyIonsAnalysis.Configuration.analysisProducers_cff")
process.load("HeavyIonsAnalysis.Configuration.analysisEventContent_cff")

process.output = cms.OutputModule("PoolOutputModule",
                                  process.jetTrkSkimContent,
                                  fileName = cms.untracked.string(ivars.output)
                                  )

process.prod = cms.Path(process.allTracks +
                        process.heavyIon +
                        process.makeHeavyIonJets
                        )

process.out_step = cms.EndPath(process.output)

