import FWCore.ParameterSet.Config as cms


from CmsHi.JetAnalysis.JetAnalyzers_cff import *

icPu5JetAnalyzer.fillGenJets   = cms.untracked.bool(True)
akPu1PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu2PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu3PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu4PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu5PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu6PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)

ak1PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak2PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak3PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak4PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak5PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak6PFJetAnalyzer.fillGenJets = cms.untracked.bool(True)

akPu1CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu2CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu3CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu4CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu5CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
akPu6CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)

ak1CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak2CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak3CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak4CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak5CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)
ak6CaloJetAnalyzer.fillGenJets = cms.untracked.bool(True)

icPu5JetAnalyzer.genPtMin = cms.untracked.double(15)
akPu1PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akPu2PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akPu3PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akPu4PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akPu5PFJetAnalyzer.genPtMin = cms.untracked.double(15)
akPu6PFJetAnalyzer.genPtMin = cms.untracked.double(15)

icPu5JetAnalyzer.isMC   = cms.untracked.bool(True)
akPu1PFJetAnalyzer.isMC = cms.untracked.bool(True)
akPu2PFJetAnalyzer.isMC = cms.untracked.bool(True)
akPu3PFJetAnalyzer.isMC = cms.untracked.bool(True)
akPu4PFJetAnalyzer.isMC = cms.untracked.bool(True)
akPu5PFJetAnalyzer.isMC = cms.untracked.bool(True)
akPu6PFJetAnalyzer.isMC = cms.untracked.bool(True)
ak1PFJetAnalyzer.isMC = cms.untracked.bool(True)
ak2PFJetAnalyzer.isMC = cms.untracked.bool(True)
ak3PFJetAnalyzer.isMC = cms.untracked.bool(True)
ak4PFJetAnalyzer.isMC = cms.untracked.bool(True)
ak5PFJetAnalyzer.isMC = cms.untracked.bool(True)
ak6PFJetAnalyzer.isMC = cms.untracked.bool(True)
akPu1CaloJetAnalyzer.isMC = cms.untracked.bool(True)
akPu2CaloJetAnalyzer.isMC = cms.untracked.bool(True)
akPu3CaloJetAnalyzer.isMC = cms.untracked.bool(True)
akPu4CaloJetAnalyzer.isMC = cms.untracked.bool(True)
akPu5CaloJetAnalyzer.isMC = cms.untracked.bool(True)
akPu6CaloJetAnalyzer.isMC = cms.untracked.bool(True)
ak1CaloJetAnalyzer.isMC = cms.untracked.bool(True)
ak2CaloJetAnalyzer.isMC = cms.untracked.bool(True)
ak3CaloJetAnalyzer.isMC = cms.untracked.bool(True)
ak4CaloJetAnalyzer.isMC = cms.untracked.bool(True)
ak5CaloJetAnalyzer.isMC = cms.untracked.bool(True)
ak6CaloJetAnalyzer.isMC = cms.untracked.bool(True)




jetAnalyzers = cms.Sequence(
    
                            akPu1PFJetAnalyzer +
                            akPu2PFJetAnalyzer +
                            akPu3PFJetAnalyzer +
                            akPu4PFJetAnalyzer +
                            akPu5PFJetAnalyzer +
                            akPu6PFJetAnalyzer +
                            
                            akPu1CaloJetAnalyzer +
                            akPu2CaloJetAnalyzer +
                            akPu3CaloJetAnalyzer +
                            akPu4CaloJetAnalyzer +
                            akPu5CaloJetAnalyzer +
                            akPu6CaloJetAnalyzer +
                            
                            ak1PFJetAnalyzer +
                            ak2PFJetAnalyzer +
                            ak3PFJetAnalyzer +
                            ak4PFJetAnalyzer +
                            ak5PFJetAnalyzer +
                            ak6PFJetAnalyzer +
                            
                            ak1CaloJetAnalyzer +
                            ak2CaloJetAnalyzer +
                            ak3CaloJetAnalyzer +
                            ak4CaloJetAnalyzer +
                            ak5CaloJetAnalyzer +
                            ak6CaloJetAnalyzer
                            
                            )



jetAnalyzers2to5 = cms.Sequence(
    
                            akPu2PFJetAnalyzer +
                            akPu3PFJetAnalyzer +
                            akPu4PFJetAnalyzer +
                            akPu5PFJetAnalyzer +
                            
                            akPu2CaloJetAnalyzer +
                            akPu3CaloJetAnalyzer +
                            akPu4CaloJetAnalyzer +
                            akPu5CaloJetAnalyzer +

                            ak2PFJetAnalyzer +
                            ak3PFJetAnalyzer +
                            ak4PFJetAnalyzer +
                            ak5PFJetAnalyzer +
                            
                            ak2CaloJetAnalyzer +
                            ak3CaloJetAnalyzer +
                            ak4CaloJetAnalyzer +
                            ak5CaloJetAnalyzer
                            
                            )



