import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Generator_cff import *
from CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff import *

inclusiveJetAnalyzer.eventInfoTag = cms.InputTag("hiSignal")
inclusiveJetAnalyzer.useCentrality   = cms.untracked.bool(False)

icPu5JetAnalyzer = inclusiveJetAnalyzer.clone()

akPu5PFJetAnalyzer = icPu5JetAnalyzer.clone(
    jetTag = 'akPu5PFpatJets',
    genjetTag = 'ak5HiGenJets',
    matchTag = 'akPu5CalopatJets',
    rParam = 0.5
    )

akPu3PFJetAnalyzer = icPu5JetAnalyzer.clone(
    jetTag = 'akPu3PFpatJets',
    genjetTag = 'ak3HiGenJets',
    matchTag = 'akPu3CalopatJets',    
    rParam = 0.3    
    )

akPu5CaloJetAnalyzer = icPu5JetAnalyzer.clone(
    jetTag = 'akPu5patJets',
    genjetTag = 'ak5HiGenJets',
    matchTag = 'akPu5PFpatJets',    
    rParam = 0.5
    
    )

akPu3CaloJetAnalyzer = icPu5JetAnalyzer.clone(
    jetTag = 'akPu3patJets',
    genjetTag = 'ak3HiGenJets',
    matchTag = 'akPu3PFpatJets',    
    rParam = 0.3
    
    )

icPu5JetAnalyzer.jetTag = cms.InputTag("icPu5patJets")

akPu1PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("akPu1PFpatJets"),genjetTag = 'ak1HiGenJets',rParam = 0.1,matchTag = 'akPu1CalopatJets')
akPu2PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("akPu2PFpatJets"),genjetTag = 'ak2HiGenJets',rParam = 0.2,matchTag = 'akPu2CalopatJets')
akPu3PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("akPu3PFpatJets"),genjetTag = 'ak3HiGenJets',rParam = 0.3,matchTag = 'akPu3CalopatJets')
akPu4PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("akPu4PFpatJets"),genjetTag = 'ak4HiGenJets',rParam = 0.4,matchTag = 'akPu4CalopatJets')
akPu5PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("akPu5PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.5,matchTag = 'akPu5CalopatJets')
akPu6PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("akPu6PFpatJets"),genjetTag = 'ak6HiGenJets',rParam = 0.6,matchTag = 'akPu6CalopatJets')

akPu1CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("akPu1CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.1,matchTag = 'akPu1PFpatJets')
akPu2CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("akPu2CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.2,matchTag = 'akPu2PFpatJets')
akPu3CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("akPu3CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.3,matchTag = 'akPu3PFpatJets')
akPu4CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("akPu4CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.4,matchTag = 'akPu4PFpatJets')
akPu5CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("akPu5CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.5,matchTag = 'akPu5PFpatJets')
akPu6CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("akPu6CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.6,matchTag = 'akPu6PFpatJets')

ak1PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak1PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.1,matchTag = 'ak1CalopatJets')
ak2PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak2PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.2,matchTag = 'ak2CalopatJets')
ak3PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak3PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.3,matchTag = 'ak3CalopatJets')
ak4PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak4PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.4,matchTag = 'ak4CalopatJets')
ak5PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak5PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.5,matchTag = 'ak5CalopatJets')
ak6PFJetAnalyzer = akPu3PFJetAnalyzer.clone(jetTag = cms.InputTag("ak6PFpatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.6,matchTag = 'ak6CalopatJets')
ak3CaloJetAnalyzer = akPu3PFJetAnalyzer.clone()

ak1CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak1CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.1,matchTag = 'ak1PFpatJets')
ak2CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak2CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.2,matchTag = 'ak2PFpatJets')
ak3CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak3CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.3,matchTag = 'ak3PFpatJets')
ak4CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak4CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.4,matchTag = 'ak4PFpatJets')
ak5CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak5CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.5,matchTag = 'ak5PFpatJets')
ak6CaloJetAnalyzer = akPu3CaloJetAnalyzer.clone(jetTag = cms.InputTag("ak6CalopatJets"),genjetTag = 'ak5HiGenJets',rParam = 0.6,matchTag = 'ak6PFpatJets')

jetAnalyzers = cms.Sequence(icPu5JetAnalyzer+

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



