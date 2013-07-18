import FWCore.ParameterSet.Config as cms

# Build the Objects from AOD (Jets, Muons, Electrons, METs, Taus)
from PhysicsTools.PatAlgos.patSequences_cff import *

# jec
patJetCorrFactors.levels = ['L2Relative','L3Absolute']

# turn off btagging related
patJets.addJetCharge = False
patJets.addBTagInfo = cms.bool(False)
patJets.addDiscriminators = cms.bool(False)
patJets.addGenPartonMatch   = False
patJets.embedGenPartonMatch = False
patJets.genPartonMatch      = ''
patJets.addGenJetMatch      = False
patJets.genJetMatch         = ''
patJets.getJetMCFlavour     = False
patJets.JetPartonMapSource  = ''
patJets.embedCaloTowers			= False
								
# Select jets
selectedPatJets.cut = cms.string('pt > 10')

# clean jets
cleanedPatJets = selectedPatJets.clone(
  src = cms.InputTag("selectedPatJets"),
  cut = cms.string('jetID().n90Hits>1 && jetID().fHPD<0.98 && emEnergyFraction>0.01')
  )

# Define Sequence
patAnaSequence = cms.Sequence(
   patJetCorrections *
	 patJets *
   selectedPatJets
)
