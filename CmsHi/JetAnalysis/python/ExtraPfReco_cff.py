import FWCore.ParameterSet.Config as cms

from RecoHI.Configuration.Reconstruction_hiPF_cff import *

pfTrack.TkColList = cms.VInputTag("hiSelectedTracks")

# PF uisng a grid of pseudo-towers to emulate Russian style subtraction
from RecoHI.HiJetAlgos.ParticleTowerProducer_cff import *
# Needed to produce "HcalSeverityLevelComputerRcd" used by CaloTowersCreator/towerMakerPF
from RecoLocalCalo.Configuration.hcalLocalReco_cff import *
PFTowers = particleTowerProducer.clone()

particleFlowTmp.postMuonCleaning = False
particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)


particleFlowBlock.RecMuons = 'muons'
particleFlowTmp.postMuonCleaning = False
pfTrack.UseQuality = True
pfTrack.GsfTracksInEvents = cms.bool(False)
HiParticleFlowReco.remove(electronsWithPresel)
#HiParticleFlowReco.remove(electronsCiCLoose)

HiParticleFlowReco.remove(mvaElectrons)

particleFlowTmp.usePFElectrons = cms.bool(False)
particleFlowTmp.useEGammaElectrons = cms.bool(False)



HiParticleFlowRecoNoJets = cms.Sequence(
	particleFlowCluster
	* pfTrack
	* particleFlowReco
	* PFTowers
	)



