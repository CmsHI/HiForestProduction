import FWCore.ParameterSet.Config as cms

# reco
from RecoHI.HiEgammaAlgos.HiEgamma_cff import *

# detector responce
from CmsHi.PhotonAnalysis.isoConeInspector_cfi import *
from CmsHi.PhotonAnalysis.ecalHistProducer_cfi import *
from CmsHi.PhotonAnalysis.SpikeInspector_cfi import *

# spike cleaner for island superclsters
cleanPhotons.maxHoverEBarrel = cms.double(100)
  
# clean photon filter
goodPhotons = cms.EDFilter("PhotonSelector",
	src = cms.InputTag("cleanPhotons"),
	cut = cms.string("pt > 18 & hadronicOverEm < 0.2 & abs(eta) < 1.45"),
	)
filterGoodPhotons = cms.EDFilter("PhotonCountFilter",
	src = cms.InputTag("goodPhotons"),
	minNumber = cms.uint32(0), # tmp solution for both gamma and jets
	maxNumber = cms.uint32(999999),
	)

barrelPhotonFilter = cms.Sequence ( goodPhotons * filterGoodPhotons )

photon_extra_reco = cms.Sequence(
	hiEcalClusteringSequence *
	hiPhotonCleaningSequence *
	barrelPhotonFilter
	)
