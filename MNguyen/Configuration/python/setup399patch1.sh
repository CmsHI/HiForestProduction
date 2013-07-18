#Particle Flow, only check out if you want to fix muons
#addpkg RecoParticleFlow/PFProducer
#addpkg RecoParticleFlow/PFRootEvent
#addpkg DataFormats/ParticleFlowCandidate

## hiGoodMergedTracks and the track analyzer
cvs co -d       edwenger/HiVertexAnalyzer UserCode/edwenger/HiVertexAnalyzer
cvs co -d       edwenger/HiTrkEffAnalyzer UserCode/edwenger/HiTrkEffAnalyzer
#cvs co -d      Appeltel/PixelTrackAnalysis  UserCode/Appeltel/PixelTrackAnalysis
#cvs co -d      Appeltel/PixelTracksRun2010  UserCode/Appeltel/PixelTracksRun2010
cvs co -d      edwenger/HiTrkEffAnalyzer/python/ UserCode/edwenger/HiTrkEffAnalyzer/python/

# pat macros
cvs co -d      MNguyen/patMacrosForMC  UserCode/MNguyen/patMacrosForMC
cvs co -d      MNguyen/patMacrosForDataSkims  UserCode/MNguyen/patMacrosForDataSkims
cvs co -d      MNguyen/Configuration  UserCode/MNguyen/Configuration
cvs co         HeavyIonsAnalysis/Configuration
#
## Text files no longer needed
##cvs co         CondFormats/JetMETObjects/dat
##cvs co         UserCode/CmsHi/JetAnalysis/data
##cp             UserCode/CmsHi/JetAnalysis/data/*txt CondFormats/JetMETObjects/data/ 
#
#
##cp              /afs/cern.ch/user/m/mnguyen/scratch0/CMSSW_3_9_1_patch1/src/RecoHI/HiJetAlgos/plugins/BuildFile* RecoHI/HiJetAlgos/plugins/ 
#
#
# Jet and HI Software
cvs co -d       MNguyen/InclusiveJetAnalyzer UserCode/MNguyen/InclusiveJetAnalyzer

cvs co -r       V05-01-09 RecoJets/JetProducers
cvs co          RecoHI/HiJetAlgos
addpkg           RecoHI/Configuration

#Centrality
cvs co          CondFormats/HIObjects
cvs co          DataFormats/HeavyIonEvent
cvs co          RecoHI/HiCentralityAlgos
cvs co -d       CmsHi/Analysis2010 UserCode/CmsHi/Analysis2010

# Jet plus tracks
cvs co -r       V00-02-00 RecoJets/JetPlusTracks  

#ecal and hcal cleaning

#hcal
addpkg RecoLocalCalo/HcalRecAlgos
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/src/HBHETimingShapedFlag.cc
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/interface/HBHETimingShapedFlag.h
cvs co -r V00-00-18 JetMETAnalysis/HcalReflagging
cvs co -d HcalFilter UserCode/MTonjes/HcalFilter

#ecal requires some code from the photon analyzer, 'rm' commands are to truncate the endless chain of packages
# obviously, this should be cleaned up
cvs co -d CmsHi/PhotonAnalysis UserCode/CmsHi/PhotonAnalysis
rm CmsHi/PhotonAnalysis/plugins/MultiPhotonAnalyzer.*
rm CmsHi/PhotonAnalysis/plugins/SinglePhotonAnalyzer.*
cvs co -d UserCode/HafHistogram UserCode/CmsHi/HafHistogram

scram b -j4
