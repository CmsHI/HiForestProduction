#Centrality
cvs co          CondFormats/HIObjects
cvs co          DataFormats/HeavyIonEvent
cvs co          RecoHI/HiCentralityAlgos
cvs co -d       CmsHi/Analysis2010 UserCode/CmsHi/Analysis2010

#ecal and hcal cleaning
#hcal
addpkg RecoLocalCalo/HcalRecAlgos
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/src/HBHETimingShapedFlag.cc
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/interface/HBHETimingShapedFlag.h
cvs co -r V00-00-18 JetMETAnalysis/HcalReflagging
cvs co -d HcalFilter UserCode/MTonjes/HcalFilter

#ecal requires some code from the photon analyzer, 'rm' commands are to truncate the endless chain of packages obviously, this should be cleaned up
cvs co -d CmsHi/PhotonAnalysis UserCode/CmsHi/PhotonAnalysis
rm CmsHi/PhotonAnalysis/plugins/MultiPhotonAnalyzer.*
rm CmsHi/PhotonAnalysis/plugins/SinglePhotonAnalyzer.*
cvs co -d UserCode/HafHistogram UserCode/CmsHi/HafHistogram

# openhlt analyzers
cvs co -r V03-03-08 HLTrigger/HLTanalyzers
cvs co -d CmsHi/HiHLTAlgos UserCode/CmsHi/HiHLTAlgos
# private code: clean up output
cp CmsHi/HiHLTAlgos/modifiedFiles/HLTAnalyzer.h HLTrigger/HLTanalyzers/interface/HLTAnalyzer.h
cp CmsHi/HiHLTAlgos/modifiedFiles/HLTAnalyzer.cc HLTrigger/HLTanalyzers/src/HLTAnalyzer.cc

scram b -j4
