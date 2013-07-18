#############################
#  HiForest2 Configuration  #
#############################

# Simulations
cvs co -r V00-01-14 GeneratorInterface/HiGenCommon
cvs co -r V02-02-26 GeneratorInterface/PyquenInterface

# Track studies
cvs co -r HiForest_V02_00 -d       edwenger/EvtSelAnalyzer                  UserCode/edwenger/EvtSelAnalyzer
cvs co -r HiForest_V02_00 -d       edwenger/HiVertexAnalyzer                UserCode/edwenger/HiVertexAnalyzer
cvs co -r HiForest_V02_01 -d       edwenger/HiTrkEffAnalyzer                UserCode/edwenger/HiTrkEffAnalyzer
cvs co -r HiForest_V02_00 -d       edwenger/TrkEffAnalyzer                  UserCode/edwenger/TrkEffAnalyzer
cvs co -r HiForest_V02_00 -d       edwenger/TrackSpectraAnalyzer            UserCode/edwenger/TrackSpectraAnalyzer
cvs co -r HiForest_V02_00 -d       edwenger/Skims                           UserCode/edwenger/Skims
cvs co -r HiForest_V02_00 -d       edwenger/VertexAnalyzer                  UserCode/edwenger/VertexAnalyzer
cvs co -r HiForest_V02_00 -d       SpectraAna/JetEtaPhiFilter               UserCode/ASYoon/SpectraAna/JetEtaPhiFilter
cvs co -r HiForest_V02_00 -d       SpectraAna/TrackSpectraAnalyzer          UserCode/ASYoon/SpectraAna/TrackSpectraAnalyzer
cvs co -r HiForest_V02_16 -d       MitHig/PixelTrackletAnalyzer             UserCode/MitHig/PixelTrackletAnalyzer
cvs co -r HiForest_V02_00          UserCode/FerencSiklerVertexing

# pat macros
cvs co -r HiForest_V02_00 -d       MNguyen/patMacrosForMC  UserCode/MNguyen/patMacrosForMC
cvs co -r HiForest_V02_00 -d       MNguyen/patMacrosForDataSkims  UserCode/MNguyen/patMacrosForDataSkims
cvs co -r HiForest_V02_00 -d       MNguyen/Configuration  UserCode/MNguyen/Configuration
cvs co -r pPbProd_v07              HeavyIonsAnalysis/Configuration

# pi0 analyzer (optional)
cvs co -r HiForest_V02_00 -d       MNguyen/pi0Analyzer  UserCode/MNguyen/pi0Analyzer

cvs co -r pPbProd_v05              RecoHI/HiCentralityAlgos
cvs co -r pPbProd_v04              DataFormats/HeavyIonEvent

#cvs co -r V_11573_001 -d HcalNoise/DumpHcalNoiseInformationFromRECO UserCode/YiChen/HcalNoise/DumpHcalNoiseInformationFromRECO

# Track validation histograms
cvs co -r HiForest_V02_00 -d       PbPbTrackingTools/HiTrackValidator UserCode/ASYoon/PbPbTrackingTools/HiTrackValidator

# Vertex constraint producer (Yen-Jie where do we need this?)
cvs co -r HiForest_V02_00 -d       PbPbTrackingTools/VertexConstraintProducer 

# Centrality and event plane
#      Yen-Jie: Evt plane configuration need check!
cvs co -r HiForest_V02_01 -d       CmsHi/PhotonAnalysis UserCode/CmsHi/PhotonAnalysis
cvs co -r branch_44x -d       CmsHi/JetAnalysis UserCode/CmsHi/JetAnalysis

cvs co -r backport44X RecoHI/HiJetAlgos

cvs co -r HiForest_V02_00 -d       HcalFilter UserCode/MTonjes/HcalFilter

#PixelTrack
cvs co -r HiForest_V02_00 -d       Appeltel/PixelTrackAnalysis UserCode/Appeltel/PixelTrackAnalysis
#cvs co -r HiForest_V02_01 -d       Appeltel/PixelTracksRun2010 UserCode/Appeltel/PixelTracksRun2010
cvs co -r HiForest_V02_0 -d       Appeltel/PixelTracksRun2010 UserCode/Appeltel/PixelTracksRun2010


#ecal requires some code from the photon analyzer, 'rm' commands are to truncate the endless chain of packages obviously, 

# Top Level Forest cfg, and a consistent version of this script
cvs co -r HEAD -d		   CmsHi/JetAnalysis/test UserCode/CmsHi/JetAnalysis/test
cvs co -r HEAD -d                  CmsHi/JetAnalysis/scripts UserCode/CmsHi/JetAnalysis/scripts

# trigger analyzers
#cvs co -r hi538_01 HLTrigger/HLTanalyzers
cvs co -r hi44X_02 HLTrigger/HLTanalyzers


# event analyzer
cvs co -r  HiForest_V02_02 -d      CmsHi/HiHLTAlgos UserCode/CmsHi/HiHLTAlgos

# muon analyzer
cvs co -r HiForest_V02_05 -d       MuTrig/HLTMuTree UserCode/yenjie/HLTMuTree/

# X analyzer
cvs co -r format_1_11  -d          XbFrame/Xb_frame UserCode/NTUHEP/Xb_frame
# Build!

scram build -c
scram b -j4

