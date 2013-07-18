#/bin/bash -
# l1
cvs co UserCode/L1TriggerDPG
cvs co RecoLuminosity/LumiDB
addpkg L1Trigger/Configuration

# HLT
cvs co -r V12-02-45 HLTrigger/Configuration
cvs co -r V03-08-14 HLTrigger/HLTanalyzers

# get hlt menu
hltGetConfiguration --cff --offline --data /users/frankma/devCMSSW_4_2_0/HIonV976 --type HIon --unprescale > $CMSSW_BASE/src/HLTrigger/Configuration/python/HLT_HIon_data_cff.py

# hiGoodMergedTracks and the track analyzer
cvs co -d       edwenger/HiVertexAnalyzer UserCode/edwenger/HiVertexAnalyzer
cvs co -d       edwenger/HiTrkEffAnalyzer UserCode/edwenger/HiTrkEffAnalyzer

# jet analyzers
cvs co          RecoHI/HiJetAlgos
cvs co         HeavyIonsAnalysis/Configuration
cvs co -d       CmsHi/JetAnalysis UserCode/CmsHi/JetAnalysis
rm CmsHi/JetAnalysis/python/EGammaAnalyzers_cff.py # tmp fix
rm CmsHi/JetAnalysis/python/ExtraEGammaReco_cff.py # tmp fix
rm CmsHi/JetAnalysis/python/EventSelection_cff.py # tmp fix

# trigger analyzers
cvs co -d       CmsHi/HiHLTAlgos UserCode/CmsHi/HiHLTAlgos

# custom l1 menu from xml
addpkg L1TriggerConfig/L1GtConfigProducers
cvs co L1TriggerConfig/L1GtConfigProducers/data/Luminosity/startup/L1Menu_CollisionsHeavyIons2011_v0_L1T_Scales_20101224_Imp0_0x1026.xml
cp -v CmsHi/HiHLTAlgos/modifiedFiles/L1Trigger_custom.py L1Trigger/Configuration/python/

scram build -c
