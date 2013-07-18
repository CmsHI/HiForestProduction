#!/bin/sh

cd $CMSSW_BASE/src

# event selection and centrality 
cvs co -r pPbProd_v04 DataFormats/HeavyIonEvent
cvs co -r pPbProd_v10 RecoHI/HiCentralityAlgos
cvs co -r pPbProd_v09 HeavyIonsAnalysis/Configuration

# jets 
cvs co -d CmsHi/JetAnalysis -r HiForest_V02_78 UserCode/CmsHi/JetAnalysis 

# forest
#cvs co -d HiForest -r HiForest_V2_03_01 UserCode/CmsHi/HiForest/V2

scram b
