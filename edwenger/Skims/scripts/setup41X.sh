#!/bin/sh
## This setup is for 2.76 analysis

cd $CMSSW_BASE/src

# skiming code
cvs co UserCode/edwenger/Skims

# analyzers
cvs co UserCode/edwenger/EvtSelAnalyzer
cvs co UserCode/edwenger/VertexAnalyzer
cvs co UserCode/edwenger/TrackSpectraAnalyzer
cvs co UserCode/edwenger/TrkEffAnalyzer

mv UserCode/edwenger .

# jet analyzers
cvs co -d       MNguyen/InclusiveJetAnalyzer UserCode/MNguyen/InclusiveJetAnalyzer

# trigger analyzers
cvs co -r CMSSW_4_1_2_patch2 HLTrigger/HLTanalyzers
# some output cleanup
sed -i 's|std::cout <<EOF
 " Beginning HLTAnalyzer Analysis|//std::cout <<EOF
 " Beginning HLTAnalyzer Analysis|' HLTrigger/HLTanalyzers/src/HLTAnalyzer.cc
sed -i 's|errMax(){return 100;}|errMax(){return 0;}|' HLTrigger/HLTanalyzers/interface/HLTAnalyzer.h

# option to remove TrackerGeomtry dependency
#cvs co UserCode/ASYoon/MISC/python/ExtraVertex_cff.py
#cp UserCode/ASYoon/MISC/python/ExtraVertex_cff.py edwenger/Skims/python/
#cvs co UserCode/ASYoon/MISC/src/TrkEffAnalyzer.cc
#cp UserCode/ASYoon/MISC/src/TrkEffAnalyzer.cc edwenger/TrkEffAnalyzer/src/

# agglomerative vertex
cvs co -r V00-00-07 UserCode/FerencSiklerVertexing
cvs co -r 1.8 UserCode/FerencSiklerVertexing/python  # EDFilter to EDProducer fixed

# vertex constraint for refitting
cvs co UserCode/ASYoon/PPTrackingTools/VertexConstraintProducer

# event selector with specified jet eta-phi
cvs co UserCode/ASYoon/SpectraAna/JetEtaPhiFilter

# jet sorting based on raw jet energy
cvs co RecoJets/JetAlgorithms/interface
cvs co UserCode/ASYoon/MISC/src/JetAlgoHelper.h
cp UserCode/ASYoon/MISC/src/JetAlgoHelper.h RecoJets/JetAlgorithms/interface

# for GEN spectra production
cvs co UserCode/ASYoon/MCGeneration
mv UserCode/ASYoon/MCGeneration .
cvs co Configuration/Generator/python/PythiaUEZ2Settings_cfi.py
cvs co Configuration/Generator/python/PythiaUECW900ASettings_cfi.py
cvs co Configuration/Generator/python/PythiaUED6TSettings_cfi.py
cvs co Configuration/Generator/python/PythiaUEP0Settings_cfi.py
cvs co Configuration/Generator/python/PythiaUEProQ20Settings_cfi.py
cvs co Configuration/Generator/python/MinBias_7TeV_pythia8_cff.py

# cleanup
mv UserCode/ASYoon/* .
rm -r UserCode/ASYoon

# further clean up to keep it light
rm edwenger/TrackSpectraAnalyzer/macros/ -r
rm edwenger/TrkEffAnalyzer/macros/ -r

# needed to re-reconstruct jet excluding calo-cell
cvs co -d Saved/CaloTools UserCode/SavedFMa/Saved/CaloTools

# FW update for 36X
sed -i 's/DEFINE_SEAL_MODULE();//' UserCode/FerencSiklerVertexing/plugins/module.cc
sed -i 's/DEFINE_ANOTHER_FWK_MODULE/DEFINE_FWK_MODULE/' UserCode/FerencSiklerVertexing/plugins/module.cc

# ======== IMPORTANT ===========
# Running PixelTrackProducer require Ferenc's package (not in CVS)
# ~sikler/public/CMSSW_4_1_2_patch1/src/RecoPixelVertexing/PixelLowPtUtilities/

# To convert to XML format
scram build -c
