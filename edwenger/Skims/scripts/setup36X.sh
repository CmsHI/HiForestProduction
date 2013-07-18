#!/bin/sh

cd $CMSSW_BASE/src

# skiming code
cvs co UserCode/edwenger/Skims
cvs co -r 1.4 UserCode/edwenger/Skims/python/patAna_cff.py

# analyzers
cvs co UserCode/edwenger/EvtSelAnalyzer
cvs co UserCode/edwenger/VertexAnalyzer
cvs co UserCode/edwenger/TrackSpectraAnalyzer
cvs co UserCode/edwenger/TrkEffAnalyzer
cvs co -r Y1JAna_0_1_4  -d Saved UserCode/SavedFMa/Saved

# agglomerative vertex
cvs co -r V00-00-07 UserCode/FerencSiklerVertexing

# vertex constraint for refitting
cvs co UserCode/ASYoon/PPTrackingTools/VertexConstraintProducer

# needed for dijetana
cvs co -r CMSSW_3_9_1 RecoHI/HiCentralityAlgos
cvs co -r CMSSW_3_9_1 DataFormats/HeavyIonEvent
cvs co -r CMSSW_3_9_1 HeavyIonsAnalysis/Configuration
cvs co -r CMSSW_3_9_1 RecoHI/HiCentralityAlgos
cvs co -r CMSSW_3_9_1 RecoHI/HiMuonAlgos

# for GEN spectra production
#cvs co UserCode/ASYoon/MCGeneration
#mv UserCode/ASYoon/MCGeneration .
#cvs co Configuration/Generator/python/PythiaUEZ2Settings_cfi.py
#cvs co Configuration/Generator/python/PythiaUECW900ASettings_cfi.py
#cvs co Configuration/Generator/python/PythiaUED6TSettings_cfi.py
#cvs co Configuration/Generator/python/PythiaUEP0Settings_cfi.py
#cvs co Configuration/Generator/python/PythiaUEProQ20Settings_cfi.py
#cvs co Configuration/Generator/python/MinBias_7TeV_pythia8_cff.py


# cleanup
mv UserCode/edwenger .
mv UserCode/ASYoon/PPTrackingTools/ .
rm -r UserCode/ASYoon

# further clean up to keep it light
#rm edwenger/TrackSpectraAnalyzer/macros/ -r
#rm edwenger/TrkEffAnalyzer/macros/ -r

# FW update for 36X
sed -i 's/DEFINE_SEAL_MODULE();//' UserCode/FerencSiklerVertexing/plugins/module.cc
sed -i 's/DEFINE_ANOTHER_FWK_MODULE/DEFINE_FWK_MODULE/' UserCode/FerencSiklerVertexing/plugins/module.cc
