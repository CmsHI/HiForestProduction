#!/bin/bash -
tag=V0531_v4

# use histo
#inFile=$scratch/data/MinimumBias/MB-C10-PR9-TRKANASKIM-v2_job3/all/trkhists_histOnly_all.root
#inFile=$scratch/data/MinimumBias/MB-C10-PR9-TRKANASKIM-v3/all/trkhists_histOnly_all.root
#root -b -q anaHltJetEff_simple.C+\(\"$inFile\",\"plots/${tag}_fromHist\",true\)

# use nt
#inFile=$scratch/data/MinimumBias/MB-C10-PR9-TRKANASKIM-v2_job3/trkhists_trkAnaSkimAOD_*_1.root
inFile=$scratch/data/MinimumBias/MB-C10-PR9-TRKANASKIM-v3_proc0531/trkhists_*.root
root -b -q anaHltJetEff_simple.C+\(\"$inFile\",\"plots/${tag}_fromNt\",false\)

# use hist from nt
inFile=plots/${tag}_fromNt/anahlt.root
root -b -q anaHltJetEff_simple.C+\(\"$inFile\",\"plots/${tag}_histFromNt\",true,\"\",1\)
