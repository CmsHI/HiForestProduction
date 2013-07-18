{
  // include paths
  //gSystem->AddIncludePath("-I$CMSSW_BASE/src/");

  // helper classes
  gROOT->Macro("CPlot.cc+");

  // styles
  gROOT->Macro("rootlogon_dNdEtaBase.C");
  //gStyle->SetOptStat(1111);
}
