{
  // include paths
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");

  // helper classes
  gROOT->Macro("$CMSSW_BASE/src/edwenger/TrackSpectraAnalyzer/macros/hlt_eff_simple/CPlot.cc+");

  // styles
  gROOT->Macro("$CMSSW_BASE/src/edwenger/TrackSpectraAnalyzer/macros/hlt_eff_simple/rootlogon_dNdEtaBase.C");
  //gStyle->SetOptStat(1111);
}
