
#include "SpectraAna/TrackSpectraAnalyzer/interface/TrackSpectraAnalyzerHist.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"

#include <iostream>
#include <cmath>
using namespace std;


TrackSpectraAnalyzerHist::TrackSpectraAnalyzerHist(const edm::ParameterSet& pset)
{

   isGEN_ = pset.getUntrackedParameter<bool>("isGEN");

   string outputFileLabel = pset.getParameter<string>("outputFile");
   outputFile = new TFile(outputFileLabel.c_str(),"recreate");

   etaBinW_ = pset.getUntrackedParameter<double>("etaBinW");
   etaMax_ = pset.getUntrackedParameter<double>("etaMax");

   ptBinW_ = pset.getUntrackedParameter<double>("ptBinW");
   ptMin_ = pset.getUntrackedParameter<double>("ptMin");
   ptMax_ = pset.getUntrackedParameter<double>("ptMax");

   jetaBinW_ = pset.getUntrackedParameter<double>("jetaBinW");
   jetaMax_ = pset.getUntrackedParameter<double>("jetaMax");

   jetBinW_ = pset.getUntrackedParameter<double>("jetBinW");
   jetMin_ = pset.getUntrackedParameter<double>("jetMin");
   jetMax_ = pset.getUntrackedParameter<double>("jetMax");

}

TrackSpectraAnalyzerHist::~TrackSpectraAnalyzerHist()
{
}

void
TrackSpectraAnalyzerHist::prepareHist()
{

   // eta bins

   static double etaMax = etaMax_;
   static double etaMin = -1.0*etaMax_;
   static double etaWidth =  etaBinW_;

   for(float eta = etaMin; eta < etaMax + etaWidth/2; eta += etaWidth)
      etaBins.push_back(eta);

   // pt bins 

   static float ptMin   =  ptMin_;
   static float ptMax   =  ptMax_;
   static float ptWidth =  ptBinW_;

   for(float pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
      ptBins.push_back(pt);

   // jet eta bins 

   static float jetaMax   = jetaMax_;
   static float jetaMin   =  -1.0*jetaMax_;
   static float jetaWidth =  jetaBinW_;

   for(float jeta = jetaMin; jeta < jetaMax + jetaWidth/2; jeta += jetaWidth)
      jetaBins.push_back(jeta);

   // jet et bins

   static float jetMin   =  jetMin_;
   static float jetMax   =  jetMax_;
   static float jetWidth =  jetBinW_;

   for(float jet = jetMin; jet < jetMax + jetWidth/2; jet += jetWidth)
      jetBins.push_back(jet);

   // multiplicity bins                                                                                                                                    
   static float multMin = 0.0;
   static float multMax = 200;
   static float multWidth = 1;

   for(float mult = multMin; mult < multMax + multWidth/2; mult += multWidth)
      multBins.push_back(mult);



   // spectra 

   hTrkPtEta = new TH2F("hTrkPtEta","eta vs pt;#eta;p_{T} (GeV/c)",
                        etaBins.size()-1, &etaBins[0],
                        ptBins.size()-1, &ptBins[0]);

   if(isGEN_) hGenTrkPtEta = new TH2F("hGenTrkPtEta","eta vs pt;#eta;p_{T} (GeV/c)",
				      etaBins.size()-1, &etaBins[0],
				      ptBins.size()-1, &ptBins[0]);

   // jet spectra 2D

   hJetEtEtaMB = new TH2F("hJetEtEtaMB","jeta vs jet;#eta;E_{T} (GeV/c)",
			  jetaBins.size()-1, &jetaBins[0],
			  jetBins.size()-1, &jetBins[0]);

   hJetEtEta6U = new TH2F("hJetEtEta6U","jeta vs jet;#eta;E_{T} (GeV/c)",
                          jetaBins.size()-1, &jetaBins[0],
                          jetBins.size()-1, &jetBins[0]);

   hJetEtEta15U = new TH2F("hJetEtEta15U","jeta vs jet;#eta;E_{T} (GeV/c)",
                          jetaBins.size()-1, &jetaBins[0],
                          jetBins.size()-1, &jetBins[0]);

   hJetEtEta30U = new TH2F("hJetEtEta30U","jeta vs jet;#eta;E_{T} (GeV/c)",
                          jetaBins.size()-1, &jetaBins[0],
                          jetBins.size()-1, &jetBins[0]);

   hJetEtEta50U = new TH2F("hJetEtEta50U","jeta vs jet;#eta;E_{T} (GeV/c)",
			   jetaBins.size()-1, &jetaBins[0],
			   jetBins.size()-1, &jetBins[0]);

			   
   // spectra (with jet et)

   hTrkPtEtaMB = new TH3F("hTrkPtEtaMB","spectra ;#eta;p_{T} (GeV/c);E_{T} (GeV/c)",
			  etaBins.size()-1, &etaBins[0],
                          ptBins.size()-1, &ptBins[0],
			  jetBins.size()-1, &jetBins[0]);


   hTrkPtEta6U = new TH3F("hTrkPtEta6U","spectra ;#eta;p_{T} (GeV/c);E_{T} (GeV/c)",
			  etaBins.size()-1, &etaBins[0],
                          ptBins.size()-1, &ptBins[0],
			  jetBins.size()-1, &jetBins[0]); 


   hTrkPtEta15U = new TH3F("hTrkPtEta15U","spectra ;#eta;p_{T} (GeV/c);E_{T} (GeV/c)",
                          etaBins.size()-1, &etaBins[0],
                          ptBins.size()-1, &ptBins[0],
                          jetBins.size()-1, &jetBins[0]);


   hTrkPtEta30U = new TH3F("hTrkPtEta30U","spectra ;#eta;p_{T} (GeV/c);E_{T} (GeV/c)",
                          etaBins.size()-1, &etaBins[0],
                          ptBins.size()-1, &ptBins[0],
                          jetBins.size()-1, &jetBins[0]);


   hTrkPtEta50U = new TH3F("hTrkPtEta50U","spectra ;#eta;p_{T} (GeV/c);E_{T} (GeV/c)",
                          etaBins.size()-1, &etaBins[0],
                          ptBins.size()-1, &ptBins[0],
                          jetBins.size()-1, &jetBins[0]);


}

void
TrackSpectraAnalyzerHist::fillTrack2DHist(float eta, float pt, int type){
   hTrkPtEta->Fill(eta,pt);
   if(type==1 && isGEN_) hTrkPtEta->Fill(eta,pt);
}

void
TrackSpectraAnalyzerHist::fillJet2DHist(float jeta, float jet, int tbit){

   if(tbit==0)      hJetEtEtaMB->Fill(jeta,jet);
   else if(tbit==1) hJetEtEta6U->Fill(jeta,jet);
   else if(tbit==2) hJetEtEta15U->Fill(jeta,jet);
   else if(tbit==3) hJetEtEta30U->Fill(jeta,jet);
   else if (tbit==4) hJetEtEta50U->Fill(jeta,jet);

}

void
TrackSpectraAnalyzerHist::fillTrackJet3DHist(float eta, float pt, float jet, int tbit){

   if(tbit==0)      hTrkPtEtaMB->Fill(eta,pt,jet);
   else if(tbit==1) hTrkPtEta6U->Fill(eta,pt,jet);
   else if(tbit==2) hTrkPtEta15U->Fill(eta,pt,jet);
   else if(tbit==3) hTrkPtEta30U->Fill(eta,pt,jet);
   else if (tbit==4) hTrkPtEta50U->Fill(eta,pt,jet);

}

void
TrackSpectraAnalyzerHist::writeHist()
{
   outputFile->cd();
   hTrkPtEta->Write();
   if(isGEN_) hGenTrkPtEta->Write();
   hJetEtEtaMB->Write();
   hJetEtEta6U->Write();
   hJetEtEta15U->Write();
   hJetEtEta30U->Write();
   hJetEtEta50U->Write();
   hTrkPtEtaMB->Write();
   hTrkPtEta6U->Write();
   hTrkPtEta15U->Write();
   hTrkPtEta30U->Write();
   hTrkPtEta50U->Write();
   outputFile->Close();
}
