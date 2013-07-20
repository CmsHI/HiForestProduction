#include "edwenger/TrkEffAnalyzer/interface/TrkEffHistograms.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"

#include <iostream>
#include <cmath>
using namespace std;

TrkEffHistograms::TrkEffHistograms(const edm::ParameterSet& pset)
{
  fillHistograms         = pset.getParameter<bool>("fillHistograms");
  fillNtuples            = pset.getParameter<bool>("fillNtuples");
  constPtBins            = pset.getParameter<bool>("constPtBins");
  lowPtMode              = pset.getParameter<bool>("lowPtMode");
  mode900GeV             = pset.getParameter<bool>("mode900GeV");
}


TrkEffHistograms::~TrkEffHistograms()
{
}

void 
TrkEffHistograms::declareHistograms()
{
  if(fillNtuples) {
    
    TString leafStr;
    
    trackTrees.push_back(f->make<TTree>("simTrackTree","simTrackTree"));
    leafStr = "ids/I:etas/F:pts/F:hits/I:status/I:acc/I:nrec/I:ptr/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I";
    trackTrees[0]->Branch("simTrackValues", &simTrackValues, leafStr.Data());
    
    trackTrees.push_back(f->make<TTree>("recTrackTree","recTrackTree"));
    leafStr = "charge/I:etar/F:ptr/F:phir/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I:nsim/I:status/I:ids/I:parids/I:etas/F:pts/F";
    trackTrees[1]->Branch("recTrackValues", &recTrackValues, leafStr.Data());
    
  }
  
  if(fillHistograms) {

    // pt bins
    if(!constPtBins){
       const double small = 1e-3;
       double pt;

       // simple rebinning possible with a rebinning facto n = 2, 3, 6 !
       for(pt =   0  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt); // 24 bins
       for(pt =   1.2; pt <   2.4-small; pt +=  0.1 ) ptBins.push_back(pt); // 12 bins
       for(pt =   2.4; pt <   7.2-small; pt +=  0.4 ) ptBins.push_back(pt); // 12 bins
       for(pt =   7.2; pt <  14.4-small; pt +=  1.2 ) ptBins.push_back(pt); // 6 bins
       for(pt =  14.4; pt <  28.8-small; pt +=  2.4 ) ptBins.push_back(pt); // 6 bins 
       for(pt =  28.8; pt <  48.0-small; pt +=  3.2 ) ptBins.push_back(pt); // 6 bins
       for(pt =  48.0; pt <  86.4-small; pt +=  6.4 ) ptBins.push_back(pt); // 6 bins
       for(pt =  86.4; pt < 189.6-small; pt +=  8.6 ) ptBins.push_back(pt); // 6 bins
       ptBins.push_back(189.6);

    }else if(lowPtMode){

       static float ptMin   =  0.0;
       static float ptMax   =  2.0;
       static float ptWidth =  0.02;

       for(double pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
	  ptBins.push_back(pt);

    }else{

       static float ptMin   =  0.0;
       static float ptMax   =  360.0;
       static float ptWidth =  0.5;

       for(double pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
	  ptBins.push_back(pt);

    }
    
    // eta bins
    static float etaMin   = -3.0;
    static float etaMax   =  3.0;
    static float etaWidth =  0.2;

    for(double eta = etaMin; eta < etaMax + etaWidth/2; eta += etaWidth)
      etaBins.push_back(eta);


    // jet et bins
    static float jetMin;
    static float jetMax; // good to be matched with ana 
    static float jetWidth;

    if(!mode900GeV){
       jetMin = 0, jetMax = 2400, jetWidth = 20;
    }else{
       jetMin = 0, jetMax = 300, jetWidth = 5;
    }


    for(double jet = jetMin; jet < jetMax + jetWidth/2; jet += jetWidth)
       jetBins.push_back(jet);


    /*
    double jet;
    const double small = 1e-3;
    for(jet =    0; jet <   10-small; jet +=  10 ) jetBins.push_back(jet);
    for(jet =   10; jet <   1000-small; jet +=  20 ) jetBins.push_back(jet);
    jetBins.push_back(1010);
    */

    // generated
    hgen = f->make<TH2F>("hgen","Gen Tracks;#eta;p_{T} (GeV/c)",
			 etaBins.size()-1, &etaBins[0],
			 ptBins.size()-1, &ptBins[0]);

    // simulated
    hsim = f->make<TH2F>("hsim","Sim Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // accepted
    hacc = f->make<TH2F>("hacc","Accepted Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // efficiency
    heff = f->make<TH2F>("heff","Effic Rec Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // multiply reconstructed
    hmul = f->make<TH2F>("hmul","Mult Rec Tracks;#eta;p_{T} (GeV/c)",
			 etaBins.size()-1, &etaBins[0],
			 ptBins.size()-1, &ptBins[0]);

    // reconstructed
    hrec = f->make<TH2F>("hrec","Rec Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // fakes
    hfak = f->make<TH2F>("hfak","Fake Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // secondary
    hsec = f->make<TH2F>("hsec","Secondary Tracks;#eta;p_{T} (GeV/c)",
			 etaBins.size()-1, &etaBins[0],
			 ptBins.size()-1, &ptBins[0]);
    

    // simulated 3D 
    hsim3D = f->make<TH3F>("hsim3D","Sim Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // efficiency  3D 
    heff3D = f->make<TH3F>("heff3D","Effic Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // multiply reconstructed 3D 
    hmul3D = f->make<TH3F>("hmul3D","Mult Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                           etaBins.size()-1, &etaBins[0],
                           ptBins.size()-1, &ptBins[0],
                           jetBins.size()-1, &jetBins[0]);


    // reconstructed 3D 
    hrec3D = f->make<TH3F>("hrec3D","Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // fakes 3D 
    hfak3D = f->make<TH3F>("hfak3D","Fake Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // secondary
    hsec3D = f->make<TH3F>("hsec3D","Secondary Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
			   etaBins.size()-1, &etaBins[0],
			   ptBins.size()-1, &ptBins[0],
			   jetBins.size()-1, &jetBins[0]);

    // mom resolution (Sim to Rec) 
    hresStoR3D = f->make<TH3F>("hresStoR3D","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               ptBins.size()-1, &ptBins[0]);

    /*
    hresStoR3D_etaS = f->make<TH3F>("hresStoR3D_etaS","Momentum resolution (sim to rec);#jet E_{T} (GeV/c);sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
				    jetBins.size()-1, &jetBins[0],
				    ptBins.size()-1, &ptBins[0],
				    ptBins.size()-1, &ptBins[0]);

    hresStoR3D_etaL = f->make<TH3F>("hresStoR3D_etaL","Momentum resolution (sim to rec);#jet E_{T} (GeV/c);sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
                                    jetBins.size()-1, &jetBins[0],
                                    ptBins.size()-1, &ptBins[0],
                                    ptBins.size()-1, &ptBins[0]);
    */

    // mom resolution (Sim to Rec) v2      
    //hresStoR3D_v2 = f->make<TH3F>("hresStoR3D_v2","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
    //etaBins.size()-1, &etaBins[0],
    //ptBins.size()-1, &ptBins[0],
    //ptBins.size()-1, &ptBins[0]);


  }

}

void
TrkEffHistograms::fillGenHistograms(const GenTrack_t & g)
{
   hgen->Fill(g.etag, g.ptg);
}

void 
TrkEffHistograms::fillSimHistograms(const SimTrack_t & s)
{

  if(fillNtuples && s.status>0){
    simTrackValues = s;
    trackTrees[0]->Fill();
  }

  if(fillHistograms && s.status>0) {
    hsim->Fill(s.etas, s.pts);
    hsim3D->Fill(s.etas, s.pts, s.jetr);
    if(s.acc)    hacc->Fill(s.etas, s.pts);
    if(s.nrec==1) {
       hresStoR3D->Fill(s.etas, s.pts, s.ptr);
       //if(fabs(s.etas)<1.0) hresStoR3D_etaS->Fill(s.jetr, s.pts, s.ptr);
       //if(fabs(s.etas)<2.4) hresStoR3D_etaL->Fill(s.jetr, s.pts, s.ptr);
    }
    if(s.nrec>0) heff->Fill(s.etas, s.pts), heff3D->Fill(s.etas, s.pts, s.jetr);
    if(s.nrec>1) hmul->Fill(s.etas, s.pts), hmul3D->Fill(s.etas, s.pts, s.jetr);
  }

}

void 
TrkEffHistograms::fillRecHistograms(const RecTrack_t & r)
{

  if(fillNtuples){
    recTrackValues = r;
    trackTrees[1]->Fill();
  }

  if(fillHistograms) {
    hrec->Fill(r.etar, r.ptr);
    hrec3D->Fill(r.etar, r.ptr, r.jetr);
    if(!r.nsim) hfak->Fill(r.etar, r.ptr), hfak3D->Fill(r.etar, r.ptr, r.jetr);
    if(r.nsim>0 && r.status<0) hsec->Fill(r.etar, r.ptr), hsec3D->Fill(r.etar, r.ptr, r.jetr); // nsim>0 redudant?
  }

}

void 
TrkEffHistograms::writeHistograms()
{

}
