#include "edwenger/HiTrkEffAnalyzer/interface/HiTrkEffHistograms.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"

#include <iostream>
#include <cmath>
using namespace std;

HiTrkEffHistograms::HiTrkEffHistograms(const edm::ParameterSet& pset)
{
  fillHistograms         = pset.getParameter<bool>("fillHistograms");
  fillNtuples            = pset.getParameter<bool>("fillNtuples");
  ptBinScheme            = pset.getParameter<Int_t>("ptBinScheme");
  lowPtMode              = pset.getParameter<bool>("lowPtMode");
  trkPtMin		 = pset.getParameter<Double_t>("trkPtMin");
  neededCentBins         = pset.getUntrackedParameter<std::vector<Int_t> >("neededCentBins");
}


HiTrkEffHistograms::~HiTrkEffHistograms()
{
}

void 
HiTrkEffHistograms::declareHistograms()
{
  if(fillNtuples) {
    
    TString leafStr;
    
    trackTrees.push_back(f->make<TTree>("simTrackTree","simTrackTree"));
    leafStr = "ids/I:etas/F:pts/F:hits/I:status/I:acc/I:nrec/I:ptr/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I:jetr/F:cbin/I";
    leafStr += ":jetar/F:jrdr/F:jrind/I:jrflavor/I";
    trackTrees[0]->Branch("simTrackValues", &simTrackValues, leafStr.Data());
    
    trackTrees.push_back(f->make<TTree>("recTrackTree","recTrackTree"));
    leafStr = "charge/I:etar/F:ptr/F:phir/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I:nsim/I:status/I:ids/I:parids/I:etas/F:pts/F:jetr/F:cbin/I";
    leafStr += ":jetar/F:jrdr/F:jrind/I:jrflavor/I";
    trackTrees[1]->Branch("recTrackValues", &recTrackValues, leafStr.Data());
    
  }
  
  if(fillHistograms) {

     const Double_t small = 1e-3;

    // pt bins
    if(ptBinScheme==1){
       Double_t pt;

       // simple rebinning possible with a rebinning facto n = 2, 3, 6 !
       for(pt =   0  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt); // 24 bins 
       for(pt =   1.2; pt <   2.4-small; pt +=  0.1 ) ptBins.push_back(pt); // 12 bins
       for(pt =   2.4; pt <   7.2-small; pt +=  0.4 ) ptBins.push_back(pt); // 12 bins
       for(pt =   7.2; pt <  16.8-small; pt +=  1.6 ) ptBins.push_back(pt); // 6 bins 
       for(pt =  16.8; pt <  55.2-small; pt +=  6.4 ) ptBins.push_back(pt); // 6 bins
       for(pt =  55.2; pt < 170.4-small; pt +=  9.6 ) ptBins.push_back(pt); // 12 bins
       ptBins.push_back(170.4);

    } else if(ptBinScheme==2){
       Double_t pt;

       // simple rebinning possible with a rebinning facto n = 2, 3, 6 !
       for(pt =   0  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt); // 24 bins 
       for(pt =   1.2; pt <   2.4-small; pt +=  0.1 ) ptBins.push_back(pt); // 12 bins
       for(pt =   2.4; pt <   7.2-small; pt +=  0.4 ) ptBins.push_back(pt); // 12 bins
       for(pt =   7.2; pt <  20.0-small; pt +=  1.6 ) ptBins.push_back(pt); // 8 bins 
       for(pt =  20.0; pt <  24-small;	 pt +=  2.0 ) ptBins.push_back(pt); // 2 bins
       for(pt =  24.0; pt <  30-small;	 pt +=  3.0 ) ptBins.push_back(pt); // 2 bins
       for(pt =  30.0; pt <  40-small;	 pt +=  5.0 ) ptBins.push_back(pt); // 2 bins
       for(pt =  40.0; pt <  60-small;	 pt += 10.0 ) ptBins.push_back(pt); // 2 bins
       for(pt =  60.0; pt < 300-small;	 pt += 30.0 ) ptBins.push_back(pt); // 8 bins
       ptBins.push_back(400); // total: 72 bins

    } else if(ptBinScheme==3){
      // coarse binning for jet analysis
      const Int_t numPPtBins=24;
      Float_t pptBins[numPPtBins+1] = {0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10,12,15,20,25,30,45,60,90,120,180,300,500};
      vector<Double_t> ptBinsB2(pptBins,pptBins+numPPtBins+1);
      ptBins = ptBinsB2;

    } else if(lowPtMode){

       static Double_t ptMin   =  0.0;
       static Double_t ptMax   =  6.0;
       static Double_t ptWidth =  0.02;

       for(Double_t pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
	  ptBins.push_back(pt);

    }else{

       static Double_t ptMin   =  0.0;
       static Double_t ptMax   =  200.0;
       static Double_t ptWidth =  0.2;

       for(Double_t pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
	  ptBins.push_back(pt);

    }
    
    // eta bins
    static Double_t etaMin   = -2.6;
    static Double_t etaMax   =  2.6;
    static Double_t etaWidth =  0.4;

    for(Double_t eta = etaMin; eta < etaMax + etaWidth/2; eta += etaWidth)
      etaBins.push_back(eta);


    // jet et bins
    /*
    Double_t jet;
    for(jet =    0; jet <   50-small; jet +=  50 ) jetBins.push_back(jet);
    for(jet =   50; jet <   80-small; jet +=  30 ) jetBins.push_back(jet); 
    for(jet =   80; jet < 1000-small; jet +=  25 ) jetBins.push_back(jet);
    jetBins.push_back(1005);
    */

    static Float_t jetMin = 0.0;
    static Float_t jetMax = 1000; // good to be matched with ana 
    static Float_t jetWidth = 20;

    for(Double_t jet = jetMin; jet < jetMax + jetWidth/2; jet += jetWidth)
       jetBins.push_back(jet);

    if(ptBinScheme==3){
      // coarse binning for jet analysis
      const Int_t numJetBins=10;
      Float_t jBins[numJetBins+1] = {0,20,40,60,80,120,160,200,250,500,1000};
      vector<Double_t> jetBinsB2(jBins,jBins+numJetBins+1);
      jetBins = jetBinsB2;
    }

    // simulated
    hsim = f->make<TH2F>("hsim","Sim Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    hsimSPt1D = f->make<TH1F>("hsimSPt1D","Sim Tracks with sim pT;p_{T} (GeV/c)",
			      ptBins.size()-1, &ptBins[0]);
     
    hsimRPt1D = f->make<TH1F>("hsimRPt1D","Sim Tracks with rec pT;p_{T} (GeV/c)",
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

    // eta vs phi
    hetaphi = f->make<TH2F>("hetaphi","Track eta vs phi;#eta;#phi",40,-2.65,2.65, 80,-3.3,3.3);

    for(unsigned i=0;i<neededCentBins.size()-1;i++){
       vhsim3D.push_back(f->make<TH3F>("","Sim Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
				       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], jetBins.size()-1, &jetBins[0]) );
       vheff3D.push_back(f->make<TH3F>("","Effic Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                                       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], jetBins.size()-1, &jetBins[0]) );
       vhmul3D.push_back(f->make<TH3F>("","Mult Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                                       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], jetBins.size()-1, &jetBins[0]) );
       vhrec3D.push_back(f->make<TH3F>("","Rec Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
				      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], jetBins.size()-1, &jetBins[0]) );
       vhfak3D.push_back(f->make<TH3F>("","Fake Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                                      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], jetBins.size()-1, &jetBins[0]) );
       vhsec3D.push_back(f->make<TH3F>("","Secondary Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
				       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], jetBins.size()-1, &jetBins[0]) );
       vhresStoR3D.push_back(f->make<TH3F>("","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
					   etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], ptBins.size()-1, &ptBins[0]) );
       if(i==0) {
	  vhsim3D[i]->SetName(Form("hsim3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
	  vheff3D[i]->SetName(Form("heff3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
	  vhrec3D[i]->SetName(Form("hrec3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
	  vhfak3D[i]->SetName(Form("hfak3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
	  vhmul3D[i]->SetName(Form("hmul3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
	  vhsec3D[i]->SetName(Form("hsec3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
	  vhresStoR3D[i]->SetName(Form("hresStoR3D_cbin%dto%d",neededCentBins[i],neededCentBins[i+1]));
       }else{
	  vhsim3D[i]->SetName(Form("hsim3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
	  vheff3D[i]->SetName(Form("heff3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
          vhrec3D[i]->SetName(Form("hrec3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
          vhfak3D[i]->SetName(Form("hfak3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
	  vhmul3D[i]->SetName(Form("hmul3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
          vhsec3D[i]->SetName(Form("hsec3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
          vhresStoR3D[i]->SetName(Form("hresStoR3D_cbin%dto%d",neededCentBins[i]+1,neededCentBins[i+1]));
       }

    }

  }

}

void 
HiTrkEffHistograms::fillSimHistograms(const SimTrack_t & s)
{

  if(fillNtuples && s.status>0){
    if (trkPtMin>0&&s.pts>=trkPtMin) {
      simTrackValues = s;
      trackTrees[0]->Fill();
    }
  }

  if(fillHistograms && s.status>0) {
    hsim->Fill(s.etas, s.pts);
    hsim3D->Fill(s.etas, s.pts, s.jetr);
    if(s.acc)    hacc->Fill(s.etas, s.pts);
    if(s.nrec==1) {
       if(fabs(s.etas)<1.0) hsimSPt1D->Fill(s.pts), hsimRPt1D->Fill(s.ptr);
       hresStoR3D->Fill(s.etas, s.pts, s.ptr);
    }
    if(s.nrec>0) heff->Fill(s.etas, s.pts), heff3D->Fill(s.etas, s.pts, s.jetr);
    if(s.nrec>1) hmul->Fill(s.etas, s.pts), hmul3D->Fill(s.etas, s.pts, s.jetr);


    // filling histogram in vector 
    for(unsigned i=0;i<neededCentBins.size()-1;i++){
       if(i==0){
          if(s.cbin<=neededCentBins[i+1]){
             vhsim3D[i]->Fill(s.etas, s.pts, s.jetr);
	     if(s.nrec>0) vheff3D[i]->Fill(s.etas, s.pts, s.jetr);
	     if(s.nrec==1) vhresStoR3D[i]->Fill(s.etas, s.pts, s.ptr);
	     if(s.nrec>1) vhmul3D[i]->Fill(s.etas, s.pts, s.jetr);
          }
       }else{
          if(s.cbin>neededCentBins[i] && s.cbin<=neededCentBins[i+1]){
	     vhsim3D[i]->Fill(s.etas, s.pts, s.jetr);
             if(s.nrec>0) vheff3D[i]->Fill(s.etas, s.pts, s.jetr);
	     if(s.nrec==1) vhresStoR3D[i]->Fill(s.etas, s.pts, s.ptr);
             if(s.nrec>1) vhmul3D[i]->Fill(s.etas, s.pts, s.jetr);
          }
       }
    } // end of vector loop 

  } // end of (fillHistogram && s.status) loop 

}

void 
HiTrkEffHistograms::fillRecHistograms(const RecTrack_t & r)
{

  if(fillNtuples){
    if (trkPtMin>0&&r.ptr>=trkPtMin) {
      recTrackValues = r;
      trackTrees[1]->Fill();
    }
  }

  if(fillHistograms) {
    hrec->Fill(r.etar, r.ptr);
    hrec3D->Fill(r.etar, r.ptr, r.jetr);
    if(!r.nsim) hfak->Fill(r.etar, r.ptr), hfak3D->Fill(r.etar, r.ptr, r.jetr);
    if(r.nsim>0 && r.status<0) hsec->Fill(r.etar, r.ptr), hsec3D->Fill(r.etar, r.ptr, r.jetr); // nsim>0 redudant?
    hetaphi->Fill(r.etar,r.phir);

    // filling histogram in vector
    for(unsigned i=0;i<neededCentBins.size()-1;i++){
       if(i==0){
          if(r.cbin<=neededCentBins[i+1]){
	     vhrec3D[i]->Fill(r.etar, r.ptr, r.jetr);
	     if(!r.nsim) vhfak3D[i]->Fill(r.etar, r.ptr, r.jetr);
	     if(r.nsim>0 && r.status<0) vhsec3D[i]->Fill(r.etar, r.ptr, r.jetr);
	  }
       }else{
          if(r.cbin>neededCentBins[i] && r.cbin<=neededCentBins[i+1]){
	     vhrec3D[i]->Fill(r.etar, r.ptr, r.jetr);
	     if(!r.nsim) vhfak3D[i]->Fill(r.etar, r.ptr, r.jetr);
	     if(r.nsim>0 && r.status<0) vhsec3D[i]->Fill(r.etar, r.ptr, r.jetr);
	  }
       }
    } // end of vector loop

  } // end of (fillHistogram) loop

}

void 
HiTrkEffHistograms::writeHistograms()
{

}
