using namespace std;
//#include "utilities.h"
#include "math.h"

#include "hiForest.h"
#include "RpPbTrackingForestAnalyzer_pilot.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>


void RpPbTrackingForestAnalyzer_pilot(char *infName = "/mnt/hadoop/cms/store/user/luck/pA2013_pilot/forest200kHz.root")
{

  // init forest

  HiForest *c = new HiForest(infName,"",cPPb);
  c->doTrackCorrections=0;
  c->InitTree();

  // init output file

  TFile * f = new TFile("trackingAnalysis_pilot202792.root","RECREATE");
  f->mkdir("trkAna_etaFull");
  f->cd("trkAna_etaFull");   

  // init histograms

  std::map<std::string,TH1F*> evtPerf_;
  std::map<std::string,TH2F*> evtPerf2D_;
  std::map<std::string,TH1F*> trkPerf_;
  std::map<std::string,TH2F*> trkPerf2D_;
  std::map<std::string,TH1F*> vtxPerf_;
  std::map<std::string,TH2F*> vtxPerf2D_;

  TH1F* events_;
  TH1F* vertices_;
  TH1F* tracks_;
  TH1F* trackseta_;

  events_ = new TH1F("events","",1,0,1);
  tracks_ = new TH1F("tracks","",1,0,1);
  trackseta_ = new TH1F("trackseta","",1,0,1);
  vertices_ = new TH1F("vertices","",1,0,1);

  evtPerf_["Ntrk"] = new TH1F("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["Nvtx"] = new TH1F("evtNvtx","Primary Vertices per event",10,0,10);

  vtxPerf_["Ntrk"] = new TH1F("vtxNtrk","Tracks per vertex",50,0,200);
  vtxPerf_["x"] = new TH1F("vtxX","Vertex x position",1000,-1,1);
  vtxPerf_["y"] = new TH1F("vtxY","Vertex y position",1000,-1,1);
  vtxPerf_["z"] = new TH1F("vtxZ","Vertex z position",100,-30,30);

  vtxPerf2D_["Ntrk2D"] = new TH2F("vtxNtrk2D","Tracks per vertex;vertex (sorted by Ntrk);Ntrk"
                                ,10,0,10,200,0,200);
  
  trkPerf_["Nhit"] = new TH1F("trkNhit", "Tracks by Number of Valid Hits;N hits",    35,  0,35);
  trkPerf_["pt"] = new TH1F("trkPt", "Track p_{T} Distribution;p_{T} [GeV/c]",100,0,6);
  trkPerf_["eta"] = new TH1F("trkEta", "Track Pseudorapidity Distribution;#eta",50,-2.5,2.5);
  trkPerf_["phi"] = new TH1F("trkPhi", "Track Azimuthal Distribution;#phi",100,-3.15,3.15);
  trkPerf_["chi2"] = new TH1F("trkChi2", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",60,0,6);
  trkPerf_["pterr"] = new TH1F("trkPterr", "Track p_{T} error;#delta p_{T} / p_{T}",50,0,0.2);
  trkPerf_["dxyErr"] = new TH1F("trkDxyErr", "Transverse DCA Significance;dxy / #sigma_{dxy}",100,-8,8);
  trkPerf_["dzErr"] = new TH1F("trkDzErr", "Longitudinal DCA Significance;dz / #sigma_{dz}",100,-8,8);

  trkPerf2D_["etaphi"] = new TH2F("trkEtaPhi","Track Eta-Phi Map;#eta;#phi",50,-2.5,2.5,100,-3.15,3.15);

  // event loop

//   for (int i=0;i<c->GetEntries();i++) {
  for (int i=0;i<1000;i++) 
  {
    c->GetEntry(i);
    if (i % 2000 == 0) cout <<i<<" / "<<c->GetEntries()<<endl;
 
    //event selection
    bool event_accepted = true;
    if(!(c->skim.phfPosFilter1 && c->skim.phfNegFilter1
        && c->skim.phltPixelClusterShapeFilter
        && c->skim.pprimaryvertexFilter
        && c->hlt.HLT_PAZeroBiasPixel_SingleTrack_v1)) event_accepted = false;
    //if(fabs(c->evt.vz)>15.) event_accepted = false;
    if(c->evt.run!=202792) event_accepted = false;
    if(c->evt.lumi<317 || c->evt.lumi>1014) event_accepted = false;

    if(event_accepted == false)
      continue;
    
    events_->Fill(0.5);
    evtPerf_["Ntrk"]->Fill(c->track.nTrk);

    // Vertices
    if( c->track.vz[1] < -90. ) continue; 
    //vtxPerf_["Ntrk"]->Fill(vi.tracksSize());
     vtxPerf_["x"]->Fill(c->track.vx[1]);
     vtxPerf_["y"]->Fill(c->track.vy[1]);
     vtxPerf_["z"]->Fill(c->track.vz[1]);
        //vtxPerf2D_["Ntrk2D"]->Fill(vcount,vi.tracksSize());
  //  vertices_->Fill(0.5);        

//    evtPerf_["Nvtx"]->Fill(nVtx);

    // Tracks
    for(int j=0;j<c->track.nTrk;j++) 
    {
      trkPerf_["eta"]->Fill(c->track.trkEta[j] );
      trkPerf2D_["etaphi"]->Fill( c->track.trkEta[j], c->track.trkPhi[j] );
      tracks_->Fill(0.5);
      trackseta_->Fill(0.5);
      trkPerf_["Nhit"]->Fill(c->track.trkNHit[j]); 
      trkPerf_["pt"]->Fill(c->track.trkPt[j]); 
      trkPerf_["phi"]->Fill(c->track.trkPhi[j]); 
      trkPerf_["dxyErr"]->Fill(c->track.trkDxy1[j]/c->track.trkDxyError1[j]); 
      trkPerf_["dzErr"]->Fill(c->track.trkDz1[j]/c->track.trkDzError1[j]); 
      trkPerf_["chi2"]->Fill(c->track.trkChi2[j]/c->track.trkNdof[j]);
      trkPerf_["pterr"]->Fill(c->track.trkPtError[j]/c->track.trkPt[j] );
    }
  }

  f->Write();
  f->Close();

}

