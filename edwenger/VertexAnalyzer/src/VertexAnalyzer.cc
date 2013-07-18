// Original Author:  Edward Allen Wenger,32 4-A06,+41227676980,
//         Created:  Fri May  7 13:11:39 CEST 2010
// $Id: VertexAnalyzer.cc,v 1.20 2011/02/08 16:18:39 sungho Exp $
//

#include "edwenger/VertexAnalyzer/interface/VertexAnalyzer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "edwenger/VertexAnalyzer/interface/VertexComparator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& iConfig)
:
  vtxlabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxlabel")),
  trklabel_(iConfig.getUntrackedParameter<edm::InputTag>("trklabel")),
  jetlabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetlabel")),
  dzcut_(iConfig.getUntrackedParameter<double>("dzcut")),
  minJetEt_(iConfig.getUntrackedParameter<double>("minJetEt")),
  nthleadingMax_(iConfig.getUntrackedParameter<int>("nthleadingMax")),
  nTracksBins_(iConfig.getUntrackedParameter<int>("nTracksBins")),
  jetTrkVerticesCorr_(iConfig.getUntrackedParameter<bool>("jetTrkVerticesCorr"))

{

}


// ------------ method called to for each event  ------------
void
VertexAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace std;
   using namespace reco;
   
   //------- Vertices ---------------------
   edm::Handle<reco::VertexCollection> vtxsH;
   iEvent.getByLabel(vtxlabel_,vtxsH);

   reco::VertexCollection vtxs = *vtxsH;
   if(jetTrkVerticesCorr_) std::sort(vtxs.begin(),vtxs.end(),VertexHigherPtSquared()); // sorted by sum pt^2 for jet-vtx correlation
   else std::sort(vtxs.begin(),vtxs.end(),MoreTracksThenLowerChi2<reco::Vertex>());

   //------- Vertex analyzer ---------------
   hVtxSize->Fill(vtxs.size());

   if(!vtxs.size()) return;  // this protects from an event with an empty vertex

   hVtxTracks->Fill(vtxs[0].tracksSize());
   hVtxZ->Fill(vtxs[0].z());

   for(unsigned i=0; i<vtxs.size(); i++)
      hAllVtxZ->Fill(vtxs[i].z());

   if(vtxs.size()<2) return; // downstream below is gauranteed to deal with at least 2 vertices

   hMultVtxTracks->Fill(vtxs[0].tracksSize(),vtxs[1].tracksSize());
   hMultVtxZ->Fill(vtxs[0].z(),vtxs[1].z());
      
   for(unsigned i=0; i<vtxs.size(); i++)
      edm::LogVerbatim("VertexAnalyzer") << "vtx #" << i << " has "
                                         << vtxs[i].tracksSize() << " tracks";

   
   // -------- Vertex,Track,Jet analyzer --------
   if(jetTrkVerticesCorr_){  // this is to check jet,track,vertices correlation

      //------ Jets ---------------------------- 
      edm::Handle<reco::CandidateView> jets;
      iEvent.getByLabel(jetlabel_,jets);

      std::vector<const reco::Candidate *> sortedJets;

      for(unsigned it=0; it<jets->size(); ++it){
         const reco::Candidate* jet = &((*jets)[it]);
         sortedJets.push_back(jet);
         sortByEtRef (&sortedJets);
      }

      // Get Leading jet energy 
      double leadJetEt_=0.0, leadJetEta_=0.0, leadJetPhi_=0.0;
      double sleadJetEt_=0.0, sleadJetEta_=0.0, sleadJetPhi_=0.0;

      unsigned index = 0;

      if(sortedJets.size()==0){
	 leadJetEt_ = 0,leadJetEta_ = 100, leadJetPhi_ = 0; // so that dr > 10 for events with no pat jet
      }else{
	 leadJetEt_ = sortedJets[index]->et(), leadJetEta_ = sortedJets[index]->eta(), leadJetPhi_ = sortedJets[index]->phi();
	 if(sortedJets.size()>1) sleadJetEt_ = sortedJets[1]->et(), sleadJetEta_ = sortedJets[1]->eta(), sleadJetPhi_ = sortedJets[1]->phi();
      }
	 

      if(minJetEt_>leadJetEt_) return;  // if jet 
 
      if(sortedJets.size()>1) 
	 hLeadnSLeadJetEta->Fill(leadJetEta_,sortedJets[1]->eta());
      else 
	 hLeadMonoJetEta->Fill(leadJetEta_); // events with mono-jet

      //------ Tracks -------------------------- 
      edm::Handle<std::vector<reco::Track> > tracks;
      iEvent.getByLabel(trklabel_, tracks);

      const std::string qualityString = "highPurity";

      std::vector<const reco::Track *> sortedTrks;

      for(unsigned it=0; it<tracks->size(); ++it){
         const reco::Track & trk = (*tracks)[it];

         // track selections - do we need eta cut? 
         if(!trk.quality(reco::TrackBase::qualityByName(qualityString))) continue;

         sortedTrks.push_back(&(*tracks)[it]);
         sortByPtRef (&sortedTrks);

         hVtxZTrkEta->Fill(vtxs[0].z(),trk.eta());
      }

      // Tracks/Vertices corrleation ========================
      double dzPV_leading=0, dzSV_leading=0;
      double dzPV_sleading=0, dzSV_sleading=0;

      if(sortedTrks.size()>0){
	 dzPV_leading = sortedTrks[0]->dz(vtxs[0].position());
	 dzSV_leading = sortedTrks[0]->dz(vtxs[1].position());
      }

      if(sortedTrks.size()>1){
	 dzPV_sleading = sortedTrks[1]->dz(vtxs[0].position());
	 dzSV_sleading = sortedTrks[1]->dz(vtxs[1].position());
      }

      unsigned int nthleading = 0.0;
      double zPV=vtxs[0].z(), zSV=vtxs[1].z();
      double trkSizePV=vtxs[0].tracksSize(), trkSizeSV=vtxs[1].tracksSize();
      double dzPV_mostCorr=0, dzSV_mostCorr=0;
      double ptSum_PV=0, ptSum_SV=0;
      double least_dr_corr = 999.;

      for(unsigned j=0; j<sortedTrks.size(); j++){
	 
	 double dzPV = sortedTrks[j]->dz(vtxs[0].position());
	 double dzSV = sortedTrks[j]->dz(vtxs[1].position());
	 
	 // correlation with leading jet
	 double deta = leadJetEta_ - sortedTrks[j]->eta();
	 double dphi = leadJetPhi_ - sortedTrks[j]->phi();

	 if(fabs(dphi)>(TMath::Pi()) && dphi>0) dphi = dphi - 2.*TMath::Pi();
         if(fabs(dphi)>(TMath::Pi()) && dphi<0) dphi = dphi + 2.*TMath::Pi();

	 double dr = TMath::Sqrt(deta*deta + dphi*dphi); // dr < pi when deta = 0

	 // correlation with sub-leading jet
	 double deta_sub = sleadJetEta_ - sortedTrks[j]->eta();
         double dphi_sub = sleadJetPhi_ - sortedTrks[j]->phi();

         if(fabs(dphi_sub)>(TMath::Pi()) && dphi_sub>0) dphi_sub = dphi_sub - 2.*TMath::Pi();
         if(fabs(dphi_sub)>(TMath::Pi()) && dphi_sub<0) dphi_sub = dphi_sub + 2.*TMath::Pi();

         double dr_sub = TMath::Sqrt(deta_sub*deta_sub + dphi_sub*dphi_sub); 
	 
	 if(sortedJets.size()<2) dr_sub = 999.9; // if no second jet, assign large value

	 double dr_corr = min(dr,dr_sub);

	 if(j<nthleadingMax_){
            if(dr_corr<least_dr_corr)
	       least_dr_corr = dr_corr, nthleading = j, dzPV_mostCorr = dzPV, dzSV_mostCorr = dzSV;
	 }

	 if(j==0) {
	    dEtadPhi_leading->Fill(deta,dphi);
	    hLeadingTrkPVdZ->Fill(dzPV), hLeadingTrkSVdZ->Fill(dzSV), hLeadingPVdZdR->Fill(dzPV,dr), hLeadingSVdZdR->Fill(dzSV,dr); // leading
	    hLeadingdR_eitherJets->Fill(dr_corr), hLeadingdR_leadingJets->Fill(dr);
	    if(sortedJets.size()>1) hLeadingdR_sleadingJets->Fill(dr_sub);
	    if(fabs(leadJetEta_)<2.0) {
	       hLeadingTrkPVdZ_narrowEta->Fill(dzPV), hLeadingTrkSVdZ_narrowEta->Fill(dzSV); 
	       if(sortedJets.size()>1) hLeadingTrkPVdZ_tightest->Fill(dzPV), hLeadingTrkSVdZ_tightest->Fill(dzSV);
	    }
	 }
	 if(j==1) {
	    hSLeadingTrkPVdZ->Fill(dzPV), hSLeadingTrkSVdZ->Fill(dzSV), hSLeadingPVdZdR->Fill(dzPV,dr), hSLeadingSVdZdR->Fill(dzSV,dr);
	    hSLeadingdR_eitherJets->Fill(dr_corr), hSLeadingdR_leadingJets->Fill(dr);
	    if(sortedJets.size()>1) hSLeadingdR_sleadingJets->Fill(dr_sub);
	    hLeadingAndSLeadingTrkPVdZ->Fill(dzPV_leading,dzPV), hLeadingAndSLeadingTrkSVdZ->Fill(dzSV_leading,dzSV);
	 }
	 if(j==2) {
	    hSSLeadingTrkPVdZ->Fill(dzPV), hSSLeadingTrkSVdZ->Fill(dzSV), hSSLeadingPVdZdR->Fill(dzPV,dr), hSSLeadingSVdZdR->Fill(dzSV,dr);
	    hSSLeadingdR_eitherJets->Fill(dr_corr), hSSLeadingdR_leadingJets->Fill(dr);
	    if(sortedJets.size()>1) hSSLeadingdR_sleadingJets->Fill(dr_sub);
	    hLeadingAndSSLeadingTrkPVdZ->Fill(dzPV_leading,dzPV), hLeadingAndSSLeadingTrkSVdZ->Fill(dzSV_leading,dzSV);
	    hSLeadingAndSSLeadingTrkPVdZ->Fill(dzPV_sleading,dzPV), hSLeadingAndSSLeadingTrkSVdZ->Fill(dzSV_sleading,dzSV);
	 }
	 
	 hTrkPtFromAV->Fill(sortedTrks[j]->pt());

	 if(fabs(dzPV)<=dzcut_) { // PV originated tracks
	    hTrkPtFromPV->Fill(sortedTrks[j]->pt());
	    if(j==0) {
	       hLeadingJetEtaFromPV->Fill(leadJetEta_), hLeadingJetEtFromPV->Fill(leadJetEt_), hLeadingdRFromPV->Fill(dr); 
	    }
	    if(j==1){
	       hSLeadingJetEtaFromPV->Fill(leadJetEta_), hSLeadingJetEtFromPV->Fill(leadJetEt_), hSLeadingdRFromPV->Fill(dr);
	    }
	    if(j==2){
	       hSSLeadingJetEtaFromPV->Fill(leadJetEta_), hSSLeadingJetEtFromPV->Fill(leadJetEt_), hSSLeadingdRFromPV->Fill(dr);
	    }
	    ptSum_PV = ptSum_PV + sortedTrks[j]->pt();
	 }else{ // any non-PV originatd tracks
	    hTrkPtFromSV->Fill(sortedTrks[j]->pt());
	    if(j==0) {
	       hLeadingJetEtaFromSV->Fill(leadJetEta_), hLeadingJetEtFromSV->Fill(leadJetEt_), hLeadingdRFromSV->Fill(dr);
	    }
            if(j==1) {
	       hSLeadingJetEtaFromSV->Fill(leadJetEta_), hSLeadingJetEtFromSV->Fill(leadJetEt_), hSLeadingdRFromSV->Fill(dr);
	    }
            if(j==2) {
	       hSSLeadingJetEtaFromSV->Fill(leadJetEta_), hSSLeadingJetEtFromSV->Fill(leadJetEt_), hSSLeadingdRFromSV->Fill(dr);
	    }
	    ptSum_SV = ptSum_SV + sortedTrks[j]->pt();
	 }
      }

      hTrkIndex->Fill(nthleading);
      dRofMostCorrTrk->Fill(least_dr_corr);
      hMostCorrTrkPVdZ->Fill(dzPV_mostCorr), hMostCorrTrkSVdZ->Fill(dzSV_mostCorr);
      
      // events with 0, 1 jet only --> leading jet eta>2.4
      // events with more then 1 jets --> both leading and sub-leading jet eta > 2.4 
      if((sortedJets.size()<2 && fabs(leadJetEta_)>= 2.4) || 
	 (sortedJets.size()>1 && fabs(leadJetEta_)>= 2.4 && fabs(sleadJetEta_)>=2.4))
	 dRofMostCorrTrk_nonAcc->Fill(least_dr_corr), hMostCorrTrkPVdZ_nonAcc->Fill(dzPV_mostCorr), hMostCorrTrkSVdZ_nonAcc->Fill(dzSV_mostCorr);
      else
	 dRofMostCorrTrk_Acc->Fill(least_dr_corr), hMostCorrTrkPVdZ_Acc->Fill(dzPV_mostCorr), hMostCorrTrkSVdZ_Acc->Fill(dzSV_mostCorr);

      hJetEtTrkPtSumPV->Fill(ptSum_PV,leadJetEt_);
      hJetEtTrkPtSumSV->Fill(ptSum_SV,leadJetEt_);

   } // end of if(jetTrkVerticesCorr_)
}


// ------------ method called once each job just before starting event loop  ------------
void 
VertexAnalyzer::beginJob()
{

   // pt bins
   const double small = 1e-3;
   double ptb;

   // simple rebinning possible with a rebinning factor n = 2, 3, 4 !
   for(ptb =   0  ; ptb <   1.2-small; ptb +=  0.05) ptBins.push_back(ptb); // 24 bins
   for(ptb =   1.2; ptb <   2.4-small; ptb +=  0.1 ) ptBins.push_back(ptb); // 12 bins
   for(ptb =   2.4; ptb <   7.2-small; ptb +=  0.2 ) ptBins.push_back(ptb); // 24 bins
   for(ptb =   7.2; ptb <  13.2-small; ptb +=  0.5 ) ptBins.push_back(ptb); // 12 bins
   for(ptb =  13.2; ptb <  25.2-small; ptb +=  1.0 ) ptBins.push_back(ptb); // 12 bins
   for(ptb =  25.2; ptb <  61.2-small; ptb +=  3.0 ) ptBins.push_back(ptb); // 12 bins
   for(ptb =  61.2; ptb < 121.2-small; ptb +=  5.0 ) ptBins.push_back(ptb); // 12 bins
   for(ptb = 121.2; ptb < 361.2-small; ptb += 10.0 ) ptBins.push_back(ptb); // 24 bins
   ptBins.push_back(361.2);

   // constant pt bins
   static float ptMin;
   static float ptMax;
   static float ptWidth;

   ptMin = 0, ptMax = 1200, ptWidth = 10;

   for(double pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
      cptBins.push_back(pt);


   // jet et bins
   static float jetMin;
   static float jetMax; 
   static float jetWidth;

   jetMin = 0, jetMax = 2000, jetWidth = 20;

   for(double jet = jetMin; jet < jetMax + jetWidth/2; jet += jetWidth)
      jetBins.push_back(jet);


  hVtxSize = f->make<TH1D>("hVtxSize","number of reconstructed vertices",10,-0.5,9.5);
  hVtxTracks = f->make<TH1D>("hVtxTracks","number of tracks fitted to vertex",nTracksBins_,-0.5,(float)nTracksBins_-0.5);
  hVtxZ = f->make<TH1D>("hVtxZ","z position of best reconstructed vertex",80,-20.0,20.0);
  hMultVtxTracks = f->make<TH2D>("hMultVtxTracks","number of tracks fitted to vertex; most populated; 2nd most populated",
				 nTracksBins_,-0.5,(float)nTracksBins_-0.5,nTracksBins_,-0.5,(float)nTracksBins_-0.5);
  hMultVtxZ = f->make<TH2D>("hMultVtxZ","z position of reconstructed vertex; most populated; 2nd most populated",80,-20.0,20.0,80,-20.0,20.0);
  hAllVtxZ = f->make<TH1D>("hAllVtxZ","z position of all reconstructed vertices",80,-20.0,20.0);

  if(jetTrkVerticesCorr_){
     hLeadingTrkPVdZ = f->make<TH1F>("hLeadingTrkPVdZ","dz of leading track with primary vertex",300,-30,30);
     hSLeadingTrkPVdZ = f->make<TH1F>("hSLeadingTrkPVdZ","dz of sub-leading track with primary vertex",300,-30,30);
     hSSLeadingTrkPVdZ = f->make<TH1F>("hSSLeadingTrkPVdZ","dz of sub-sub-leading track with primary vertex",300,-30,30);
     
     hLeadingTrkSVdZ = f->make<TH1F>("hLeadingTrkSVdZ","dz of leading track with 2nd vertex",300,-30,30);
     hSLeadingTrkSVdZ = f->make<TH1F>("hSLeadingTrkSVdZ","dz of sub-leading track with 2nd vertex",300,-30,30);
     hSSLeadingTrkSVdZ = f->make<TH1F>("hSSLeadingTrkSVdZ","dz of sub-sub-leading track with 2nd vertex",300,-30,30);

     hLeadingTrkPVdZ_tightest = f->make<TH1F>("hLeadingTrkPVdZ_tightest","dz of leading track with primary vertex",300,-30,30);
     hLeadingTrkSVdZ_tightest = f->make<TH1F>("hLeadingTrkSVdZ_tightest","dz of leading track with 2nd vertex",300,-30,30);

     hMostCorrTrkPVdZ = f->make<TH1F>("hMostCorrTrkPVdZ","dz of most leading jet-correlated tracks with PV",300,-30,30);
     hMostCorrTrkSVdZ = f->make<TH1F>("hMostCorrTrkSVdZ","dz of most leading jet-correlated tracks with non-PV",300,-30,30);

     hMostCorrTrkPVdZ_Acc = f->make<TH1F>("hMostCorrTrkPVdZ_Acc","dz of most leading jet-correlated tracks with PV when jets inside acc"
                                             ,300,-30,30);
     hMostCorrTrkSVdZ_Acc =f->make<TH1F>("hMostCorrTrkSVdZ_Acc","dz of most leading jet-correlated tracks with non-PV when jets inside acc",
                                            300,-30,30);
     hMostCorrTrkPVdZ_nonAcc = f->make<TH1F>("hMostCorrTrkPVdZ_nonAcc","dz of most leading jet-correlated tracks with PV when jets outside acc"
					     ,300,-30,30);
     hMostCorrTrkSVdZ_nonAcc =f->make<TH1F>("hMostCorrTrkSVdZ_nonAcc","dz of most leading jet-correlated tracks with non-PV when jets outside acc",
					    300,-30,30);
     hTrkIndex = f->make<TH1F>("hTrkIndex","nth leading tracks that is most correlated wit jet",10,-0.5,9.5);
     dRofMostCorrTrk = f->make<TH1F>("dRofMostCorrTrk","dr of most jet-correlated leading track in event",80,0.0,8.0);
     dRofMostCorrTrk_nonAcc = f->make<TH1F>("dRofMostCorrTrk_nonAcc","dr of most jet-correlated leading track in event",90,0.0,8.0);
     dRofMostCorrTrk_Acc = f->make<TH1F>("dRofMostCorrTrk_Acc","dr of most jet-correlated leading track in event",80,0.0,9.0);
     dEtadPhi_leading = f->make<TH2F>("dEtadPhi_leading","deta vs dphi of leading track", 60,-9.0,9.0, 22,0.0,2.*TMath::Pi());

     hLeadingAndSLeadingTrkPVdZ = f->make<TH2F>("hLeadingAndSLeadingTrkPVdZ","dz of leading track vs dz of sub-leading tracks",60,-30,30, 60,-30,30);
     hLeadingAndSSLeadingTrkPVdZ = f->make<TH2F>("hLeadingAndSSLeadingTrkPVdZ","dz of leading track vs dz of ssub-leading tracks",60,-30,30, 60,-30,30);
     hSLeadingAndSSLeadingTrkPVdZ = f->make<TH2F>("hSLeadingAndSSLeadingTrkPVdZ","dz of sub-leading track vs dz of ssub-leading tracks",60,-30,30, 60,-30,30);

     hLeadingAndSLeadingTrkSVdZ = f->make<TH2F>("hLeadingAndSLeadingTrkSVdZ","dz of leading track vs dz of sub-leading tracks",60,-30,30, 60,-30,30);
     hLeadingAndSSLeadingTrkSVdZ = f->make<TH2F>("hLeadingAndSSLeadingTrkSVdZ","dz of leading track vs dz of ssub-leading tracks",60,-30,30, 60,-30,30);
     hSLeadingAndSSLeadingTrkSVdZ = f->make<TH2F>("hSLeadingAndSSLeadingTrkSVdZ","dz of sub-leading track vs dz of ssub-leading tracks",60,-30,30, 60,-30,30);

     hLeadingTrkPVdZ_narrowEta = f->make<TH1F>("hLeadingTrkPVdZ_narrowEta","dz of leading track with primary vertex",300,-30,30);
     hLeadingTrkSVdZ_narrowEta = f->make<TH1F>("hLeadingTrkSVdZ_narrowEta","dz of leading track with 2nd vertex",300,-30,30);

     hLeadingJetEtaFromPV = f->make<TH1F>("hLeadingJetEtaFromPV","leading jet eta whean l-track is from PV",60,-5.0,5.0);
     hLeadingJetEtaFromSV = f->make<TH1F>("hLeadingJetEtaFromSV","leading jet eta whean l-track is from SV",60,-5.0,5.0);
     hSLeadingJetEtaFromPV = f->make<TH1F>("hSLeadingJetEtaFromPV","leading jet eta whean sl-track is from PV",60,-5.0,5.0);
     hSLeadingJetEtaFromSV = f->make<TH1F>("hSLeadingJetEtaFromSV","leading jet eta whean sl-track is from SV",60,-5.0,5.0);
     hSSLeadingJetEtaFromPV = f->make<TH1F>("hSSLeadingJetEtaFromPV","leading jet eta whean ssl-track is from PV",60,-5.0,5.0);
     hSSLeadingJetEtaFromSV = f->make<TH1F>("hSSLeadingJetEtaFromSV","leading jet eta whean ssl-track is from SV",60,-5.0,5.0);

     hLeadingJetEtFromPV = f->make<TH1F>("hLeadingJetEtFromPV","leading jet et whean l-track is from PV", jetBins.size()-1, &jetBins[0]);
     hLeadingJetEtFromSV = f->make<TH1F>("hLeadingJetEtFromSV","leading jet et whean l-track is from SV", jetBins.size()-1, &jetBins[0]);
     hSLeadingJetEtFromPV = f->make<TH1F>("hSLeadingJetEtFromPV","leading jet et whean sl-track is from PV", jetBins.size()-1, &jetBins[0]);
     hSLeadingJetEtFromSV = f->make<TH1F>("hSLeadingJetEtFromSV","leading jet et whean sl-track is from SV", jetBins.size()-1, &jetBins[0]);
     hSSLeadingJetEtFromPV = f->make<TH1F>("hSSLeadingJetEtFromPV","leading jet et whean ssl-track is from PV", jetBins.size()-1, &jetBins[0]);
     hSSLeadingJetEtFromSV = f->make<TH1F>("hSSLeadingJetEtFromSV","leading jet et whean ssl-track is from SV", jetBins.size()-1, &jetBins[0]);
     
     hLeadingdRFromPV = f->make<TH1F>("hLeadingdRFromPV","dr of leading track w.r.t leading jet when from PV", 90,0.0,9.0);
     hLeadingdRFromSV =f->make<TH1F>("hLeadingdRFromSV","dr of leading track w.r.t leading jet when from SV", 90,0.0,9.0);
     hSLeadingdRFromPV =f->make<TH1F>("hSLeadingdRFromPV","dr of s-leading track w.r.t leading jet when from PV", 90,0.0,9.0);
     hSLeadingdRFromSV =f->make<TH1F>("hSLeadingdRFromSV","dr of s-leading track w.r.t leading jet when from SV", 90,0.0,9.0);
     hSSLeadingdRFromPV =f->make<TH1F>("hSSLeadingdRFromPV","dr of ss-leading track w.r.t leading jet when from PV", 90,0.0,9.0);
     hSSLeadingdRFromSV =f->make<TH1F>("hSSLeadingdRFromSV","dr of ss-leading track w.r.t leading jet when from SV", 90,0.0,9.0);
     
     hLeadingdR_eitherJets = f->make<TH1F>("hLeadingdR_eitherJets","dr of leading track w.r.t leading or sub leading jet", 80,0.0,8.0);
     hSLeadingdR_eitherJets = f->make<TH1F>("hSLeadingdR_eitherJets","dr of leading track w.r.t leading or sub leading jet", 80,0.0,8.0);
     hSSLeadingdR_eitherJets = f->make<TH1F>("hSSLeadingdR_eitherJets","dr of leading track w.r.t leading or sub leading jet", 80,0.0,8.0);

     hLeadingdR_leadingJets = f->make<TH1F>("hLeadingdR_leadingJets","dr of leading track w.r.t leading jet", 80,0.0,8.0);
     hSLeadingdR_leadingJets = f->make<TH1F>("hSLeadingdR_leadingJets","dr of leading track w.r.t leading jet", 80,0.0,8.0);
     hSSLeadingdR_leadingJets = f->make<TH1F>("hSSLeadingdR_leadingJets","dr of leading track w.r.t leading jet", 80,0.0,8.0);
     
     hLeadingdR_sleadingJets = f->make<TH1F>("hLeadingdR_sleadingJets","dr of leading track w.r.t sleading jet", 80,0.0,8.0);
     hSLeadingdR_sleadingJets = f->make<TH1F>("hSLeadingdR_sleadingJets","dr of leading track w.r.t sleading jet", 80,0.0,8.0);
     hSSLeadingdR_sleadingJets = f->make<TH1F>("hSSLeadingdR_sleadingJets","dr of leading track w.r.t sleading jet", 80,0.0,8.0);

     hTrkPtFromAV = f->make<TH1F>("hTrkPtFromAV","track p_{T} distribution from all vertices;p_{T} (GeV/c)",ptBins.size()-1, &ptBins[0]);
     hTrkPtFromPV = f->make<TH1F>("hTrkPtFromPV","track p_{T} distribution from PV;p_{T} (GeV/c)",ptBins.size()-1, &ptBins[0]);
     hTrkPtFromSV = f->make<TH1F>("hTrkPtFromSV","track p_{T} distribution from SV;p_{T} (GeV/c)",ptBins.size()-1, &ptBins[0]);

     hVtxZTrkEta = f->make<TH2F>("hVtxZTrkEta","vtx vz vs track eta",80,-20.,20.,24,-3.0,3.0);
     hLeadMonoJetEta = f->make<TH1F>("hLeadMonoJetEta","leading jet eta if monojet", 60,-5.0,5.0);
     hLeadnSLeadJetEta = f->make<TH2F>("hLeadnSLeadJetEta","leading jet eta vs sleading jet eta", 60,-5.0,5.0, 60,-5.0,5.0);

     hLeadingPVdZdR = f->make<TH2F>("hLeadingPVdZdR","dz vs dr of leading track when from PV",60,-30,30, 50,0.0,3.5);
     hSLeadingPVdZdR = f->make<TH2F>("hSLeadingPVdZdR","dz vs dr of s-leading track when from PV",60,-30,30, 50,0.0,3.5);
     hSSLeadingPVdZdR = f->make<TH2F>("hSSLeadingPVdZdR","dz vs dr ss-leading track when from PV",60,-30,30, 50,0.0,3.5);
     hLeadingSVdZdR = f->make<TH2F>("hLeadingSVdZdR","dz vs dr of leading track when from SV",60,-30,30, 50,0.0,3.5);
     hSLeadingSVdZdR = f->make<TH2F>("hSLeadingSVdZdR","dz vs dr of s-leading track when from SV",60,-30,30, 50,0.0,3.5);
     hSSLeadingSVdZdR = f->make<TH2F>("hSSLeadingSVdZdR","dz vs dr ss-leading track when from SV",60,-30,30, 50,0.0,3.5);

     hJetEtTrkPtSumPV = f->make<TH2F>("hJetEtTrkPtSumPV","track sum pt (PV originated) vs leading jet et",cptBins.size()-1, &cptBins[0], jetBins.size()-1, &jetBins[0]);
     hJetEtTrkPtSumSV = f->make<TH2F>("hJetEtTrkPtSumSV","track sum pt (SV originated) vs leading jet et",cptBins.size()-1, &cptBins[0], jetBins.size()-1, &jetBins[0]);
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexAnalyzer);
