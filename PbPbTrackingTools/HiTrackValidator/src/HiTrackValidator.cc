#include "PbPbTrackingTools/HiTrackValidator/interface/HiTrackValidator.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

HiTrackValidator::HiTrackValidator(const edm::ParameterSet& iConfig)
:
   bslabel_(iConfig.getUntrackedParameter<edm::InputTag>("bslabel")),
   vtxlabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxlabel")),
   trklabel_(iConfig.getUntrackedParameter<edm::InputTag>("trklabel")),
   jetlabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetlabel")),
   simtrklabel_(iConfig.getUntrackedParameter<edm::InputTag>("simtrklabel")),
   towerlabel_(iConfig.getUntrackedParameter<edm::InputTag>("towerlabel")),
   associatorMap_(iConfig.getUntrackedParameter<edm::InputTag>("associatorMap")),
   etaMax_(iConfig.getUntrackedParameter<double>("etaMax")),
   jetEtMin_(iConfig.getUntrackedParameter<double>("jetEtMin")),
   hasSimInfo_(iConfig.getUntrackedParameter<bool>("hasSimInfo")),
   selectFake_(iConfig.getUntrackedParameter<bool>("selectFake")),
   hasCaloMat_(iConfig.getUntrackedParameter<bool>("hasCaloMat")),
   towerPtMin_(iConfig.getUntrackedParameter<double>("towerPtMin",5.0)),
   useQaulityStr_(iConfig.getUntrackedParameter<bool>("useQaulityStr")),
   qualityString_(iConfig.getUntrackedParameter<std::string>("qualityString")),
   fiducialCut_(iConfig.getUntrackedParameter<bool>("fiducialCut",false)),
   funcDeltaRTowerMatch_(iConfig.getParameter<std::string>("funcDeltaRTowerMatch")),
   funcCaloComp_(iConfig.getParameter<std::string>("funcCaloComp")),
   neededCentBins_(iConfig.getUntrackedParameter<std::vector<int> >("neededCentBins")),
   centrality_(0)
{
   // pt dependence of delta R matching requirement
   fDeltaRTowerMatch = new TF1("fDeltaRTowerMatch",funcDeltaRTowerMatch_.c_str(),0,200); 
   // pt dependance of calo compatibility, i.e., minimum sum Calo Et vs. track pT
   fCaloComp = new TF1("fCaloComp",funcCaloComp_.c_str(),0,200); // a parameterization of pt dependent cu
}


// ------------ method called to for each event  ------------
void
HiTrackValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace std;
   using namespace reco;

   //------- Jets --------------------------
   double jetet=0.0;

   if(jetEtMin_>0.0){ // run the validator for events with jet > jetEtMin

      edm::Handle<reco::CandidateView> jets;
      iEvent.getByLabel(jetlabel_,jets);

      vector<const reco::Candidate *> sortedJets;

      for(unsigned it=0; it<jets->size(); ++it){
	 const reco::Candidate* jet = &((*jets)[it]);

	 sortedJets.push_back(jet);
	 sortByEtRef (&sortedJets);
      }

      if(sortedJets.size()!=0) jetet = sortedJets[0]->et(); // leading jet
   }

   if(jetet<jetEtMin_) return;  // skip if leading jet et <jetEtMin_


   //------- Beam spot --------------------
   edm::Handle<reco::BeamSpot> beamsp;
   iEvent.getByLabel(bslabel_, beamsp);
   reco::BeamSpot vertexBeamSpot = *beamsp;

   //------- Vertices ---------------------
   edm::Handle<reco::VertexCollection> vtxsH;
   iEvent.getByLabel(vtxlabel_,vtxsH);

   reco::VertexCollection vtxs = *vtxsH;

   bool isPV = false;
   if(vtxs.size()!=0) isPV = true;
   hVtxSize->Fill(vtxs.size());

   
   //------- Centrality information --------
   centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup);
   int cbin = centrality_->getBin();

   //------ Sim tracks -------------------
   edm::Handle<TrackingParticleCollection>  TPCollectionHfake;
   edm::Handle<edm::View<reco::Track> >  trackCollection;
   edm::ESHandle<TrackAssociatorBase> theAssociator;
   const TrackAssociatorByHits *theAssociatorByHits;
   reco::RecoToSimCollection recSimColl;

   iEvent.getByLabel(trklabel_,trackCollection); // track collection ref base 

   if(hasSimInfo_) {
      iEvent.getByLabel(simtrklabel_,TPCollectionHfake);
      iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
      theAssociatorByHits = (const TrackAssociatorByHits*) theAssociator.product();  
      recSimColl= theAssociatorByHits->associateRecoToSim(trackCollection,TPCollectionHfake,&iEvent); // to find fake
   }

   //------- Calo tower ---------------------
   edm::Handle<CaloTowerCollection> towers;
   bool isTowerThere = iEvent.getByLabel(towerlabel_, towers);

   //------- Get tracker geometry -------------
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
   theTracker = tracker.product();
   
   //------- Tracks -----------------------
   edm::Handle<std::vector<reco::Track> > tracks;
   iEvent.getByLabel(trklabel_, tracks);

   //const std::string qualityString = "highPurity";
   

   for(unsigned it=0; it<tracks->size(); ++it){
      

      const reco::Track & trk = (*tracks)[it];

      if(fiducialCut_ && hitDeadPXF(trk)) continue; // if track hits the dead region, igonore it;
      if(useQaulityStr_ && !trk.quality(reco::TrackBase::qualityByName(qualityString_))) continue;

      // rec to sim association
      edm::RefToBase<reco::Track> track(trackCollection, it);
      std::vector<std::pair<TrackingParticleRef, double> > tp;
      bool isFake = false;

      if(hasSimInfo_ && recSimColl.find(track) != recSimColl.end()){
	 tp = recSimColl[track];
      } else {
	 isFake = true;
      }

      if(selectFake_ && !isFake) continue; // fill histograms for fake tracks only

      // basic kinematic variable
      double pt = trk.pt(), eta = trk.eta(), phi = trk.phi();

      hEtaPhi->Fill(eta,phi);
      hEta->Fill(eta);
      hPhi->Fill(phi);
      hPt->Fill(pt);
      if(!isFake) hPtReal->Fill(pt);

      // centrality binned
      for(unsigned i=0;i<neededCentBins_.size();i++){
	 if(i==0){
	    if(cbin<=neededCentBins_[i+1]){
	       //hEtaPhi_Cent[i]->Fill(eta,phi);
	       hEta_Cent[i]->Fill(eta);
	       hPhi_Cent[i]->Fill(phi);
	       hPt_Cent[i]->Fill(pt);
	    }
	 }else{
	    if(cbin>neededCentBins_[i] && cbin<=neededCentBins_[i+1]){
	       //hEtaPhi_Cent[i]->Fill(eta,phi);
	       hEta_Cent[i]->Fill(eta);
               hPhi_Cent[i]->Fill(phi);
	       hPt_Cent[i]->Fill(pt);
	    }
	 }
      }

      if(etaMax_<eta) continue; // only for a given eta range

      // track - calo cell matching
      if(hasCaloMat_ && isTowerThere){
	 double matchPt = 0.0;
	 double matchDr;
	 matchByDrAllowReuse(trk,towers,matchDr,matchPt);
	 float matchConeRadius_pt = (fDeltaRTowerMatch->Eval(trk.pt())!=fDeltaRTowerMatch->Eval(trk.pt())) ? 0 : fDeltaRTowerMatch->Eval(trk.pt()); 
	 float caloComp_pt = (fCaloComp->Eval(trk.pt())!=fCaloComp->Eval(trk.pt())) ? 0 : fCaloComp->Eval(trk.pt()); 
	 hdR->Fill(matchDr);
	 hdRdPt->Fill(matchDr,pt);
	 hdRdCalE->Fill(matchDr,matchPt);
	 hdPtdCalE->Fill(pt,matchPt);
	 if(matchPt/pt>caloComp_pt && matchDr<=matchConeRadius_pt){
	    hPtCaloMat->Fill(pt);
	    if(!isFake) hPtRealCaloMat->Fill(pt);
	 }
      }

      // basic hit level quality varialbes
      uint32_t nlayers     = trk.hitPattern().trackerLayersWithMeasurement();
      uint32_t nlayers3D   = trk.hitPattern().pixelLayersWithMeasurement() + 
	 trk.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
      uint32_t nlayersLost = trk.hitPattern().trackerLayersWithoutMeasurement();
      uint32_t nhits = trk.numberOfValidHits();

      // basic track level quality variables
      double relpterr = trk.ptError()/pt, chi2n =  trk.normalizedChi2();

      // compatibility variables
      // if no primary vertex is present, calculate dz(d0) w.r.t the beam spot position 
      double dz = (isPV) ? trk.dz(vtxs[0].position()) : trk.dz(vertexBeamSpot.position());
      double d0 = (isPV) ? trk.dxy(vtxs[0].position()) : trk.dxy(vertexBeamSpot.position());

      double dzE =  trk.dzError();
      double d0E =  trk.d0Error();

      // parameterized d0, dz resolution
      double nomd0E = sqrt(99999.*99999.+(99999./max(pt,1e-9))*(99999./max(pt,1e-9)));
      double nomdzE = nomd0E*(std::cosh(eta));
      
      // again if no PV, error from BS
      double vzE = (isPV) ? vtxs[0].zError() : vertexBeamSpot.z0Error();
      double vtE = (isPV) ? sqrt(vtxs[0].xError()*vtxs[0].yError()) 
	 : sqrt(vertexBeamSpot.x0Error()*vertexBeamSpot.y0Error());
      
      double dzErrPV = sqrt(dzE*dzE+vzE*vzE);
      double d0ErrPV = sqrt(d0E*d0E+vtE*vtE);
      
      double dzOverdZErr = dz/dzE;
      double d0Overd0Err = d0/d0E;

      double dzOverdZErrPV = dz/dzErrPV;
      double d0Overd0ErrPV = d0/d0ErrPV;

      // histogram filling
      hNlayers->Fill(nlayers);
      hNlayers3D->Fill(nlayers3D);
      hNlayersLost->Fill(nlayersLost);
      hNvalidHits->Fill(nhits);
      hChi2n->Fill(chi2n);
      hRelPtErr->Fill(relpterr);
      hdzPV->Fill(dz);
      hd0PV->Fill(d0);
      hdzErr->Fill(dzE);
      hd0Err->Fill(d0E);
      hdzErrPV->Fill(dzErrPV);
      hd0ErrPV->Fill(d0ErrPV);
      hvzErr->Fill(vzE);
      hvtErr->Fill(vtE);
      hdzOverdzErr->Fill(dzOverdZErr); 
      hd0Overd0Err->Fill(d0Overd0Err);
      hdzOverdzErrPV->Fill(dzOverdZErrPV);
      hd0Overd0ErrPV->Fill(d0Overd0ErrPV);

      // fill quality bits
      for ( int i = 0; i<reco::TrackBase::qualitySize; ++i){
	 if(trk.quality(reco::TrackBase::qualityByName(reco::TrackBase::qualityName( reco::TrackBase::TrackQuality(i) ).c_str()))){
	    hQualityName->Fill(reco::TrackBase::qualityName( reco::TrackBase::TrackQuality(i) ).c_str(),1);
	 }
      } 

      // histogram filling 2D
      hNlayersdPt->Fill(pt,nlayers);
      hNlayers3DdPt->Fill(pt,nlayers3D);
      hNlayersLostdPt->Fill(pt,nlayersLost);
      hNvalidHitsdPt->Fill(pt,nhits);
      hChi2ndPt->Fill(pt,chi2n);
      hRelPtErrdPt->Fill(pt,relpterr);
      hdzPVdPt->Fill(pt,dz);
      hd0PVdPt->Fill(pt,d0);
      hdzErrdPt->Fill(pt,dzE);
      hd0ErrdPt->Fill(pt,d0E);
      hdzErrPVdPt->Fill(pt,dzErrPV);
      hd0ErrPVdPt->Fill(pt,d0ErrPV);
      hvzErrdPt->Fill(pt,vzE);
      hvtErrdPt->Fill(pt,vtE);
      hdzOverdzErrdPt->Fill(pt,dzOverdZErr);
      hd0Overd0ErrdPt->Fill(pt,d0Overd0Err);
      hdzOverdzErrPVdPt->Fill(pt,dzOverdZErrPV);
      hd0Overd0ErrPVdPt->Fill(pt,d0Overd0ErrPV);

      // correlation
      //hdzErrd0Err->Fill(dzE,d0E);
      //hdzErrd0ErrPV->Fill(dzErrPV,d0ErrPV);
      //hdzOverdzErrd0Err->Fill(dzOverdZErr,d0Overd0Err);
      //hdzOverdzErrd0ErrPV->Fill(dzOverdZErrPV,d0Overd0ErrPV);


      // centrality binned
      for(unsigned i=0;i<neededCentBins_.size();i++){
         if(i==0){
            if(cbin<=neededCentBins_[i+1]){
	       hNlayers_Cent[i]->Fill(nlayers);
	       hNlayers3D_Cent[i]->Fill(nlayers3D);
	       hNlayersLost_Cent[i]->Fill(nlayersLost);
	       hNvalidHits_Cent[i]->Fill(nhits);
	       hChi2n_Cent[i]->Fill(chi2n);
	       hRelPtErr_Cent[i]->Fill(relpterr);
	       hdzPV_Cent[i]->Fill(dz);
               hd0PV_Cent[i]->Fill(d0);
	       hdzErrPV_Cent[i]->Fill(dzErrPV);
	       hd0ErrPV_Cent[i]->Fill(d0ErrPV);
	       hdzOverdzErrPV_Cent[i]->Fill(dzOverdZErrPV);
	       hd0Overd0ErrPV_Cent[i]->Fill(d0Overd0ErrPV);
	       hNvalidHitsdPt_Cent[i]->Fill(pt,nhits);
	       hChi2ndPt_Cent[i]->Fill(pt,chi2n);
	       hRelPtErrdPt_Cent[i]->Fill(pt,relpterr);
	       hdzPVdPt_Cent[i]->Fill(pt,dz);
               hd0PVdPt_Cent[i]->Fill(pt,d0);
	       hdzErrPVdPt_Cent[i]->Fill(pt,dzErrPV);
	       hd0ErrPVdPt_Cent[i]->Fill(pt,d0ErrPV);
	       hdzOverdzErrPVdPt_Cent[i]->Fill(pt,dzOverdZErrPV);
	       hd0Overd0ErrPVdPt_Cent[i]->Fill(pt,d0Overd0ErrPV);
	       //hdzOverdzErrd0Err_Cent[i]->Fill(dzOverdZErr,d0Overd0Err);
	       //hdzOverdzErrd0ErrPV_Cent[i]->Fill(dzOverdZErrPV,d0Overd0ErrPV);
	    }
         }else{
            if(cbin>neededCentBins_[i] && cbin<=neededCentBins_[i+1]){
	       hNlayers_Cent[i]->Fill(nlayers);
	       hNlayers3D_Cent[i]->Fill(nlayers3D);
	       hNlayersLost_Cent[i]->Fill(nlayersLost);
	       hNvalidHits_Cent[i]->Fill(nhits);
	       hChi2n_Cent[i]->Fill(chi2n);
	       hRelPtErr_Cent[i]->Fill(relpterr);
	       hdzPV_Cent[i]->Fill(dz);
               hd0PV_Cent[i]->Fill(d0);
	       hdzErrPV_Cent[i]->Fill(dzErrPV);
	       hd0ErrPV_Cent[i]->Fill(d0ErrPV);
	       hdzOverdzErrPV_Cent[i]->Fill(dzOverdZErrPV);
	       hd0Overd0ErrPV_Cent[i]->Fill(d0Overd0ErrPV);
	       hNvalidHitsdPt_Cent[i]->Fill(pt,nhits);
	       hChi2ndPt_Cent[i]->Fill(pt,chi2n);
	       hRelPtErrdPt_Cent[i]->Fill(pt,relpterr);
	       hdzPVdPt_Cent[i]->Fill(pt,dz);
               hd0PVdPt_Cent[i]->Fill(pt,d0);
	       hdzErrPVdPt_Cent[i]->Fill(pt,dzErrPV);
	       hd0ErrPVdPt_Cent[i]->Fill(pt,d0ErrPV);
	       hdzOverdzErrPVdPt_Cent[i]->Fill(pt,dzOverdZErrPV);
	       hd0Overd0ErrPVdPt_Cent[i]->Fill(pt,d0Overd0ErrPV);
	       //hdzOverdzErrd0Err_Cent[i]->Fill(dzOverdZErr,d0Overd0Err);
	       //hdzOverdzErrd0ErrPV_Cent[i]->Fill(dzOverdZErrPV,d0Overd0ErrPV);
	    }
         }
      }


   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
HiTrackValidator::beginJob()
{

   // pt bins
   const double small = 1e-3;
   double ptb;

   // simple rebinning possible with a rebinning facto n = 2, 3, 6 !
   for(ptb =   0  ; ptb <   1.2-small; ptb +=  0.05) ptBins.push_back(ptb); // 24 bins
   for(ptb =   1.2; ptb <   2.4-small; ptb +=  0.1 ) ptBins.push_back(ptb); // 12 bins
   for(ptb =   2.4; ptb <   7.2-small; ptb +=  0.4 ) ptBins.push_back(ptb); // 12 bins
   for(ptb =   7.2; ptb <  14.4-small; ptb +=  1.2 ) ptBins.push_back(ptb); // 6 bins
   for(ptb =  14.4; ptb <  28.8-small; ptb +=  2.4 ) ptBins.push_back(ptb); // 6 bins 
   for(ptb =  28.8; ptb <  48.0-small; ptb +=  3.2 ) ptBins.push_back(ptb); // 6 bins
   for(ptb =  48.0; ptb <  86.4-small; ptb +=  6.4 ) ptBins.push_back(ptb); // 6 bins
   for(ptb =  86.4; ptb < 189.6-small; ptb +=  8.6 ) ptBins.push_back(ptb); // 6 bins
   ptBins.push_back(189.6);

   //histograms
   int nmaxhits = 30;

   hVtxSize = f->make<TH1F>("hVtxSize","number of veritces; N_{vertex}",5,-0.5,4.5);
   hNlayers = f->make<TH1F>("hNlayers","number of layers with tracker hits; N_{hits}", nmaxhits,0.,(double)nmaxhits);
   hNlayers3D = f->make<TH1F>("hNlayers3D","number of layers with pixel and tracker mono and streohits; N_{hits}", 
			      nmaxhits,0.,(double)nmaxhits);
   hNlayersLost = f->make<TH1F>("hNlayersLost","number of layers with no tracker hits; N_{hits}", nmaxhits,0.,(double)nmaxhits);
   hNvalidHits = f->make<TH1F>("hNvalidHits","number of valid hits; N_{hits}", nmaxhits,0.,(double)nmaxhits);
   hChi2n = f->make<TH1F>("hChi2n","normalized track chi2;chi^{2}/ndof", 60,0.0,5.0);
   hRelPtErr = f->make<TH1F>("hRelPtErr","relative track p_{T} error; p_{T} err/p_{T}", 100,0.0,0.4);
   hdzPV = f->make<TH1F>("hdzPV","dz ; dz", 100,-1.5,1.5);
   hd0PV = f->make<TH1F>("hd0PV","d0 ; d0", 100,-1.5,1.5);
   hdzErr = f->make<TH1F>("hdzErr","dz error; dz error", 80,0.0,0.8);
   hd0Err = f->make<TH1F>("hd0Err","d0 error; d0 error", 100,0.0,0.2);
   hdzErrPV = f->make<TH1F>("hdzErrPV","dz error with vz error summed ; dz error", 80,0.0,0.8);
   hd0ErrPV = f->make<TH1F>("hd0ErrPV","d0 error with vz error summed ; d0 error", 100,0.0,0.2);
   hvzErr = f->make<TH1F>("hvzErr","vz error; vz error", 100,0.0, 0.08);
   hvtErr = f->make<TH1F>("hvtErr","vxy error; vxy error", 200,0.0, 0.08);
   hdzOverdzErr = f->make<TH1F>("hdzOverdzErr","dz/dzError; dz/dzError", 80,-10.0,10.0);
   hd0Overd0Err = f->make<TH1F>("hd0Overd0Err","d0/d0Error; d0/d0Error", 80,-10.0,10.0);
   hdzOverdzErrPV =  f->make<TH1F>("hdzOverdzErrPV","dz/dzError with PV error; dz/dzError", 80,-10.0,10.0);
   hd0Overd0ErrPV = f->make<TH1F>("hd0Overd0ErrPV","d0/d0Error with PV error; d0/d0Error", 80,-10.0,10.0);

   // as a function of pt
   hNlayersdPt = f->make<TH2F>("hNlayersdPt","number of layers with tracker hits vs p_{T}; p_{T};N_{hits}", 150,0.0,150.0, nmaxhits,0.,(double)nmaxhits);
   hNlayers3DdPt = f->make<TH2F>("hNlayers3DdPt","number of layers with tracker hits vs p_{T}; p_{T};N_{hits}", 150,0.0,150.0, nmaxhits,0.,(double)nmaxhits);
   hNlayersLostdPt = f->make<TH2F>("hNlayersLostdPt","number of layers with no tracker hits vs p_{T}; p_{T}; N_{hits}", 150,0.0,150.0, nmaxhits,0.,(double)nmaxhits);
   hNvalidHitsdPt = f->make<TH2F>("hNvalidHitsdPt","number of valid hits vs p_{T}; p_{T}; N_{hits}", 150,0.0,150.0, nmaxhits,0.,(double)nmaxhits);
   hChi2ndPt = f->make<TH2F>("hChi2ndPt","normalized track chi2 vs p_{T}; p_{T};chi^{2}/ndofnormalized track", 150,0.0,150.0, 60,0.0,5.0);
   hRelPtErrdPt = f->make<TH2F>("hRelPtErrdPt","relative track p_{T} error vs p_{T}; p_{T}; p_{T} err/p_{T}", 150,0.0,150.0, 100,0.0,0.4);
   hdzPVdPt = f->make<TH2F>("hdzPVdPt","dz vs p_{T}; p_{T}; dz", 150,0.0,150.0, 100,-1.5,1.5);
   hd0PVdPt = f->make<TH2F>("hd0PVdPt","d0 vs p_{T}; p_{T}; d0", 150,0.0,150.0, 100,-1.5,1.5);
   hdzErrdPt = f->make<TH2F>("hdzErrdPt","dz error vs p_{T}; p_{T}; dz error", 150,0.0,150.0, 80,0.0,0.8);
   hd0ErrdPt = f->make<TH2F>("hd0ErrdPt","d0 error vs p_{T}; p_{T}; d0 error", 150,0.0,150.0, 100,0.0,0.2);
   hdzErrPVdPt = f->make<TH2F>("hdzErrPVdPt","dz error with vz error summed vs p_{T}; p_{T}; dz error", 150,0.0,150.0, 80,0.0,0.8);
   hd0ErrPVdPt = f->make<TH2F>("hd0ErrPVdPt","d0 error with vz error summed vs p_{T}; p_{T}; d0 error", 150,0.0,150.0, 100,0.0,0.2);
   hvzErrdPt =  f->make<TH2F>("hvzErrdPt","vz error vs p_{T}; p_{T}; vz error", 150,0.0,150.0, 100,0.0,0.08);
   hvtErrdPt = f->make<TH2F>("hvtErrdPt","vxy error vs p_{T}; p_{T}; vxy error", 150,0.0,150.0, 200,0.0,0.08);
   hdzOverdzErrdPt = f->make<TH2F>("hdzOverdzErrdPt","dz/dzError vs p_{T}; p_{T}; dz/dzError", 150,0.0,150.0, 80,-10.0,10.0);
   hd0Overd0ErrdPt = f->make<TH2F>("hd0Overd0ErrdPt","d0/d0Error vs p_{T}; p_{T}; d0/d0Error", 150,0.0,150.0, 80,-10.0,10.0);
   hdzOverdzErrPVdPt = f->make<TH2F>("hdzOverdzErrPVdPt","dz/dzError with PV error vs p_{T}; p_{T}; dz/dzError", 150,0.0,150.0, 80,-10.0,10.0);
   hd0Overd0ErrPVdPt = f->make<TH2F>("hd0Overd0ErrPVdPt","d0/d0Error with PV error vs p_{T}; p_{T}; d0/d0Error", 150,0.0,150.0, 80,-10.0,10.0);
   
   // correlation
   //hdzErrd0Err = f->make<TH2D>("hdzErrd0Err","dz error vs d0 error;dz error;d0 error", 80,0.0,0.8, 80,0.0,0.8);
   //hdzErrd0ErrPV = f->make<TH2D>("hdzErrd0ErrPVt","dz error vs d0 error (vz error summed);dz error;d0 error", 80,0.0,0.8, 80,0.0,0.8);
   //hdzOverdzErrd0Err = f->make<TH2D>("hdzOverdzErrd0Err","dz/dzError vs d0/d0Error", 80,-10.0,10.0, 80,-10.0,10.0);
   //hdzOverdzErrd0ErrPV = f->make<TH2D>("hdzOverdzErrd0ErrPV","dz/dzError with PV error vs d0/d0Error with PV error", 80,-10.0,10.0, 80,-10.0,10.0);

   // kinematic distributions
   hEtaPhi = f->make<TH2F>("hEtaPhi","eta vs phi;#eta;#phi", 40,-2.65,2.65, 80,-1.05*TMath::Pi(),1.05*TMath::Pi());
   hEta = f->make<TH1F>("hEta","eta distribution; #eta", 60,-2.65,2.65);
   hPhi = f->make<TH1F>("hPhi","phi distribution; #phi", 80,-1.05*TMath::Pi(),1.05*TMath::Pi());
   hPt = f->make<TH1F>("hPt","pt distribution; p_{T} [GeV/c]", ptBins.size()-1, &ptBins[0]);

   hQualityName = f->make<TH1F>("hQualityName","Quality Bits",10,0.0,10);

   // calo-matching 
   hdR = f->make<TH1F>("hdR","dr between track and calo; dr", 100,0.0,4.0);
   hPtCaloMat = f->make<TH1F>("hPtCaloMat","pt distribution of calo-matched track; p_{T} [GeV/c]", ptBins.size()-1, &ptBins[0]);
   hPtReal = f->make<TH1F>("hPtReal","pt distribution of real track; p_{T} [GeV/c]", ptBins.size()-1, &ptBins[0]);
   hPtRealCaloMat = f->make<TH1F>("hPtRealCaloMat","pt distribution of calo-matched real track; p_{T} [GeV/c]", ptBins.size()-1, &ptBins[0]);
   hdRdPt = f->make<TH2F>("hdRdPt","dr vs pt;dr;p_{T} [GeV/c]", 100,0.0,4.0, 100,0.0,200);
   hdRdCalE = f->make<TH2F>("hdRdCalE","dr vs calo cell pt;dr;Calo tower p_{T} [GeV/c]", 100,0.0,4.0, 100,0.0,200);
   hdPtdCalE = f->make<TH2F>("hdPtdCalE","pt vs calo cell pt;p_{T} [GeV/c]; Calo tower  p_{T} [GeV/c]", 100,0.0,200, 100,0.0,200);


   // centrality binned histogram 
   for(unsigned i=0;i<neededCentBins_.size()-1;i++){

      hNlayers_Cent.push_back(f->make<TH1F>("","number of layers with tracker hits; N_{hits}", nmaxhits,0.,(double)nmaxhits));
      hNlayers3D_Cent.push_back(f->make<TH1F>("","number of layers with pixel and tracker mono and streohits; N_{hits}",nmaxhits,0.,(double)nmaxhits));
      hNlayersLost_Cent.push_back(f->make<TH1F>("","number of layers with no tracker hits; N_{hits}", nmaxhits,0.,(double)nmaxhits));
      hNvalidHits_Cent.push_back(f->make<TH1F>("","number of valid hits; N_{hits}", nmaxhits,0.,(double)nmaxhits));
      hChi2n_Cent.push_back(f->make<TH1F>("","normalized track chi2;chi^{2}/ndof", 60,0.0,5.0));
      hRelPtErr_Cent.push_back(f->make<TH1F>("","relative track p_{T} error; p_{T} err/p_{T}", 100,0.0,0.4));
      hdzPV_Cent.push_back(f->make<TH1F>("","dz ; dz ", 100,-1.5,1.5));
      hd0PV_Cent.push_back(f->make<TH1F>("","d0 ; d0 ", 100,-1.5,1.5));
      hdzErrPV_Cent.push_back(f->make<TH1F>("","dz error with vz error summed ; dz error", 80,0.0,0.8));
      hd0ErrPV_Cent.push_back(f->make<TH1F>("","d0 error; d0 error", 100,0.0,0.2));
      hdzOverdzErrPV_Cent.push_back(f->make<TH1F>("","dz/dzError with PV error; dz/dzError", 80,-10.0,10.0));
      hd0Overd0ErrPV_Cent.push_back(f->make<TH1F>("","d0/d0Error with PV error; d0/d0Error", 80,-10.0,10.0));

      hNvalidHitsdPt_Cent.push_back(f->make<TH2F>("","number of valid hits vs p_{T}; p_{T}; N_{hits}",150,0.0,150.0, nmaxhits,0.,(double)nmaxhits));
      hChi2ndPt_Cent.push_back(f->make<TH2F>("","normalized track chi2 vs p_{T}; p_{T};chi^{2}/ndofnormalized track", 150,0.0,150.0, 60,0.0,5.0));
      hRelPtErrdPt_Cent.push_back(f->make<TH2F>("","relative track p_{T} error vs p_{T}; p_{T}; p_{T} err/p_{T}", 150,0.0,150.0, 100,0.0,0.4));
      hdzPVdPt_Cent.push_back(f->make<TH2F>("","dz vs p_{T}; p_{T}; dz", 150,0.0,150.0, 100,-1.5,1.5));
      hd0PVdPt_Cent.push_back(f->make<TH2F>("","d0 vs p_{T}; p_{T}; d0r", 150,0.0,150.0, 100,-1.5,1.5));
      hdzErrPVdPt_Cent.push_back(f->make<TH2F>("","dz error with vz error summed vs p_{T}; p_{T}; dz error", 150,0.0,150.0, 80,0.0,0.8));
      hd0ErrPVdPt_Cent.push_back(f->make<TH2F>("","d0 error with vz error summed vs p_{T}; p_{T}; d0 error", 150,0.0,150.0, 100,0.0,0.2));
      hdzOverdzErrPVdPt_Cent.push_back(f->make<TH2F>("","dz/dzError with PV error vs p_{T}; p_{T}; dz/dzError", 150,0.0,150.0, 80,-10.0,10.0));
      hd0Overd0ErrPVdPt_Cent.push_back(f->make<TH2F>("","d0/d0Error with PV error vs p_{T}; p_{T}; d0/d0Error", 150,0.0,150.0, 80,-10.0,10.0));

      //hdzOverdzErrd0Err_Cent.push_back(f->make<TH2D>("","dz/dzError vs d0/d0Error", 80,-10.0,10.0, 80,-10.0,10.0));
      //hdzOverdzErrd0ErrPV_Cent.push_back(f->make<TH2D>("","dz/dzError with PV error vs d0/d0Error with PV error", 80,-10.0,10.0, 80,-10.0,10.0));
      //hEtaPhi_Cent.push_back(f->make<TH2D>("","eta vs phi;#eta;#phi", 20,-2.65,2.65, 40,-1.05*TMath::Pi(),1.05*TMath::Pi()));

      hEta_Cent.push_back(f->make<TH1F>("","eta distribution; #eta", 60,-2.65,2.65));
      hPhi_Cent.push_back(f->make<TH1F>("","phi distribution; #phi", 80,-1.05*TMath::Pi(),1.05*TMath::Pi()));
      hPt_Cent.push_back(f->make<TH1F>("","pt distribution; p_{T} [GeV/c]",100,0.0,20));
      
      if(i==0){
	 hNlayers_Cent[i]->SetName(Form("hNlayers_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hNlayers3D_Cent[i]->SetName(Form("hNlayers3D_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hNlayersLost_Cent[i]->SetName(Form("hNlayersLost_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hNvalidHits_Cent[i]->SetName(Form("hNvalidHits_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hChi2n_Cent[i]->SetName(Form("hChi2n_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hRelPtErr_Cent[i]->SetName(Form("hRelPtErr_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hdzPV_Cent[i]->SetName(Form("hdzPV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hd0PV_Cent[i]->SetName(Form("hd0PV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hdzErrPV_Cent[i]->SetName(Form("hdzErrPV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hd0ErrPV_Cent[i]->SetName(Form("hd0ErrPV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hdzOverdzErrPV_Cent[i]->SetName(Form("hdzOverdzErrPV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hd0Overd0ErrPV_Cent[i]->SetName(Form("hd0Overd0ErrPV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 
	 hNvalidHitsdPt_Cent[i]->SetName(Form("hNvalidHitsdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hChi2ndPt_Cent[i]->SetName(Form("hChi2ndPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hRelPtErrdPt_Cent[i]->SetName(Form("hRelPtErrdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hdzPVdPt_Cent[i]->SetName(Form("hdzPVdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
         hd0PVdPt_Cent[i]->SetName(Form("hd0PVdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hdzErrPVdPt_Cent[i]->SetName(Form("hdzErrPVdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hd0ErrPVdPt_Cent[i]->SetName(Form("hd0ErrPVdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hdzOverdzErrPVdPt_Cent[i]->SetName(Form("hdzOverdzErrPVdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hd0Overd0ErrPVdPt_Cent[i]->SetName(Form("hd0Overd0ErrPVdPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));

	 //hdzOverdzErrd0Err_Cent[i]->SetName(Form("hdzOverdzErrd0Err_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 //hdzOverdzErrd0ErrPV_Cent[i]->SetName(Form("hdzOverdzErrd0ErrPV_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 //hEtaPhi_Cent[i]->SetName(Form("hEtaPhi_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));

	 hEta_Cent[i]->SetName(Form("hEta_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hPhi_Cent[i]->SetName(Form("hPhi_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hPt_Cent[i]->SetName(Form("hPt_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));

      }else{
	 hNlayers_Cent[i]->SetName(Form("hNlayers_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hNlayers3D_Cent[i]->SetName(Form("hNlayers3D_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hNlayersLost_Cent[i]->SetName(Form("hNlayersLost_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hNvalidHits_Cent[i]->SetName(Form("hNvalidHits_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hChi2n_Cent[i]->SetName(Form("hChi2n_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hRelPtErr_Cent[i]->SetName(Form("hRelPtErr_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hdzPV_Cent[i]->SetName(Form("hdzPV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hd0PV_Cent[i]->SetName(Form("hd0PV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hdzErrPV_Cent[i]->SetName(Form("hdzErrPV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hd0ErrPV_Cent[i]->SetName(Form("hd0ErrPV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hdzOverdzErrPV_Cent[i]->SetName(Form("hdzOverdzErrPV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hd0Overd0ErrPV_Cent[i]->SetName(Form("hd0Overd0ErrPV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));

	 hNvalidHitsdPt_Cent[i]->SetName(Form("hNvalidHitsdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hChi2ndPt_Cent[i]->SetName(Form("hChi2ndPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hRelPtErrdPt_Cent[i]->SetName(Form("hRelPtErrdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hdzPVdPt_Cent[i]->SetName(Form("hdzPVdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hd0PVdPt_Cent[i]->SetName(Form("hd0PVdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hdzErrPVdPt_Cent[i]->SetName(Form("hdzErrPVdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hd0ErrPVdPt_Cent[i]->SetName(Form("hd0ErrPVdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hdzOverdzErrPVdPt_Cent[i]->SetName(Form("hdzOverdzErrPVdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hd0Overd0ErrPVdPt_Cent[i]->SetName(Form("hd0Overd0ErrPVdPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));

	 //hdzOverdzErrd0Err_Cent[i]->SetName(Form("hdzOverdzErrd0Err_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 //hdzOverdzErrd0ErrPV_Cent[i]->SetName(Form("hdzOverdzErrd0ErrPV_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 //hEtaPhi_Cent[i]->SetName(Form("hEtaPhi_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));

	 hEta_Cent[i]->SetName(Form("hEta_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hPhi_Cent[i]->SetName(Form("hPhi_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
         hPt_Cent[i]->SetName(Form("hPt_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
      }

   }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiTrackValidator::endJob() {
}


//
void 
HiTrackValidator::matchByDrAllowReuse(const reco::Track & trk, const edm::Handle<CaloTowerCollection> & towers, double & bestdr, double & bestpt)
{
   // loop over towers
   bestdr=1e10;
   bestpt=0.;
   for(unsigned int i = 0; i < towers->size(); ++i){
      const CaloTower & tower= (*towers)[i];
      double pt = tower.pt();
      if (pt<towerPtMin_) continue;
      if (fabs(tower.eta())>etaMax_) continue;
      double dr = reco::deltaR(tower,trk);
      if (dr<bestdr) {
	 bestdr = dr;
	 bestpt = pt;
      }
   }
}

//
bool
HiTrackValidator::hitDeadPXF(const reco::Track& tr){

   //-----------------------------------------------
   // For a given track, check whether this contains 
   // hits on the dead region in the forward pixel 
   //-----------------------------------------------

   
   bool hitDeadRegion = false;

   for(trackingRecHit_iterator recHit = tr.recHitsBegin();recHit!= tr.recHitsEnd(); recHit++){

      if((*recHit)->isValid()){

	 DetId detId = (*recHit)->geographicalId();
	 if(!theTracker->idToDet(detId)) continue;

	 Int_t diskLayerNum=0, bladeLayerNum=0, hcylLayerNum=0;
	  
	 unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());

	 if (subdetId == PixelSubdetector::PixelEndcap){
	        
	    PixelEndcapName pxfname(detId.rawId());
	    diskLayerNum = pxfname.diskName();
	    bladeLayerNum = pxfname.bladeName();
	    hcylLayerNum = pxfname.halfCylinder();
	        
	    // hard-coded now based on /UserCode/Appeltel/PixelFiducialRemover/pixelfiducialremover_cfg.py
	    if((bladeLayerNum==4 || bladeLayerNum==5 || bladeLayerNum==6) &&
	       (diskLayerNum==2) && (hcylLayerNum==4)) hitDeadRegion = true;
	 }
	  
      }// end of isValid
   }
   return hitDeadRegion;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiTrackValidator);
