/*
  Based on edwenger/TrkEffAnalyzer
  Modified by Matt Nguyen and Andre Yoon, November 2010

*/

#include "edwenger/HiTrkEffAnalyzer/interface/HiPFCandidateTrackAnalyzer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"


#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "edwenger/HiVertexAnalyzer/interface/HiVertexComparator.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <Math/DistFunc.h>
#include "TMath.h"



#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"


using namespace std;
using namespace edm;
using namespace reco;


HiPFCandidateTrackAnalyzer::HiPFCandidateTrackAnalyzer(const edm::ParameterSet& iConfig) :
   centrality_(0)
{

  inputTagPFCandidates_ = iConfig.getParameter<InputTag>("PFCandidates");
  inputTagVertices_ = iConfig.getParameter<InputTag>("Vertices");
  inputTagSimTracks_ = iConfig.getParameter<InputTag>("SimTracks");
  inputTagTracks_ = iConfig.getParameter<InputTag>("Tracks");

  verbose_ = iConfig.getUntrackedParameter<bool>("verbose",false);
  printBlocks_ = iConfig.getUntrackedParameter<bool>("printBlocks",false);
  thePtMin_ = iConfig.getUntrackedParameter<double>("ptMin",3.0);
  hasSimInfo_ = iConfig.getUntrackedParameter<bool>("hasSimInfo");
  funcCaloComp_ = iConfig.getParameter<std::string>("funcCaloComp");

  applyTrkQCs_ = iConfig.getUntrackedParameter<bool>("applyTrkQCs");
  minHits_ = iConfig.getUntrackedParameter<double>("minHits");
  maxPtErr_ = iConfig.getUntrackedParameter<double>("maxPtErr");
  maxD0_ = iConfig.getUntrackedParameter<double>("maxD0");
  maxDZ_ = iConfig.getUntrackedParameter<double>("maxDZ");
  maxD0Norm_ = iConfig.getUntrackedParameter<double>("maxD0Norm");
  maxDZNorm_ = iConfig.getUntrackedParameter<double>("maxDZNorm");
  pixelSeedOnly_ = iConfig.getUntrackedParameter<bool>("pixelSeedOnly");
  
  prodNtuple_ = iConfig.getUntrackedParameter<bool>("prodNtuple");
  isData_ = iConfig.getUntrackedParameter<bool>("isData");

  neededCentBins_ = iConfig.getUntrackedParameter<std::vector<int> >("neededCentBins");

  useQaulityStr_ = iConfig.getUntrackedParameter<bool>("useQaulityStr");
  qualityString_ = iConfig.getUntrackedParameter<std::string>("qualityString");

  fCaloComp = new TF1("fCaloComp",funcCaloComp_.c_str(),0,200); // a function that defines track-calo (in)compatible region

  LogDebug("HiPFCandidateTrackAnalyzer")
     <<" input collection : "<<inputTagPFCandidates_<<"\n"
     <<" funtional form : "<<funcCaloComp_<<endl;

}



HiPFCandidateTrackAnalyzer::~HiPFCandidateTrackAnalyzer() { }



void 
HiPFCandidateTrackAnalyzer::beginRun(const edm::Run& run, 
			      const edm::EventSetup & es) {}



// ------------ main analyzer run for each event ---------------------------------------- 

void 
HiPFCandidateTrackAnalyzer::analyze(const Event& iEvent, 
			     const EventSetup& iSetup) {
  
  LogDebug("HiPFCandidateTrackAnalyzer")<<"START event: "<<iEvent.id().event()
			 <<" in run "<<iEvent.id().run()<<endl;
  
  // get PF candidates
  Handle<PFCandidateCollection> pfCandidates;
  bool isPFThere = iEvent.getByLabel(inputTagPFCandidates_, pfCandidates);

  if (!isPFThere) return;  // if no PFCand in an event, skip it 

  // get centrality information
  centrality_ = new CentralityProvider(iSetup);
  centrality_->newEvent(iEvent,iSetup);
  int cbin = centrality_->getBin();

  // get vertices
  edm::Handle<reco::VertexCollection> vtxsH;
  iEvent.getByLabel(inputTagVertices_,vtxsH);
  reco::VertexCollection vtxs = *vtxsH;
  
  // do reco-to-sim association
  edm::Handle<TrackingParticleCollection>  TPCollectionHfake;
  edm::Handle<edm::View<reco::Track> >  trackCollection;
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  const TrackAssociatorByHits *theAssociatorByHits;
  reco::RecoToSimCollection recSimColl;

  // get track collection
  Handle<edm::View<reco::Track> > trackCollectionH;
  iEvent.getByLabel(inputTagTracks_,trackCollectionH);
  const edm::View<reco::Track>  tC = *(trackCollectionH.product());
  //size_t ntrk = tC.size();

  if(hasSimInfo_) {
    iEvent.getByLabel(inputTagSimTracks_,TPCollectionHfake);
    iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
    theAssociatorByHits = (const TrackAssociatorByHits*) theAssociator.product();  
    recSimColl= theAssociatorByHits->associateRecoToSim(trackCollectionH,TPCollectionHfake,&iEvent); // to find fake
  }
  
  // analyzer PF candidates
  if(verbose_) cout<<" # of pf cands: "<<pfCandidates->size()<<endl;

  for( unsigned i=0; i<pfCandidates->size(); i++ ) {
    
    max_nhits=0.0;
    max_relpterr=0.0;
    max_algo=0.0;
    max_nd0=0.0;
    max_ndz=0.0;
    max_fake=0.0;

    const reco::PFCandidate& cand = (*pfCandidates)[i];
    
    cand_type = cand.particleId();
    cand_pt = cand.pt();
    cand_eta = cand.eta();

    //--------- EDIT: Apr 21, 2010 (EAW)
    if(cand_pt < thePtMin_) continue;
    
    // Matt:  Look at only charged hardons for the moment, could be extended to study leptons
    if(!(cand_type == PFCandidate::h)) continue;

    
    // Andre: Look at only pf candidate that is associated with a track from input track collection 
    reco::TrackRef trackRef = cand.trackRef();

    int index = 0;
    bool isAssociatedTrk = false;
    
    for(edm::View<reco::Track>::size_type i=0; i<tC.size(); ++i) {
       if(i==trackRef.key()) {
	  isAssociatedTrk = true;
	  index = i;
	  break;
       }
    }

    if(!isAssociatedTrk) continue;
    
    // get track reference
    edm::RefToBase<reco::Track> track(trackCollectionH, index);
    reco::Track* tr=const_cast<reco::Track*>(track.get());

    if(useQaulityStr_ && !tr->quality(reco::TrackBase::qualityByName(qualityString_))) continue;

    //std::cout<<"associated index = "<<index<<" with pt = "<<trackRef->pt()<<endl;
    
    //----- track quality selections, MATT:  Already using high Purity hiSelectedTracks by default anyway

    double nhits = 0, relpterr = 0, algo = 0, d0 = 0, dz = 0, d0err = 0, dzerr = 0;

    nhits = trackRef->numberOfValidHits();
    relpterr = trackRef->ptError()/trackRef->pt();
    algo = trackRef->algo();
    d0 = trackRef->dxy(vtxs[0].position());
    dz = trackRef->dz(vtxs[0].position());
    d0err = sqrt(trackRef->d0Error()*trackRef->d0Error()+vtxs[0].xError()*vtxs[0].yError());
    dzerr = sqrt(trackRef->dzError()*trackRef->dzError()+vtxs[0].zError()*vtxs[0].zError());

    if(applyTrkQCs_){
       if(nhits<minHits_) continue;
       if(relpterr > maxPtErr_) continue;
       if(algo > 7 && pixelSeedOnly_) continue;
       if(fabs(d0) > maxD0_) continue;
       if(fabs(dz) > maxDZ_) continue;
       if(fabs(d0/d0err) > maxD0Norm_) continue;
       if(fabs(dz/dzerr) > maxDZNorm_) continue;
    }    


    bool fake=false;

    if(hasSimInfo_ && recSimColl.find(track) != recSimColl.end()){
       fake = false;
    } else {
       fake = true;
    }

    //-----
    if(verbose_) if(fake==true) std::cout<<" fake! "<<std::endl;

    max_trk=0.0, max_ecal=0.0, max_hcal=0.0, max_calo=0.0;
    sum_trk=0.0, sum_ecal=0.0, sum_hcal=0.0, sum_calo=0.0;

    
    for(unsigned i=0; i<cand.elementsInBlocks().size(); i++) {

      PFBlockRef blockRef = cand.elementsInBlocks()[i].first;

      
      unsigned indexInBlock = cand.elementsInBlocks()[i].second;
      const edm::OwnVector<  reco::PFBlockElement>&  elements = (*blockRef).elements();

      //This tells you what type of element it is:
      //cout<<" block type"<<elements[indexInBlock].type()<<endl;
      
      switch (elements[indexInBlock].type()) {
	
      case PFBlockElement::ECAL: {
	reco::PFClusterRef clusterRef = elements[indexInBlock].clusterRef();
	double eet = clusterRef->energy()/cosh(clusterRef->eta());
	if(verbose_)cout<<" ecal energy "<<clusterRef->energy()<<endl;
	sum_ecal+=eet;
	if(eet>max_ecal) max_ecal=eet;
	break;
      }
	
      case PFBlockElement::HCAL: {
	reco::PFClusterRef clusterRef = elements[indexInBlock].clusterRef();
	double eet = clusterRef->energy()/cosh(clusterRef->eta());
	sum_hcal+=eet;
	if(eet>max_hcal) max_hcal=eet;
	if(verbose_)cout<<" hcal energy "<<clusterRef->energy()<<endl;
	break; 
      }       
      case PFBlockElement::TRACK: {
	//This is just the reference to the track itself, since tracks can never be linked to other tracks
	break; 
      }       
      default:
	break;
      }

    } // end of elementsInBlocks()
	
    sum_calo = sum_ecal + sum_hcal; // add HCAL and ECAL cal sum
    max_calo = max_ecal + max_hcal;

    double trkpt = trackRef->pt();
    if(verbose_) cout << "pt=" << trkpt << endl;
    
    sum_trk+=trkpt;
    if(trkpt>max_trk) {
      max_trk=trkpt;
      max_nhits=nhits;
      max_relpterr=relpterr;
      max_algo=algo;
      max_nd0=d0/d0err;
      max_ndz=dz/dzerr;
      max_fake=fake;
    }
    
    
    
    if( verbose_ ) {
      cout<<cand<<endl;
      // Matt:  I think you need to protect against empty blocks here
      if (printBlocks_) printElementsInBlocks(cand);
    }   
    
    if(verbose_)cout<<" number of elements in blocks "<<cand.elementsInBlocks().size()<<endl;
    
    if(prodNtuple_) nt->Fill(cand_type,cand_pt,max_trk,sum_trk,max_ecal,sum_ecal,max_hcal,sum_hcal,max_nhits,max_relpterr,max_algo,max_nd0,max_ndz,max_fake);

    // 2D hist
    hTrkPtEcalEtSum->Fill(cand_pt,sum_ecal), hTrkPtHcalEtSum->Fill(cand_pt,sum_hcal), hTrkPtCaloEtSum->Fill(cand_pt,sum_calo);
    
    // Centrality binned 2D hist
    for(unsigned i=0;i<neededCentBins_.size();i++){
       if(i==0){
	  if(cbin<=neededCentBins_[i+1]){
	     hTrkPtEcalEtSum_Cent[i]->Fill(cand_pt,sum_ecal);
	     hTrkPtHcalEtSum_Cent[i]->Fill(cand_pt,sum_hcal);
	     hTrkPtCaloEtSum_Cent[i]->Fill(cand_pt,sum_calo);
	  }
       }else{
	  if(cbin>neededCentBins_[i] && cbin<=neededCentBins_[i+1]){
	     hTrkPtEcalEtSum_Cent[i]->Fill(cand_pt,sum_ecal);
             hTrkPtHcalEtSum_Cent[i]->Fill(cand_pt,sum_hcal);
             hTrkPtCaloEtSum_Cent[i]->Fill(cand_pt,sum_calo);
	  }
       }
    }
       
    float compatible_calo = (fCaloComp->Eval(cand_pt)!=fCaloComp->Eval(cand_pt)) ? 0 : fCaloComp->Eval(cand_pt); // protect agains NaN

    if(compatible_calo < sum_calo){
       hHitsEtaAccept->Fill(cand_eta,nhits), hPtErrEtaAccept->Fill(cand_eta,relpterr), hAlgoEtaAccept->Fill(cand_eta,algo), hD0EtaAccept->Fill(cand_eta,d0);
       hDZEtaAccept->Fill(cand_eta,dz), hD0ErrEtaAccept->Fill(cand_eta,d0err), hDZErrEtaAccept->Fill(cand_eta,dzerr);
       hD0PerErrEtaAccept->Fill(cand_eta,fabs(d0/d0err)), hDZPerErrEtaAccept->Fill(cand_eta,fabs(dz/dzerr));
    }else{
       hHitsEtaReject->Fill(cand_eta,nhits), hPtErrEtaReject->Fill(cand_eta,relpterr), hAlgoEtaReject->Fill(cand_eta,algo), hD0EtaReject->Fill(cand_eta,d0);
       hDZEtaReject->Fill(cand_eta,dz), hD0ErrEtaReject->Fill(cand_eta,d0err), hDZErrEtaReject->Fill(cand_eta,dzerr);
       hD0PerErrEtaReject->Fill(cand_eta,fabs(d0/d0err)), hDZPerErrEtaReject->Fill(cand_eta,fabs(dz/dzerr));
    }

    if(!isData_ && hasSimInfo_){  
       if(fake) hTrkPtEcalEtSum_fake->Fill(cand_pt,sum_ecal), hTrkPtHcalEtSum_fake->Fill(cand_pt,sum_hcal), hTrkPtCaloEtSum_fake->Fill(cand_pt,sum_calo);
       else hTrkPtEcalEtSum_real->Fill(cand_pt,sum_ecal), hTrkPtHcalEtSum_real->Fill(cand_pt,sum_hcal), hTrkPtCaloEtSum_real->Fill(cand_pt,sum_calo);
    }

    // 3D hist
    hTrkPtEtaEcalEtSum->Fill(cand_eta,cand_pt,sum_ecal), hTrkPtEtaHcalEtSum->Fill(cand_eta,cand_pt,sum_hcal), hTrkPtEtaCaloEtSum->Fill(cand_eta,cand_pt,sum_calo);

    if(!isData_ && hasSimInfo_){ // fake only
       if(fake) 
	  hTrkPtEtaEcalEtSum_fake->Fill(cand_eta,cand_pt,sum_ecal), hTrkPtEtaHcalEtSum_fake->Fill(cand_eta,cand_pt,sum_hcal), hTrkPtEtaCaloEtSum_fake->Fill(cand_eta,cand_pt,sum_calo);
       else 
	  hTrkPtEtaEcalEtSum_real->Fill(cand_eta,cand_pt,sum_ecal), hTrkPtEtaHcalEtSum_real->Fill(cand_eta,cand_pt,sum_hcal), hTrkPtEtaCaloEtSum_real->Fill(cand_eta,cand_pt,sum_calo);
    }

  } // end of pfCandidates loop
  

  //std::cout<<"Number of tracks = "<<ntrk<<std::endl;
 
  LogDebug("HiPFCandidateTrackAnalyzer")<<"STOP event: "<<iEvent.id().event()
			 <<" in run "<<iEvent.id().run()<<endl;
}


// ------------ define relevant histograms/tuple here ----------------------------------------

void
HiPFCandidateTrackAnalyzer::beginJob() {

   // eta bins
   static float etaMin   = -2.4;
   static float etaMax   =  2.4;
   static float etaWidth =  0.2;

   for(double eta = etaMin; eta < etaMax + etaWidth/2; eta += etaWidth)
      etaBins.push_back(eta);

   // pt bins
   static float ptMin   =  0.0;
   static float ptMax   =  200.0;
   static float ptWidth =  1.0;

   for(double pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
      ptBins.push_back(pt);
   
   // calo et sum bins
   static float cSumEtMin   =  0.0;
   static float cSumEtMax   =  400.0;
   static float cSumEtWidth =  2.0;

   for(double cSumEt = cSumEtMin; cSumEt < cSumEtMax + cSumEtWidth/2; cSumEt += cSumEtWidth)
      cEtSumBins.push_back(cSumEt);

   // nhits bins
   static float nhitMin   = 5;
   static float nhitMax   =  30;
   static float nhitWidth =  1;

   for(double nhit = nhitMin; nhit < nhitMax + nhitWidth/2; nhit += nhitWidth)
      nhitBins.push_back(nhit);


   if(prodNtuple_) nt = fs->make<TNtuple>("nt","PF Testing","type:pt:tkptmax:tkptsum:eetmax:eetsum:hetmax:hetsum:nhits:relpterr:algo:nd0:ndz:fake");

   hTrkPtEcalEtSum = fs->make<TH2F>("hTrkPtEcalEtSum","pT vs ecal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
   hTrkPtHcalEtSum = fs->make<TH2F>("hTrkPtHcalEtSum","pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
   hTrkPtCaloEtSum = fs->make<TH2F>("hTrkPtCaloEtSum","pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
   
   hTrkPtEtaEcalEtSum = fs->make<TH3F>("hTrkPtEtaEcalEtSum","eta vs pt vs ecal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
				       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
   hTrkPtEtaHcalEtSum = fs->make<TH3F>("hTrkPtEtaHcalEtSum","eta vs pt vs hcal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
				       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
   hTrkPtEtaCaloEtSum = fs->make<TH3F>("hTrkPtEtaCaloEtSum","eta vs pt vs hcal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
				       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
   
   hHitsEtaAccept = fs->make<TH2F>("hHitsEtaAccept","Nhits dist. for accepted tracks; #eta;N_{hits}", etaBins.size()-1, &etaBins[0],nhitBins.size()-1, &nhitBins[0]);
   hHitsEtaReject = fs->make<TH2F>("hHitsEtaReject","Nhits dist. for rejected tracks; #eta;N_{hits}", etaBins.size()-1, &etaBins[0],nhitBins.size()-1, &nhitBins[0]);
   
   hPtErrEtaAccept = fs->make<TH2F>("hPtErrEtaAccept","rel. pT error dist. for accepted tracks; #eta;p_{T}^{err}", etaBins.size()-1,&etaBins[0], 50,0,0.25);
   hPtErrEtaReject = fs->make<TH2F>("hPtErrEtaReject","rel. pT error dist. for rejected tracks; #eta;p_{T}^{err}", etaBins.size()-1,&etaBins[0], 50,0,0.25);
   
   hAlgoEtaAccept = fs->make<TH2F>("hAlgoEtaAccept","algo number dist. for accepted tracks; #eta;Algo", etaBins.size()-1,&etaBins[0], 10,0,10);
   hAlgoEtaReject = fs->make<TH2F>("hAlgoEtaReject","algo number dist. for rejected tracks; #eta;Algo", etaBins.size()-1,&etaBins[0], 10,0,10);

   hD0EtaAccept = fs->make<TH2F>("hD0EtaAccept","D0 dist. for accepted tracks;#eta;d_{0}", etaBins.size()-1,&etaBins[0], 100,-2.0,2.0);
   hD0EtaReject = fs->make<TH2F>("hD0EtaReject","D0 dist. for rejected tracks;#eta;d_{0}", etaBins.size()-1,&etaBins[0], 100,-2.0,2.0);

   hDZEtaAccept = fs->make<TH2F>("hDZEtaAccept","DZ dist. for accepted tracks;#eta;d_{z}", etaBins.size()-1,&etaBins[0], 200,-4,4);
   hDZEtaReject = fs->make<TH2F>("hDZEtaReject","DZ dist. for rejected tracks;#eta;d_{z}", etaBins.size()-1,&etaBins[0], 200,-4,4);

   hD0ErrEtaAccept = fs->make<TH2F>("hD0ErrEtaAccept","D0 error dist. for accepted tracks;#eta;d_{0}^{err}", etaBins.size()-1,&etaBins[0], 100,0,0.2);
   hD0ErrEtaReject = fs->make<TH2F>("hD0ErrEtaReject","D0 error dist. for rejected tracks;#eta;d_{0}^{err}", etaBins.size()-1,&etaBins[0], 100,0,0.2);

   hDZErrEtaAccept = fs->make<TH2F>("hDZErrEtaAccept","DZ error dist. for accepted tracks;#eta;d_{0}^{err}", etaBins.size()-1,&etaBins[0], 100,0,0.2);
   hDZErrEtaReject = fs->make<TH2F>("hDZErrEtaReject","DZ error dist. for rejected tracks;#eta;d_{0}^{err}", etaBins.size()-1,&etaBins[0], 100,0,0.2);

   hD0PerErrEtaAccept = fs->make<TH2F>("hD0PerErrEtaAccept","D0/D0Err dist. for accepted tracks;#eta;d_{0}/#sigma_{err}", etaBins.size()-1,&etaBins[0], 50,0,10);
   hD0PerErrEtaReject = fs->make<TH2F>("hD0PerErrEtaReject","D0/D0Err dist. for rejected tracks;#eta;d_{0}/#sigma_{err}", etaBins.size()-1,&etaBins[0], 50,0,10);

   hDZPerErrEtaAccept = fs->make<TH2F>("hDZPerErrEtaAccept","DZ/DZErr dist. for accepted tracks;#eta;d_{0}/#sigma_{err}", etaBins.size()-1,&etaBins[0],50,0,10);
   hDZPerErrEtaReject = fs->make<TH2F>("hDZPerErrEtaReject","DZ/DZErr dist. for rejected tracks;#eta;d_{0}/#sigma_{err}", etaBins.size()-1,&etaBins[0],50,0,10);

   //centrality binned histograms
   for(unsigned i=0;i<neededCentBins_.size()-1;i++){
      hTrkPtEcalEtSum_Cent.push_back(fs->make<TH2F>("","pT vs ecal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], 
						    cEtSumBins.size()-1, &cEtSumBins[0]));
      hTrkPtHcalEtSum_Cent.push_back(fs->make<TH2F>("","pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], 
						    cEtSumBins.size()-1, &cEtSumBins[0]));
      hTrkPtCaloEtSum_Cent.push_back(fs->make<TH2F>("","pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], 
						    cEtSumBins.size()-1, &cEtSumBins[0]));
      if(i==0){
	 hTrkPtEcalEtSum_Cent[i]->SetName(Form("hTrkPtEcalEtSum_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hTrkPtHcalEtSum_Cent[i]->SetName(Form("hTrkPtHcalEtSum_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	 hTrkPtCaloEtSum_Cent[i]->SetName(Form("hTrkPtCaloEtSum_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
      }else{
	 hTrkPtEcalEtSum_Cent[i]->SetName(Form("hTrkPtEcalEtSum_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hTrkPtHcalEtSum_Cent[i]->SetName(Form("hTrkPtHcalEtSum_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	 hTrkPtCaloEtSum_Cent[i]->SetName(Form("hTrkPtCaloEtSum_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
      }
   }


   if(!isData_ && hasSimInfo_) { 
      hTrkPtEcalEtSum_fake = fs->make<TH2F>("hTrkPtEcalEtSum_fake","fake pT vs ecal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtHcalEtSum_fake = fs->make<TH2F>("hTrkPtHcalEtSum_fake","fake pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtCaloEtSum_fake = fs->make<TH2F>("hTrkPtCaloEtSum_fake","fake pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);

      hTrkPtEcalEtSum_real = fs->make<TH2F>("hTrkPtEcalEtSum_real","real pT vs ecal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtHcalEtSum_real = fs->make<TH2F>("hTrkPtHcalEtSum_real","real pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtCaloEtSum_real = fs->make<TH2F>("hTrkPtCaloEtSum_real","fake pT vs hcal et sum; p_{T} (GeV/c);E_{T} (GeV)",ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);

      hTrkPtEtaEcalEtSum_fake = fs->make<TH3F>("hTrkPtEtaEcalEtSum_fake","eta vs pt vs ecal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
					       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtEtaHcalEtSum_fake = fs->make<TH3F>("hTrkPtEtaHcalEtSum_fake","eta vs pt vs hcal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
					       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtEtaCaloEtSum_fake = fs->make<TH3F>("hTrkPtEtaCaloEtSum_fake","eta vs pt vs hcal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
					       etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);

      hTrkPtEtaEcalEtSum_real = fs->make<TH3F>("hTrkPtEtaEcalEtSum_real","eta vs pt vs ecal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
                                               etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtEtaHcalEtSum_real = fs->make<TH3F>("hTrkPtEtaHcalEtSum_real","eta vs pt vs hcal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
                                               etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);
      hTrkPtEtaCaloEtSum_real = fs->make<TH3F>("hTrkPtEtaCaloEtSum_real","eta vs pt vs hcal et sum; #eta;p_{T} (GeV/c);E_{T} (GeV)",
                                               etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], cEtSumBins.size()-1, &cEtSumBins[0]);

   }


}


// ------------ define functions/method/etc below  ---------------------------------------- 

void HiPFCandidateTrackAnalyzer::printElementsInBlocks(const PFCandidate& cand,
						ostream& out) const {
  if(!out) return;

  PFBlockRef firstRef;

  assert(!cand.elementsInBlocks().empty() );
  for(unsigned i=0; i<cand.elementsInBlocks().size(); i++) {
    PFBlockRef blockRef = cand.elementsInBlocks()[i].first;

    if(blockRef.isNull()) {
      cerr<<"ERROR! no block ref!";
      continue;
    }

    if(!i) {
      out<<(*blockRef);
      firstRef = blockRef;
    }
    else if( blockRef!=firstRef) {
      cerr<<"WARNING! This HiPFCandidate is not made from a single block"<<endl;
    }
 
    out<<"\t"<<cand.elementsInBlocks()[i].second<<endl;
  }
}


DEFINE_FWK_MODULE(HiPFCandidateTrackAnalyzer);
