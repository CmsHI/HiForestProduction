/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010

*/

#include "MNguyen/pi0Analyzer/interface/pi0Analyzer.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <Math/DistFunc.h>
#include "TMath.h"
#include "TLorentzVector.h"



#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

//#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"

using namespace std;
using namespace edm;
using namespace reco;


pi0Analyzer::pi0Analyzer(const edm::ParameterSet& iConfig) {
  
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC",false);
  verbose_ = iConfig.getUntrackedParameter<bool>("verbose",false);
  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality",false);
  pfCandidatesTag_ = iConfig.getParameter<InputTag>("pfCandidatesTag");
}



pi0Analyzer::~pi0Analyzer() { }



void 
pi0Analyzer::beginRun(const edm::Run& run, 
			      const edm::EventSetup & es) {}

void 
pi0Analyzer::beginJob() {

  centrality_ = 0;

  string pi0TagTitle = "Pi0 Analysis Tree"; 
  t = fs1->make<TTree>("pi0Tree","Pi0 Ana Tree");


  t->Branch("run",&pi0s_.run,"run/I");
  t->Branch("evt",&pi0s_.evt,"evt/I");
  t->Branch("lumi",&pi0s_.lumi,"lumi/I");
  t->Branch("b",&pi0s_.b,"b/F");
  t->Branch("vx",&pi0s_.vx,"vx/F");
  t->Branch("vy",&pi0s_.vy,"vy/F");
  t->Branch("vz",&pi0s_.vz,"vz/F");
  t->Branch("hf",&pi0s_.hf,"hf/F");

  t->Branch("bin",&pi0s_.bin,"bin/I");

  t->Branch("npi0s",&pi0s_.npi0s,"npi0s/I");

  t->Branch("invMass",pi0s_.invMass,"invMass[npi0s]/F");
  t->Branch("openAng",pi0s_.openAng,"openAng[npi0s]/F");

  t->Branch("pi0_pt",pi0s_.pi0_pt,"pi0_pt[npi0s]/F");
  t->Branch("pi0_px",pi0s_.pi0_px,"pi0_px[npi0s]/F");
  t->Branch("pi0_py",pi0s_.pi0_py,"pi0_py[npi0s]/F");
  t->Branch("pi0_pz",pi0s_.pi0_pz,"pi0_Pz[npi0s]/F");
  t->Branch("pi0_E",pi0s_.pi0_E,"pi0_E[npi0s]/F");
  t->Branch("pi0_eta",pi0s_.pi0_eta,"pi0_eta[npi0s]/F");

  t->Branch("gam1_pt",pi0s_.gam1_pt,"gam1_pt[npi0s]/F");
  t->Branch("gam1_px",pi0s_.gam1_px,"gam1_px[npi0s]/F");
  t->Branch("gam1_py",pi0s_.gam1_py,"gam1_py[npi0s]/F");
  t->Branch("gam1_pz",pi0s_.gam1_pz,"gam1_pz[npi0s]/F");
  t->Branch("gam1_E",pi0s_.gam1_E,"gam1_E[npi0s]/F");
  t->Branch("gam1_eta",pi0s_.gam1_eta,"gam1_eta[npi0s]/F");
  t->Branch("gam1_merged",pi0s_.gam1_merged,"gam1_merged[npi0s]/I");

  t->Branch("gam2_pt",pi0s_.gam2_pt,"gam2_pt[npi0s]/F");
  t->Branch("gam2_px",pi0s_.gam2_px,"gam2_px[npi0s]/F");
  t->Branch("gam2_py",pi0s_.gam2_py,"gam2_py[npi0s]/F");
  t->Branch("gam2_pz",pi0s_.gam2_pz,"gam2_pz[npi0s]/F");
  t->Branch("gam2_E",pi0s_.gam2_E,"gam2_E[npi0s]/F");
  t->Branch("gam2_eta",pi0s_.gam2_eta,"gam2_eta[npi0s]/F");
  t->Branch("gam2_merged",pi0s_.gam2_merged,"gam2_merged[npi0s]/I");


  
  
}


void 
pi0Analyzer::analyze(const Event& iEvent, 
			     const EventSetup& iSetup) {
  
  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();
  
  pi0s_.run = run;
  pi0s_.evt = event;
  pi0s_.lumi = lumi;

  LogDebug("pi0Analyzer")<<"START event: "<<event<<" in run "<<run<<endl;


 int bin = -1;
  double hf = 0.;
  double b = 999.;

  if(useCentrality_){
      if(!centrality_) centrality_ = new CentralityProvider(iSetup);      
      centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
      //double c = centrality_->centralityValue();
      const reco::Centrality *cent = centrality_->raw();
      
      hf = cent->EtHFhitSum();

      bin = centrality_->getBin();
      b = centrality_->bMean();
  }
   



   
   pi0s_.bin = bin;
   pi0s_.hf = hf;
   pi0s_.b = b;


   edm::Handle<vector<reco::Vertex> >vertex;
   iEvent.getByLabel(edm::InputTag("hiSelectedVertex"), vertex);

   if(vertex->size()>0) {
     pi0s_.vx=vertex->begin()->x();
     pi0s_.vy=vertex->begin()->y();
     pi0s_.vz=vertex->begin()->z();
   }
   

   edm::Handle<PFCandidateCollection> pfCandidates;   
   iEvent.getByLabel(pfCandidatesTag_, pfCandidates);


   
   // FILL Pi0 TREE


   
   for( unsigned ic1=0; ic1<pfCandidates->size(); ic1++ ) {

     const reco::PFCandidate& cand1 = (*pfCandidates)[ic1];
     
     if(cand1.particleId() != PFCandidate::gamma) continue;


     float pt1 = cand1.pt();
     
     if(pt1<0.4) continue;

     if(fabs(cand1.eta())>2) continue;

     for( unsigned ic2=0; ic2<pfCandidates->size(); ic2++ ) {
       
       const reco::PFCandidate& cand2 = (*pfCandidates)[ic2];

       if(ic1==ic2) continue;


       if(cand2.particleId() != PFCandidate::gamma) continue;

       float pt2 = cand2.pt();

       if(pt1<pt2) continue;
       
       if(pt2<0.4) continue;

       if(fabs(cand2.eta())>2) continue;

       TLorentzVector gam1 (cand1.px(), cand1.py(), cand1.pz(), cand1.energy());
       TLorentzVector gam2 (cand2.px(), cand2.py(), cand2.pz(), cand2.energy());

       TLorentzVector pi0 = gam1+gam2;
       



       pi0s_.openAng[pi0s_.npi0s] = acos( gam1.Vect().Unit().Dot( gam2.Vect().Unit() ) );
       
       float pi0Mass = pi0.M();
       if(pi0Mass>0.3) continue;

       pi0s_.invMass[pi0s_.npi0s] = pi0Mass;
              
       pi0s_.pi0_pt[pi0s_.npi0s] = pi0.Pt();
       pi0s_.pi0_px[pi0s_.npi0s] = pi0.Px();
       pi0s_.pi0_py[pi0s_.npi0s] = pi0.Py();
       pi0s_.pi0_pz[pi0s_.npi0s] = pi0.Pz();
       pi0s_.pi0_E[pi0s_.npi0s] = pi0.E();
       pi0s_.pi0_eta[pi0s_.npi0s] = pi0.Eta();

       pi0s_.gam1_pt[pi0s_.npi0s] = gam1.Pt();
       pi0s_.gam1_px[pi0s_.npi0s] = gam1.Px();
       pi0s_.gam1_py[pi0s_.npi0s] = gam1.Py();
       pi0s_.gam1_pz[pi0s_.npi0s] = gam1.Pz();
       pi0s_.gam1_E[pi0s_.npi0s] = gam1.E();
       pi0s_.gam1_eta[pi0s_.npi0s] = gam1.Eta();
       pi0s_.gam1_merged[pi0s_.npi0s] = (int) isMerged(cand1);
       //pi0s_.gam1_merged[pi0s_.npi0s] = 1;


       pi0s_.gam2_pt[pi0s_.npi0s] = gam2.Pt();
       pi0s_.gam2_px[pi0s_.npi0s] = gam2.Px();
       pi0s_.gam2_py[pi0s_.npi0s] = gam2.Py();
       pi0s_.gam2_pz[pi0s_.npi0s] = gam2.Pz();
       pi0s_.gam2_E[pi0s_.npi0s] = gam2.E();
       pi0s_.gam2_eta[pi0s_.npi0s] = gam2.Eta();
       pi0s_.gam2_merged[pi0s_.npi0s] = (int) isMerged(cand2);
       //pi0s_.gam2_merged[pi0s_.npi0s] = 1;
       

       pi0s_.npi0s++;
     }     
   }
   t->Fill();
   cout<<" npi0s "<<pi0s_.npi0s<<endl;
   pi0s_.npi0s=0;
   
}

  

bool pi0Analyzer::isMerged( const reco::PFCandidate& pfc ) const {
  
    
  const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();
  
  typedef PFCandidate::ElementsInBlocks::const_iterator IEB; 
  
  int iEcal = 99999999;
  bool foundEcal = false;
  bool merged = 1; // photon coming from extra energy
  for (IEB ieb=theElements.begin(); ieb<theElements.end(); ++ieb) {
    const PFBlock& block = *(ieb->first);
    const PFBlockElement& elem = block.elements()[ieb->second];
    PFBlockElement::Type type = elem.type();
    
    switch( type ) {
    case PFBlockElement::ECAL:
      if( foundEcal ) {
	//  several ECAL clusters involved in the creation of a photon PFCandidate. this is a merged photon. removing them for now.
	//cout<<"photon coming from several ECAL clusters."<<endl;
// 	cout<<pfc<<endl;
// 	cout<<block<<endl;
        continue; // to avoid photon with several Ecal clusters to be created several times
      }
      foundEcal = true;
      iEcal = ieb->second;
      break;
    default:
      continue;
    }

    // look for tracks linked to this ECAL cluster. 

    std::multimap<double, unsigned> trackElems;
    block.associatedElements( iEcal,  block.linkData(),
			      trackElems ,
			      reco::PFBlockElement::TRACK,
			      reco::PFBlock::LINKTEST_ALL );
    
    if( !trackElems.empty() ) {
      //       cout<<"warning : track linked to ECAL cluster giving a photon."<<endl;
      //       cout<<"photon "<<pfc<<endl;
      //       cout<<block<<endl;
      // continue;
      merged=true;
    }
	
    // look for HCAL cluster linked to this ECAL cluster.
    std::multimap<double, unsigned> hcalElems;
    block.associatedElements( iEcal, block.linkData(),
			      hcalElems ,
			      reco::PFBlockElement::HCAL,
			      reco::PFBlock::LINKTEST_ALL );
    if( !hcalElems.empty()) merged=true;
    
    
    // looking for associated PS clusters
    
  }

  if( !foundEcal ) {
    cout<<"photon, but ECAL cluster not found!"<<endl;
    cout<<pfc<<endl;
    
  }
  
  return merged;
  
  
}

				     					
DEFINE_FWK_MODULE(pi0Analyzer);
