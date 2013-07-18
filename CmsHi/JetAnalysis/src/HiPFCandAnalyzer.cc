// -*- C++ -*-
//
// Package:    HiPFCandAnalyzer
// Class:      HiPFCandAnalyzer
// 
/**\class HiPFCandAnalyzer HiPFCandAnalyzer.cc ana/HiPFCandAnalyzer/src/HiPFCandAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Matt, Nguyen 
//         Created:   Oct  10 2010
// 
//
//


// system include files
#include <memory>

// stl
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ana
#include "CmsHi/JetAnalysis/interface/HiPFCandAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"
#include "TStopwatch.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// constructors and destructor
//
HiPFCandAnalyzer::HiPFCandAnalyzer(const edm::ParameterSet& iConfig)
{
  // Event source
  // Event Info
  pfCandidateLabel_ = iConfig.getParameter<edm::InputTag>("pfCandidateLabel");
  genLabel_ = iConfig.getParameter<edm::InputTag>("genLabel");
  jetLabel_ = iConfig.getParameter<edm::InputTag>("jetLabel");

  pfPtMin_ = iConfig.getParameter<double>("pfPtMin");
  genPtMin_ = iConfig.getParameter<double>("genPtMin");
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");

  // debug
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity", 0);

  doJets_ = iConfig.getUntrackedParameter<bool>("doJets",0);
  doMC_ = iConfig.getUntrackedParameter<bool>("doMC",0);
  skipCharged_ = iConfig.getUntrackedParameter<bool>("skipCharged",0);
}


HiPFCandAnalyzer::~HiPFCandAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
  void
HiPFCandAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  pfEvt_.Clear();

  // Fill PF info

  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel(pfCandidateLabel_,pfCandidates);  
  const reco::PFCandidateCollection *pfCandidateColl = &(*pfCandidates);


  for(unsigned icand=0;icand<pfCandidateColl->size(); icand++) {
      const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);      

      double pt =  pfCandidate.pt();
      if(pt<pfPtMin_) continue;

      int id = pfCandidate.particleId();
      if(skipCharged_ && (abs(id) == 1 || abs(id) == 3)) continue;

      pfEvt_.pfId_[pfEvt_.nPFpart_] = id;      
      pfEvt_.pfPt_[pfEvt_.nPFpart_] = pt;      
      pfEvt_.pfEta_[pfEvt_.nPFpart_] = pfCandidate.eta();      
      pfEvt_.pfPhi_[pfEvt_.nPFpart_] = pfCandidate.phi();      
      pfEvt_.nPFpart_++;
      
  }
	

  // Fill GEN info
  if(doMC_){
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genLabel_,genParticles);     
  const reco::GenParticleCollection* genColl= &(*genParticles);

  for(unsigned igen=0;igen<genColl->size(); igen++) {    

    const reco::GenParticle gen = genColl->at(igen);    
    double eta = gen.eta();      
    double pt = gen.pt();      

    if(gen.status()==1 && fabs(eta)<3.0 && pt> genPtMin_){      
      pfEvt_.genPDGId_[pfEvt_.nGENpart_] = gen.pdgId();      
      pfEvt_.genPt_[pfEvt_.nGENpart_] = pt;      
      pfEvt_.genEta_[pfEvt_.nGENpart_] = eta;      
      pfEvt_.genPhi_[pfEvt_.nGENpart_] = gen.phi();      
      pfEvt_.nGENpart_++;
    }
  }
  }

    // Fill Jet info
  if(doJets_){
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel(jetLabel_,jets);  
  const pat::JetCollection *jetColl = &(*jets);


  for(unsigned ijet=0;ijet<jetColl->size(); ijet++) {
      const pat::Jet jet = jetColl->at(ijet);
      
      double pt =  jet.pt();
      
      if(pt>jetPtMin_){
	pfEvt_.jetPt_[pfEvt_.njets_] = pt;      
	pfEvt_.jetEta_[pfEvt_.njets_] = jet.eta();      
	pfEvt_.jetPhi_[pfEvt_.njets_] = jet.phi();      
	pfEvt_.njets_++;
      }
  }
  }	

  // All done
  pfTree_->Fill();
}

/*
void HiPFCandAnalyzer::FillEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, TreePFCandEventData & tr)
{
  // General Info
  tr.run_	  = iEvent.id().run();
  tr.evt_	  = iEvent.id().event();
  tr.lumi_	  = iEvent.luminosityBlock();

  if(!genOnly_&&sampleType_<10){
    // HI Event info
    edm::Handle<reco::Centrality> cent;
    iEvent.getByLabel(edm::InputTag("hiCentrality"),cent);
    Double_t hf	  = cent->EtHFhitSum();
    // Get Centrality bin
    cbins_ = getCentralityBinsFromDB(iSetup);
    tr.cent_ = cbins_->getBin(hf)*(100./cbins_->getNbins());
  }

  if (isMC_&&sampleType_<10) {
    edm::Handle<edm::GenHIEvent> mchievt;
    iEvent.getByLabel(edm::InputTag("heavyIon"),mchievt);
    tr.b_	  = mchievt->b();
    tr.npart_	  = mchievt->Npart();
    tr.ncoll_	  = mchievt->Ncoll();
  }
}
*/
void HiPFCandAnalyzer::beginJob()
{

  // -- trees --                                                                                                                                                                                                                        
    pfTree_ = fs->make<TTree>("pfTree","dijet tree");
    pfEvt_.SetTree(pfTree_);
    pfEvt_.doMC = doMC_;
    pfEvt_.doJets = doJets_;

    pfEvt_.SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------                                                                                                                                              
void
HiPFCandAnalyzer::endJob() {
  // ===== Done =====                                                                                                                                                                                                                   
  /*  if (verbosity_>=1) {
    cout << endl << "================ Ana Process Summaries =============" << endl;
    cout << " AnaJet: " << jetsrc_ << endl;
    if (refJetType_>=0) cout << " RefJet: " << refjetsrc_ << endl;
    cout << " AnaTrk: " << trksrc_ << endl;
    cout << "# HI Events : "<< numHiEvtSel_<< endl;
    cout << "# Base Events: "<< numEvtSel_ << endl;
    cout << "# Jet Events: "<< numJetEvtSel_<< endl;
  }
  */
}

// constructors
TreePFCandEventData::TreePFCandEventData(){

}


// set branches
void TreePFCandEventData::SetBranches()
{
  // --event level--

  // -- particle info --
  tree_->Branch("nPFpart",&(this->nPFpart_),"nPFpart/I");
  tree_->Branch("pfId",this->pfId_,"pfId[nPFpart]/I");
  tree_->Branch("pfPt",this->pfPt_,"pfPt[nPFpart]/F");
  tree_->Branch("pfEta",this->pfEta_,"pfEta[nPFpart]/F");
  tree_->Branch("pfPhi",this->pfPhi_,"pfPhi[nPFpart]/F");

  // -- jet info --
  if(doJets){
  tree_->Branch("njets",&(this->njets_),"njets/I");
  tree_->Branch("jetPt",this->jetPt_,"jetPt[njets]/F");
  tree_->Branch("jetEta",this->jetEta_,"jetEta[njets]/F");
  tree_->Branch("jetPhi",this->jetPhi_,"jetPhi[njets]/F");
  }

  // -- gen info --
  if(doMC){
  tree_->Branch("nGENpart",&(this->nGENpart_),"nGENpart/I");
  tree_->Branch("genPDGId",this->genPDGId_,"genPDGId[nGENpart]/I");
  tree_->Branch("genPt",this->genPt_,"genPt[nGENpart]/F");
  tree_->Branch("genEta",this->genEta_,"genEta[nGENpart]/F");
  tree_->Branch("genPhi",this->genPhi_,"genPhi[nGENpart]/F");
  }

}
void TreePFCandEventData::Clear()
{
  // event

  nPFpart_      = 0;
  njets_        = 0;
  nGENpart_     = 0;
}

DEFINE_FWK_MODULE(HiPFCandAnalyzer);
