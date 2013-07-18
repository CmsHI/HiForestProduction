// -*- C++ -*-
//
// Package:    TrackSpectraAnalyzer
// Class:      TrackSpectraAnalyzer
// 
/**\class TrackSpectraAnalyzer TrackSpectraAnalyzer.cc SpectraAna/TrackSpectraAnalyzer/src/TrackSpectraAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andre Yoon,32 4-A06,+41227676980,
//         Created:  Wed Apr 28 16:18:39 CEST 2010
// $Id: TrackSpectraAnalyzer.cc,v 1.9 2010/05/09 22:01:42 sungho Exp $
//
//


// system include files
#include <memory>
#include <string>
//#include <vector>
//#include <iostream>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "SpectraAna/TrackSpectraAnalyzer/interface/TrackSpectraAnalyzer.h"


TrackSpectraAnalyzer::TrackSpectraAnalyzer(const edm::ParameterSet& iConfig)
   :
   src_(iConfig.getUntrackedParameter<edm::InputTag>("src")),
   vsrc_(iConfig.getUntrackedParameter<edm::InputTag>("vsrc")),
   jsrc_(iConfig.getUntrackedParameter<edm::InputTag>("jsrc")),
   wantNtuple_(iConfig.getUntrackedParameter<bool>("wantNtuple",true)),
   doJet_(iConfig.getUntrackedParameter<bool>("doJet", true)),
   etaCut_(iConfig.getUntrackedParameter<double>("etaCut")),
   nVtxTrkCut_(iConfig.getUntrackedParameter<int>("nVtxTrkCut"))
{
   histograms = new TrackSpectraAnalyzerHist(iConfig);
}


TrackSpectraAnalyzer::~TrackSpectraAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrackSpectraAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   const int nTrigs = 5;

   const string qualityString = "highPurity";

   //-----------------------   (This part will be in an EDFilter later)  
   // get vtx collection 
   Handle<vector<Vertex> > vertices;
   iEvent.getByLabel(vsrc_, vertices);
   if(vertices->size()!=1) return; // one and only one vertex
   int numFake=0, numVtxTrk=0;
   double bestvz=-999.9;
   for(unsigned it=0; it<vertices->size(); ++it) {
      const reco::Vertex & vtx = (*vertices)[it];
      if(vtx.isFake()) numFake++;
      numVtxTrk = vtx.tracksSize();
      bestvz = vtx.z();
   }
   if(numVtxTrk<nVtxTrkCut_) return;
   if(numFake>=1) return;
   //hVtxZ->Fill(bestvz);
   //-----------------------   (This part will be in an EDFilter later)    
   
   // get hlt bit
   Handle<edm::TriggerResults> triggerResults;
   iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"), triggerResults);
   //iEvent.getByLabel("TriggerResults","","HLT",triggerResults);
   const edm::TriggerNames triggerNames = iEvent.triggerNames(*triggerResults); 
   bool accept[nTrigs];
   accept[0]=false; accept[1]=false; accept[2]=false; accept[3]=false; accept[4]=false;
   
   accept[0] = triggerResults->accept(triggerNames.triggerIndex("HLT_MinBiasBSC"));
   accept[1] = triggerResults->accept(triggerNames.triggerIndex("HLT_L1Jet6U"));
   accept[2] = triggerResults->accept(triggerNames.triggerIndex("HLT_Jet15U"));
   accept[3] = triggerResults->accept(triggerNames.triggerIndex("HLT_Jet30U"));
   accept[4] = triggerResults->accept(triggerNames.triggerIndex("HLT_Jet50U"));

   for(int i=0; i<nTrigs; i++) { 
      //cout<<triggerResults->accept(triggerNames.triggerIndex("HLT_MinBiasBSC"))<<endl;
      //accept[i] = triggerResults->accept(triggerNames.triggerIndex(hltNames[i])); 
      //if(accept[i]) hHLTPaths->Fill(hltNames[i],1); 
    } 

   
   
   //----- loop over pat jets and store in a vector -----
   Handle<std::vector<pat::Jet> > pjets;
   iEvent.getByLabel("selectedPatJets", pjets);

   vector<const pat::Jet *> sortedpJets;

   if(doJet_){
      for(unsigned it=0; it<pjets->size(); ++it){
	 const pat::Jet* pjts = &((*pjets)[it]);
	 sortedpJets.push_back( & *pjts);
	 sortByEtRef (&sortedpJets);
      }
      
      for(unsigned it=0; it<sortedpJets.size(); ++it){
	 //nt_jet->Fill(sortedpJets[it]->et(),sortedpJets[it]->eta(),
	 //accept[0],accept[1],accept[2],accept[3],accept[4]); 
	 if(accept[0]) histograms->fillJet2DHist(sortedpJets[it]->eta(),sortedpJets[it]->et(),0);
	 if(accept[1]) histograms->fillJet2DHist(sortedpJets[it]->eta(),sortedpJets[it]->et(),1);
	 if(accept[2]) histograms->fillJet2DHist(sortedpJets[it]->eta(),sortedpJets[it]->et(),2);
	 if(accept[3]) histograms->fillJet2DHist(sortedpJets[it]->eta(),sortedpJets[it]->et(),3);
	 if(accept[4]) histograms->fillJet2DHist(sortedpJets[it]->eta(),sortedpJets[it]->et(),4);
	 break;             
      }                     
   }

   // get track collection 
   Handle<vector<Track> > tracks;
   iEvent.getByLabel(src_, tracks);

   for(unsigned it=0; it<tracks->size(); ++it){
      const reco::Track & trk = (*tracks)[it];

      if(!trk.quality(reco::TrackBase::qualityByName(qualityString))) continue;

      //if(wantNtuple_)nt_dndptdeta->Fill(trk.pt(),trk.eta());
      //else hTrkPtEta->Fill(trk.eta(),trk.pt());
      histograms->fillTrack2DHist(trk.eta(),trk.pt(),0);

      // (leading jet)-track                                  
      // even if there's no jet track info saved (needed for MB) 

      if(doJet_){
         float jet_et = 0, jet_eta = 0;
         unsigned index = 0; 
         if(sortedpJets.size()==0) jet_et = 0,jet_eta = 0; 
	 else jet_et = sortedpJets[index]->et(), jet_eta = sortedpJets[index]->eta(); 
         //if(wantNtuple_) {
	 //if(fabs(jet_eta)<jetaMax_) nt_jettrack->Fill(trk.pt(),trk.eta(),jet_et,
	 //accept[0],accept[1],accept[2],accept[3],accept[4]); 
	 //}else{
	 if(fabs(jet_eta)<etaCut_){
	    if(accept[0]) histograms->fillTrackJet3DHist(trk.eta(),trk.pt(),jet_et,0);
	    if(accept[1]) histograms->fillTrackJet3DHist(trk.eta(),trk.pt(),jet_et,1);
	    if(accept[2]) histograms->fillTrackJet3DHist(trk.eta(),trk.pt(),jet_et,2);
	    if(accept[3]) histograms->fillTrackJet3DHist(trk.eta(),trk.pt(),jet_et,3);
	    if(accept[4]) histograms->fillTrackJet3DHist(trk.eta(),trk.pt(),jet_et,4);
	 }
	 //}
      }                                                       
   }
   
   if(isGEN_){
      Handle<GenParticleCollection> genParticles;
      iEvent.getByLabel("genParticles", genParticles);
      const GenParticleCollection *genCollect = genParticles.product();

      for(unsigned i=0; i<genCollect->size();i++){
	 const GenParticle & gen = (*genCollect)[i];
	 if(gen.status() != 1) continue;
	 if(gen.collisionId() != 0) continue;
	 if(gen.charge() == 0) continue;
	 if(fabs(gen.eta())>etaCut_) continue;
	 //nt_gen_dndptdeta->Fill(gen.pt(),gen.eta());
	 histograms->fillTrack2DHist(gen.eta(),gen.pt(),1);
      }
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackSpectraAnalyzer::beginJob()
{
   histograms->prepareHist();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackSpectraAnalyzer::endJob() {
   histograms->writeHist();
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackSpectraAnalyzer);
