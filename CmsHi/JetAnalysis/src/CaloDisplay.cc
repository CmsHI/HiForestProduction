// -*- C++ -*-
//
// Package:    CaloDisplay
// Class:      CaloDisplay
// 
/**\class CaloDisplay CaloDisplay.cc yetkin/CaloDisplay/src/CaloDisplay.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Wed Oct  3 08:07:18 EDT 2007
// $Id: CaloDisplay.cc,v 1.3 2010/06/03 09:11:03 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TString.h>

using namespace std;

//
// class decleration
//

class CaloDisplay : public edm::EDAnalyzer {
   public:
      explicit CaloDisplay(const edm::ParameterSet&);
      ~CaloDisplay();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      double       computeEt(const DetId &id, double energy);

      // ----------member data ---------------------------

   const CaloGeometry *geo;

   edm::Service<TFileService> fs;

   TNtuple* ntjet;
   TNtuple* ntcalo;
   TNtuple* ntmc;

   TH2D* hHad;
   TH2D* hEM;
   TH2D* hMC;

   int NoE;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CaloDisplay::CaloDisplay(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


CaloDisplay::~CaloDisplay()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CaloDisplay::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   NoE += 1;

   Handle<EcalRecHitCollection> ebhits;
   iEvent.getByLabel(InputTag("ecalRecHit","EcalRecHitsEB"),ebhits);
   for( size_t ihit = 0; ihit<ebhits->size(); ++ihit){
      const EcalRecHit & rechit = (*ebhits)[ihit];
      const DetId& detId = rechit.id();
      const GlobalPoint& hitpoint = geo->getPosition(detId);
      double phi = hitpoint.phi();
      double eta = hitpoint.eta();
      double et = computeEt(detId, rechit.energy());
      hEM->Fill(eta,phi,et); 
      ntcalo->Fill(eta,phi,et,2);
  }

   Handle<EcalRecHitCollection> eehits;
   iEvent.getByLabel(InputTag("ecalRecHit","EcalRecHitsEE"),eehits);
   for( size_t eit = 0; eit<eehits->size(); ++eit){
      const EcalRecHit & rechit = (*eehits)[eit];
      const DetId& detId = rechit.id();
      const GlobalPoint& hitpoint = geo->getPosition(detId);
      double phi = hitpoint.phi();
      double eta = hitpoint.eta();
      double et = computeEt(detId, rechit.energy());
      hEM->Fill(eta,phi,et);
      ntcalo->Fill(eta,phi,et,2);

   }
   
   Handle<HBHERecHitCollection> hhits;
   iEvent.getByLabel(InputTag("hbhereco"),hhits);
   for( size_t hhit = 0; hhit<hhits->size(); ++hhit){
      const HBHERecHit & rechit = (*hhits)[hhit];
      const DetId& detId = rechit.id();
      const GlobalPoint& hitpoint = geo->getPosition(detId);
      double phi = hitpoint.phi();
      double eta = hitpoint.eta();
      double et = computeEt(detId, rechit.energy());
      hHad->Fill(eta,phi,et);
      ntcalo->Fill(eta,phi,et,1);

   }
   
   Handle<HFRecHitCollection> hfhits;
   iEvent.getByLabel(InputTag("hfreco"),hfhits);
   for( size_t ihf = 0; ihf<hfhits->size(); ++ihf){
      const HFRecHit & rechit = (*hfhits)[ihf];
      const DetId& detId = rechit.id();
      const GlobalPoint& hitpoint = geo->getPosition(detId);
      double phi = hitpoint.phi();
      double eta = hitpoint.eta();
      double et = computeEt(detId, rechit.energy());
      hHad->Fill(eta,phi,et);
      ntcalo->Fill(eta,phi,et,1);
   }
   
   Handle<HORecHitCollection> hohits;
   iEvent.getByLabel(InputTag("horeco"),hohits);
   for( size_t iho = 0; iho<hohits->size(); ++iho){
      const HORecHit & rechit = (*hohits)[iho];
      const DetId& detId = rechit.id();
      const GlobalPoint& hitpoint = geo->getPosition(detId);
      double phi = hitpoint.phi();
      double eta = hitpoint.eta();
      double et = computeEt(detId, rechit.energy());
      hHad->Fill(eta,phi,et);
      ntcalo->Fill(eta,phi,et,1);
   }

   edm::Handle<reco::JetView> jets;
   iEvent.getByLabel("iterativeConePu5CaloJets",jets);
   for(unsigned int ijet = 0; ijet < jets->size(); ++ijet){
      const reco::Jet* jet = &((*jets)[ijet]); 
      double pt = jet->pt();
      double eta = jet->eta();
      double phi = jet->phi();
      ntjet->Fill(eta,phi,pt);
   }

   Handle<HepMCProduct> mc;
   iEvent.getByLabel("signal",mc);
   const HepMC::GenEvent* evt = mc->GetEvent();

   int all = evt->particles_size();
   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){

	 int pdg = (*it)->pdg_id();
	 float eta = (*it)->momentum().eta();
	 float phi = (*it)->momentum().phi();
	 float pt = (*it)->momentum().perp();
	 int st = (*it)->status();

	 if(st == 1) hMC->Fill(eta,phi,pt);
	 ntmc->Fill(eta,phi,pt,pdg);
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
CaloDisplay::beginJob(const edm::EventSetup& iSetup)
{


   hHad = fs->make<TH2D>("hHad","E_{T} in Hcal;#eta;#phi",60,-2.5,2.5,75,-3.15,3.15);
   hEM = fs->make<TH2D>("hEM","E_{T} in Ecal;#eta;#phi",150,-2.5,2.5,180,-3.15,3.15);
   hMC = fs->make<TH2D>("hMC","E_{T} from MC;#eta;#phi",150,-2.5,2.5,180,-3.15,3.15);

   ntcalo = fs->make<TNtuple>("ntcalo","","eta:phi:et:type");
   ntjet = fs->make<TNtuple>("ntjet","","eta:phi:et");
   ntmc = fs->make<TNtuple>("ntmc","","eta:phi:pt:pdg:st");

   NoE = 0;
   edm::ESHandle<CaloGeometry> pGeo;
   iSetup.get<CaloGeometryRecord>().get(pGeo);
   geo = pGeo.product();


}

// ------------ method called once each job just after ending the event loop  ------------
void 
CaloDisplay::endJob() {

   double n = 1/NoE;

}

double CaloDisplay::computeEt(const DetId &id, double energy){
   const GlobalPoint& pos=geo->getPosition(id);
   double et = energy*sin(pos.theta());
   return et;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloDisplay);
