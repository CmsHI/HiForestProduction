// -*- C++ -*-              
// Package:    Ate
// Class:      Ate 
//                                                                                                                                                                       
// Package:    PatAlgos
//class HiPatPhotonAna HiPatPhotonAna.cc PhysicsTools/PatAlgos/test/HiPatPhotonAna.cc
//   Description: <one line class summary>                                                                                                                                   
// Implementation:                                                                                                                                                         
// <Notes on implementation>                                                                                                                                           
//                                                                                                                                                                
// Original Author:  Yongsun Kim                                                                                                                                    
//         Created:  Sun Sep 20 05:04:32 EDT 2009                                                                                                                   
// $Id: HiPatPhotonAna.cc,v 1.3 2010/11/05 08:38:11 kimy Exp $          
//                                                                                                                                                                       
//                                                                                                                                                                       

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/ESHandle.h"



#include <string>
#include "TNtuple.h"
#include "TTree.h"
//#include <vector>
#include "CmsHi/PhotonAnalysis/interface/HiMCGammaJetSignalDef.h"
#include "CmsHi/PhotonAnalysis/interface/HiPhotonMCTruth.h"

namespace edm { using ::std::advance; }

//
// class decleration
//

class HiPatPhotonAna : public edm::EDAnalyzer {
   public:
      explicit HiPatPhotonAna(const edm::ParameterSet&);
      ~HiPatPhotonAna();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
   
      // ----------member data ---------------------------
  edm::InputTag photons_;
  edm::InputTag genParticles_;

  std::string   label_;
      enum TestMode { TestRead, TestWrite, TestExternal };
      TestMode mode_;
      TNtuple *datatemp;
      edm::Service<TFileService> fs;
   
   TTree *tr;
   
   int npho;
   float e[50];
   float rawe[50];
   float et[50];
   float eta[50];
   float phi[50];
   float r9[50];
   float r1x5[50];
   float r2x5[50];
   float e9[50];
   float e25[50];
   float maxEXtal[50];
   float sEE[50];
   float sII[50];

   float cc1[50];
   float cc2[50];
   float cc3[50];
   float cc4[50];
   float cc5[50];
   float cr1[50];
   float cr2[50];
   float cr3[50];
   float cr4[50];
   float cr5[50];

   float t11[50];
   float t12[50];
   float t13[50];
   float t14[50];
   float t21[50];
   float t22[50];
   float t23[50];
   float t24[50];
   float t31[50];
   float t32[50];
   float t33[50];
   float t34[50];
   float t41[50];
   float t42[50];
   float t43[50];
   float t44[50];

   float dr11[50];
   float dr12[50];
   float dr13[50];
   float dr14[50];
   float dr21[50];
   float dr22[50];
   float dr23[50];
   float dr24[50];
   float dr31[50];
   float dr32[50];
   float dr33[50];
   float dr34[50];
   float dr41[50];
   float dr42[50];
   float dr43[50];
   float dr44[50];
   
   int mcMatched[50];
   float ge[50];
   float get[50];
   float geta[50];
   float gphi[50];
   int iso[50];
   int dir[50];
   int mid[50];

};

//
// constructors and destructor
//
HiPatPhotonAna::HiPatPhotonAna(const edm::ParameterSet& iConfig):
   photons_(iConfig.getParameter<edm::InputTag>("photons")),
   genParticles_(iConfig.getParameter<edm::InputTag>("genParticles"))
  
{
}


HiPatPhotonAna::~HiPatPhotonAna()
{
}

// ------------ method called to for each event  ------------
void
HiPatPhotonAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   
   edm::Handle<reco::GenParticleCollection> inputHandle;
   iEvent.getByLabel(genParticles_,inputHandle);
   const GenParticleCollection *collection1 = inputHandle.product();
   // The photon signal defenition.
   HiPhotonMCTruth photonTruth(inputHandle);                                                                                  

    edm::Handle<edm::View<pat::Photon> > photons;
   iEvent.getByLabel(photons_,photons);

   std::auto_ptr<std::vector<pat::Photon> > output(new std::vector<pat::Photon>());
   
   npho = 0;
   for (edm::View<pat::Photon>::const_iterator photon = photons->begin(), end = photons->end(); photon != end; ++photon) {
      npho++;
      if ( npho >= 50 ) break;
     int  i = npho-1;
      e[i]     = photon->energy();
      rawe[i]  = photon->superCluster()->rawEnergy();
      et[i]    = photon->et();
      eta[i]   = photon->eta();
      phi[i]   = photon->phi();
      r9[i]    = photon->r9();
      r1x5[i]  = photon->r1x5();
      r2x5[i]  = photon->r2x5();
      e9[i]    = photon->e3x3();
      e25[i]   = photon->e5x5();
      maxEXtal[i]=photon->maxEnergyXtal();
      sEE[i]   = photon->sigmaEtaEta();
      sII[i]   = photon->sigmaIetaIeta();

      cc1[i]    = photon->userFloat("isoCC1");
      cc2[i]    = photon->userFloat("isoCC2");
      cc3[i]    = photon->userFloat("isoCC3");
      cc4[i]    = photon->userFloat("isoCC4");
      cc5[i]    = photon->userFloat("isoCC5");
      cr1[i]    = photon->userFloat("isoCR1");
      cr2[i]    = photon->userFloat("isoCR2");
      cr3[i]    = photon->userFloat("isoCR3");
      cr4[i]    = photon->userFloat("isoCR4");
      cr5[i]    = photon->userFloat("isoCR5");
      t11[i]    = photon->userFloat("isoT11");
      t12[i]    = photon->userFloat("isoT12");
      t13[i]    = photon->userFloat("isoT13");
      t14[i]    = photon->userFloat("isoT14");
      t21[i]    = photon->userFloat("isoT21");
      t22[i]    = photon->userFloat("isoT22");
      t23[i]    = photon->userFloat("isoT23");
      t24[i]    = photon->userFloat("isoT24");
      t31[i]    = photon->userFloat("isoT31");
      t32[i]    = photon->userFloat("isoT32");
      t33[i]    = photon->userFloat("isoT33");
      t34[i]    = photon->userFloat("isoT34");
      t41[i]    = photon->userFloat("isoT41");
      t42[i]    = photon->userFloat("isoT42");
      t43[i]    = photon->userFloat("isoT43");
      t44[i]    = photon->userFloat("isoT44");
     
      dr11[i]    = photon->userFloat("isoDR11");
      dr12[i]    = photon->userFloat("isoDR12");
      dr13[i]    = photon->userFloat("isoDR13");
      dr14[i]    = photon->userFloat("isoDR14");
      dr21[i]    = photon->userFloat("isoDR21");
      dr22[i]    = photon->userFloat("isoDR22");
      dr23[i]    = photon->userFloat("isoDR23");
      dr24[i]    = photon->userFloat("isoDR24");
      dr31[i]    = photon->userFloat("isoDR31");
      dr32[i]    = photon->userFloat("isoDR32");
      dr33[i]    = photon->userFloat("isoDR33");
      dr34[i]    = photon->userFloat("isoDR34");
      dr41[i]    = photon->userFloat("isoDR41");
      dr42[i]    = photon->userFloat("isoDR42");
      dr43[i]    = photon->userFloat("isoDR43");
      dr44[i]    = photon->userFloat("isoDR44");
      
      
      const reco::GenParticle *genPhoton = photon->genParticle();
      mcMatched[i] = 0 ;
      ge[i]   =0;
      get[i]  =0;
      geta[i] =0;
      gphi[i] =0;
      iso[i]  =0;
      dir[i]  =0;
      mid[i]  =0;

      if ( genPhoton != 0){
	 ge[i]   = genPhoton->energy();
	 get[i]  = genPhoton->et();
	 geta[i] = genPhoton->eta();
	 gphi[i] = genPhoton->phi();
	 if ( genPhoton->mother() ==0 ) 
	    mid[i] = genPhoton->mother()->pdgId();
	 if ( photonTruth.IsIsolatedPP(*genPhoton) == true ) 
	    iso[i] = 1 ;
      }
   }  
   
   tr->Fill();
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiPatPhotonAna::beginJob(const edm::EventSetup&)
{

   tr = new TTree("RecoPhotonTree","recoPhoTree");
   // RECO photon
   tr->Branch("npho",&npho,"npho/I");
   tr->Branch("e",&e,"e[npho]/F");
   tr->Branch("rawe",&rawe,"rawe[npho]/F");
   tr->Branch("et",&et,"et[npho]/F");
   tr->Branch("eta",&eta,"eta[npho]/F");
   tr->Branch("phi",&phi,"phi[npho]/F");
   
   // supercluster shape
   tr->Branch("r9",&r9,"r9[npho]/F");
   tr->Branch("r1x5",&r1x5,"r1x5[npho]/F");
   tr->Branch("r2x5",&r2x5,"r2x5[npho]/F");
   tr->Branch("e9",&e9,"e9[npho]/F");
   tr->Branch("e25",&e25,"e25[npho]/F");
   tr->Branch("maxEXtal",&maxEXtal,"maxEXtal[npho]/F");
   tr->Branch("sEE",&sEE,"sEE[npho]/F");
   tr->Branch("sII",&sII,"sII[npho]/F");
      
   // cone isolation variable
   tr->Branch("cc1",&cc1,"cc1[npho]/F");
   tr->Branch("cc2",&cc2,"cc2[npho]/F");
   tr->Branch("cc3",&cc3,"cc3[npho]/F");
   tr->Branch("cc4",&cc4,"cc4[npho]/F");
   tr->Branch("cc5",&cc5,"cc5[npho]/F");
   tr->Branch("cr1",&cr1,"cr1[npho]/F");
   tr->Branch("cr2",&cr2,"cr2[npho]/F");
   tr->Branch("cr3",&cr3,"cr3[npho]/F");
   tr->Branch("cr4",&cr4,"cr4[npho]/F");
   tr->Branch("cr5",&cr5,"cr5[npho]/F");
   tr->Branch("t11",&t11,"t11[npho]/F");
   tr->Branch("t12",&t12,"t12[npho]/F");
   tr->Branch("t13",&t13,"t13[npho]/F");
   tr->Branch("t14",&t14,"t14[npho]/F");
   tr->Branch("t21",&t21,"t21[npho]/F");
   tr->Branch("t22",&t22,"t22[npho]/F");
   tr->Branch("t23",&t23,"t23[npho]/F");
   tr->Branch("t24",&t24,"t24[npho]/F");
   tr->Branch("t31",&t31,"t31[npho]/F");
   tr->Branch("t32",&t32,"t32[npho]/F");
   tr->Branch("t33",&t33,"t33[npho]/F");
   tr->Branch("t34",&t34,"t34[npho]/F");
   tr->Branch("t41",&t41,"t41[npho]/F");
   tr->Branch("t42",&t42,"t42[npho]/F");
   tr->Branch("t43",&t43,"t43[npho]/F");
   tr->Branch("t44",&t44,"t44[npho]/F");
   tr->Branch("dr11",&dr11,"dr11[npho]/F");
   tr->Branch("dr12",&dr12,"dr12[npho]/F");
   tr->Branch("dr13",&dr13,"dr13[npho]/F");
   tr->Branch("dr14",&dr14,"dr14[npho]/F");
   tr->Branch("dr21",&dr21,"dr21[npho]/F");
   tr->Branch("dr22",&dr22,"dr22[npho]/F");
   tr->Branch("dr23",&dr23,"dr23[npho]/F");
   tr->Branch("dr24",&dr24,"dr24[npho]/F");
   tr->Branch("dr31",&dr31,"dr31[npho]/F");
   tr->Branch("dr32",&dr32,"dr32[npho]/F");
   tr->Branch("dr33",&dr33,"dr33[npho]/F");
   tr->Branch("dr34",&dr34,"dr34[npho]/F");
   tr->Branch("dr41",&dr41,"dr41[npho]/F");
   tr->Branch("dr42",&dr42,"dr42[npho]/F");
   tr->Branch("dr43",&dr43,"dr43[npho]/F");
   tr->Branch("dr44",&dr44,"dr44[npho]/F");

   // genParticle info
   tr->Branch("mcMatched",&mcMatched,"mcMatched[npho]/I");
   tr->Branch("ge",&ge,"ge[npho]/F");
   tr->Branch("get",&get,"get[npho]/F");
   tr->Branch("geta",&geta,"geta[npho]/F");
   tr->Branch("gphi",&gphi,"gphi[npho]/F");
   tr->Branch("iso",&iso,"iso[npho]/I");
   tr->Branch("dir",&dir,"dir[npho]/I");
   tr->Branch("mid",&mid,"mid[npho]/I");

   }

// ------------ method called once each job just after ending the event loop  ------------
void 
HiPatPhotonAna::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiPatPhotonAna);
