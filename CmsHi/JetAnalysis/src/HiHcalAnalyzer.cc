// -*- C++ -*-
//
// Package:    HiHcalAnalyzer
// Class:      HiHcalAnalyzer
// 
/**\class HiHcalAnalyzer HiHcalAnalyzer.cc CmsHi/HiHcalAnalyzer/src/HiHcalAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Thu Dec  1 10:28:28 EST 2011
// $Id: HiHcalAnalyzer.cc,v 1.2 2011/12/02 01:05:57 yjlee Exp $
//
//

#define versionTag "v1"

// system include files
#include <memory>
#include <vector>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TNtuple.h"

using namespace std;


//
// class declaration
//

class HiHcalAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HiHcalAnalyzer(const edm::ParameterSet&);
      ~HiHcalAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

   edm::Service<TFileService> fs;
   TTree* t;


   int filterstatus, noisetype;
   float emenergy, hadenergy, trackenergy;
   float min10, max10, rms10;
   float min25, max25, rms25;
   int cnthit10, cnthit25;
   float mine2ts, mine10ts;
   float maxe2ts, maxe10ts;
   int maxzeros;
   int maxhpdhits, maxhpdhitsnoother, maxrbxhits;
   float minhpdemf, minrbxemf;
   int nproblemRBXs;
   int nisolnoise;
   float isolnoisee, isolnoiseet;
   int nflatnoise;
   float flatnoisee, flatnoiseet;
   int nspikenoise;
   float spikenoisee, spikenoiseet;
   int ntrianglenoise;
   float trianglenoisee, trianglenoiseet;
   int nts4ts5noise;
   float ts4ts5noisee, ts4ts5noiseet;

   bool hasBadRBXTS4TS5;


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
HiHcalAnalyzer::HiHcalAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


HiHcalAnalyzer::~HiHcalAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HiHcalAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<HcalNoiseSummary> hcalSummary;
   iEvent.getByType(hcalSummary);

   filterstatus= hcalSummary->noiseFilterStatus ();
   noisetype= hcalSummary->noiseType ();
   emenergy= hcalSummary->eventEMEnergy ();
   hadenergy= hcalSummary->eventHadEnergy ();
   trackenergy= hcalSummary-> eventTrackEnergy();
   min10= hcalSummary->min10GeVHitTime ();
   max10= hcalSummary->max10GeVHitTime ();
   rms10= hcalSummary->rms10GeVHitTime ();
   min25= hcalSummary->min25GeVHitTime ();
   max25= hcalSummary->max25GeVHitTime ();
   rms25= hcalSummary->rms25GeVHitTime ();
   cnthit10= hcalSummary->num10GeVHits ();
   cnthit25= hcalSummary->num25GeVHits ();
   mine2ts= hcalSummary->minE2TS ();
   mine10ts= hcalSummary->minE10TS ();
   maxe2ts= hcalSummary->maxE2TS ();
   maxe10ts= hcalSummary->maxE10TS ();
   maxzeros= hcalSummary->maxZeros ();
   maxhpdhits= hcalSummary->maxHPDHits ();
   maxhpdhitsnoother= hcalSummary->maxHPDNoOtherHits ();
   maxrbxhits= hcalSummary->maxRBXHits ();
   minhpdemf= hcalSummary->minHPDEMF();
   minrbxemf= hcalSummary->minRBXEMF ();

   nproblemRBXs= hcalSummary->numProblematicRBXs ();
   nisolnoise= hcalSummary->numIsolatedNoiseChannels ();
   isolnoisee= hcalSummary->isolatedNoiseSumE ();
   isolnoiseet= hcalSummary->isolatedNoiseSumEt ();
   nflatnoise= hcalSummary->numFlatNoiseChannels ();
   flatnoisee= hcalSummary->flatNoiseSumE ();
   flatnoiseet= hcalSummary->flatNoiseSumEt ();
   nspikenoise= hcalSummary->numSpikeNoiseChannels ();
   spikenoisee= hcalSummary->spikeNoiseSumE ();
   spikenoiseet= hcalSummary->spikeNoiseSumEt ();
   ntrianglenoise= hcalSummary->numTriangleNoiseChannels ();
   trianglenoisee= hcalSummary->triangleNoiseSumE ();
   trianglenoiseet= hcalSummary->triangleNoiseSumEt ();
   nts4ts5noise= hcalSummary->numTS4TS5NoiseChannels ();
   ts4ts5noisee= hcalSummary->TS4TS5NoiseSumE ();
   ts4ts5noiseet= hcalSummary->TS4TS5NoiseSumEt ();
   hasBadRBXTS4TS5= hcalSummary->HasBadRBXTS4TS5 ();

   t->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HiHcalAnalyzer::beginJob()
{

   t = fs->make<TTree>("hcalNoise",versionTag);
   t->Branch("filterstatus",&filterstatus,"filterstatus/I");
   t->Branch("noisetype",&noisetype,"noisetype/I");
   t->Branch("emenergy",&emenergy,"emenergy/F");
   t->Branch("hadenergy",&hadenergy,"hadenergy/F");
   t->Branch("trackenergy",&trackenergy,"trackenergy/F");
   t->Branch("min10",&min10,"min10/F");
   t->Branch("max10",&max10,"max10/F");
   t->Branch("rms10",&rms10,"rms10/F");

   t->Branch("min25",&min25,"min25/F");
   t->Branch("max25",&max25,"max25/F");
   t->Branch("rms25",&rms25,"rms25/F");

   t->Branch("cnthit10",&cnthit10,"cnthit10/I");
   t->Branch("cnthit25",&cnthit25,"cnthit25/I");

   t->Branch("mine2ts",&mine2ts,"mine2ts/F");

   t->Branch("maxe2ts",&maxe2ts,"maxe2ts/F");
   t->Branch("maxe10ts",&maxe10ts,"maxe10ts/F");


   t->Branch("maxzeros",&maxzeros,"maxzeros/I");
   t->Branch("maxhpdhits",&maxhpdhits,"maxhpdhits/I");
   t->Branch("maxhpdhitsnoother",&maxhpdhitsnoother,"maxhpdhitsnoother/I");
   t->Branch("maxrbxhits",&maxrbxhits,"maxrbxhits/I");

   t->Branch("minhpdemf",&minhpdemf,"minhpdemf/F");
   t->Branch("minrbxemf",&minrbxemf,"minrbxemf/F");
   t->Branch("nproblemRBXs",&nproblemRBXs,"nproblemRBXs/I");
   t->Branch("nisolnoise",&nisolnoise,"nisolnoise/I");
   t->Branch("isolnoisee",&isolnoisee,"isolnoisee/F");

   t->Branch("nflatnoise",&nflatnoise,"nflatnoise/I");
   t->Branch("flatnoisee",&flatnoisee,"flatnoisee/F");
   t->Branch("flatnoiseet",&flatnoiseet,"flatnoiseet/F");
   t->Branch("nspikenoise",&nspikenoise,"nspikenoise/I");
   t->Branch("spikenoisee",&spikenoisee,"spikenoisee/F");
   t->Branch("spikenoiseet",&spikenoiseet,"spikenoiseet/F");
   t->Branch("ntrianglenoise",&ntrianglenoise,"ntrianglenoise/I");
   t->Branch("trianglenoisee",&trianglenoisee,"trianglenoisee/F");
   t->Branch("trianglenoiseet",&trianglenoiseet,"trianglenoiseet/F");
   t->Branch("nts4ts5noise",&nts4ts5noise,"nts4ts5noise/I");
   t->Branch("ts4ts5noisee",&ts4ts5noisee,"ts4ts5noisee/F");
   t->Branch("ts4ts5noiseet",&ts4ts5noiseet,"ts4ts5noiseet/F");

   t->Branch("hasBadRBXTS4TS5",&hasBadRBXTS4TS5,"hasBadRBXTS4TS5/O");



}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiHcalAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
HiHcalAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HiHcalAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HiHcalAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HiHcalAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HiHcalAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiHcalAnalyzer);
