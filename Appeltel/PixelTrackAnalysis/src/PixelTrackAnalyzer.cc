// -*- C++ -*-
//
// Package:    PixelTrackAnalysis
// Class:      PixelTrackAnalysis
// 
/**\class PixelTrackAnalysis PixelTrackAnalysis.cc Appeltel/PixelTrackAnalysis/src/PixelTrackAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric A. Appelt
//         Created:  Fri Nov 12 04:59:50 EST 2010
// $Id: PixelTrackAnalyzer.cc,v 1.5 2010/11/30 21:29:13 appeltel Exp $
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

#include <TH1.h>
#include <TH2.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"


//
// class declaration
//

class PixelTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PixelTrackAnalyzer(const edm::ParameterSet&);
      ~PixelTrackAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag pixelSrc_;
      edm::InputTag vertexSrc_;
      edm::InputTag genSrc_;

      bool doGenerator_;

      double etaCut_;
      double ptMin_;

      CentralityProvider * centrality_;



      TH1F* pixelpt_; 
      TH1F* pixelphi_;
      TH1F* pixeleta_;
      TH1F* pixeld0_; 
      TH1F* pixeld0err_; 
      TH1F* pixeldz_; 
      TH1F* pixeldzerr_; 
      TH1F* pixelchi2_;
      TH1F* pixelnhit_;
 
      TH1F* pixelptcent_[10];
      TH1F* pixelphicent_[10];
      TH1F* pixeletacent_[10];
      TH1F* pixeld0cent_[10]; 
      TH1F* pixeld0errcent_[10]; 
      TH1F* pixeldzcent_[10]; 
      TH1F* pixeldzerrcent_[10]; 
      TH1F* pixelchi2cent_[10];
      TH1F* pixelnhitcent_[10];


      TH1F* genptcent_[10];
      TH1F* genphicent_[10];
      TH1F* genetacent_[10];
      TH1I* centbins_;

};

PixelTrackAnalyzer::PixelTrackAnalyzer(const edm::ParameterSet& iConfig):
pixelSrc_(iConfig.getParameter<edm::InputTag>("pixelSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
doGenerator_(iConfig.getParameter<bool>("doGenerator")),
etaCut_(iConfig.getParameter<double>("etaCut")),
ptMin_(iConfig.getParameter<double>("ptMin"))
{

  edm::Service<TFileService> fs;

  char centStrings[10][256] = { "0-10%", "10-20%", "20-30%", "30-40%",
       "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };

  pixelpt_  = fs->make<TH1F>("pixelpt", "Pixel Track p_{T} Spectrum",    200,  0.,10.);
  pixelphi_ = fs->make<TH1F>("pixelphi", "Pixel Track #phi Distribution", 200, -3.15, 3.15);
  pixeleta_ = fs->make<TH1F>("pixeleta", "Pixel Track #eta Distribution", 200, -2.5, 2.5);
  pixeld0_ = fs->make<TH1F>("pixeld0", "Pixel Track Transverse DCA Distribution", 200, -0.25, 0.25);
  pixeld0err_ = fs->make<TH1F>("pixeld0err", "Pixel Track Transverse DCA Significance", 200, 0, 5);
  pixeldz_ = fs->make<TH1F>("pixeldz", "Pixel Track Longitudinal DCA Distribution", 200, -0.25, 0.25);
  pixeldzerr_ = fs->make<TH1F>("pixeldzerr", "Pixel Track Longitudinal DCA Significance", 200, 0, 5);
  pixelchi2_ = fs->make<TH1F>("pixelchi2", "Pixel Track #chi^2 / n.d.o.f. Distribution",200,0.,100.);
  pixelnhit_ = fs->make<TH1F>("pixelnhit", "Pixel Track Valid Hit Distribution",30.,0.,30.);

  for( int i = 0; i<10; i++)
  {
    pixelptcent_[i] = fs->make<TH1F>(Form("pixelptcent%d",i),
        Form("Pixel Track p_{T} Spectrum Centrality Range %s",centStrings[i]), 200, 0., 10. );
    pixelphicent_[i] = fs->make<TH1F>(Form("pixelphicent%d",i),
        Form("Pixel Track #phi Distribution Centrality Range %s",centStrings[i]), 200, -3.15, 3.15 );
    pixeletacent_[i] = fs->make<TH1F>(Form("pixeletacent%d",i),
        Form("Pixel Track #eta Distribution Centrality Range %s",centStrings[i]), 200, -2.5, 2.5 );
    pixeld0cent_[i] = fs->make<TH1F>(Form("pixeld0cent%d",i),
        Form("Pixel Track Transverse DCA  Distribution Centrality Range %s",centStrings[i]), 200, -0.25, 0.25 );
    pixeld0errcent_[i] = fs->make<TH1F>(Form("pixeld0errcent%d",i),
        Form("Pixel Track Transverse DCA  Significance Centrality Range %s",centStrings[i]), 200, 0, 5 );
    pixeldzcent_[i] = fs->make<TH1F>(Form("pixeldzcent%d",i),
        Form("Pixel Track Longitudinal DCA  Distribution Centrality Range %s",centStrings[i]), 200, -0.25, 0.25 );
    pixeldzerrcent_[i] = fs->make<TH1F>(Form("pixeldzerrcent%d",i),
        Form("Pixel Track Longitudinal DCA  Significance Centrality Range %s",centStrings[i]), 200, 0, 5 );
    pixelchi2cent_[i] = fs->make<TH1F>(Form("pixelchi2cent%d",i), 
        Form("Pixel Track #chi^2 / n.d.o.f. Distribution Centrality Range %s",centStrings[i]),200,0.,100.);
    pixelnhitcent_[i] = fs->make<TH1F>(Form("pixelnhitcent%d",i), 
        Form("Pixel Track Valid Hit Distribution Centrality Range %s",centStrings[i]),30.,0.,30.);
    genptcent_[i] = fs->make<TH1F>(Form("genptcent%d",i),
        Form("Charged GenParticle p_{T} Spectrum Centrality Range %s",centStrings[i]), 200, 0., 10. );
    genphicent_[i] = fs->make<TH1F>(Form("genphicent%d",i),
        Form("Charged GenParticle #phi Distribution Centrality Range %s",centStrings[i]), 200, -3.15, 3.15 );
    genetacent_[i] = fs->make<TH1F>(Form("genetacent%d",i),
        Form("Charged GenParticle #eta Distribution Centrality Range %s",centStrings[i]), 200, -2.5, 2.5 );

  }
  centbins_ = fs->make<TH1I>("centbins","Centrality Distribution", 40, 0, 39);

  // safety
  centrality_ = 0;
}


PixelTrackAnalyzer::~PixelTrackAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PixelTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<reco::Track> > tracks;
   iEvent.getByLabel(pixelSrc_, tracks);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   
   // Get Centrality bin for the event

   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
   int bin = centrality_->getBin();

   centbins_->Fill(bin);

   // find the vertex point and error

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
  
   if(vertex->size()>0) {
     vtxPoint=vertex->begin()->position();
     vzErr=vertex->begin()->zError();
     vxErr=vertex->begin()->xError();
     vyErr=vertex->begin()->yError();
   }


   // loop over tracks, fill histos

   std::vector<reco::Track>::const_iterator tk = tracks->begin();
   for ( ; tk != tracks->end(); ++tk )
   {
     if ( fabs( tk->eta() ) <= etaCut_ && tk->pt() >= ptMin_ )
     {
       pixelpt_->Fill( tk->pt() );
       pixeleta_->Fill( tk->eta() );
       pixelphi_->Fill( tk->phi() );
       pixelchi2_->Fill( tk->normalizedChi2() );
       pixelnhit_->Fill( tk->numberOfValidHits() );

       pixelptcent_[bin/4]->Fill( tk->pt() );
       pixeletacent_[bin/4]->Fill( tk->eta() );
       pixelphicent_[bin/4]->Fill( tk->phi() );
       pixelchi2cent_[bin/4]->Fill( tk->normalizedChi2() );
       pixelnhitcent_[bin/4]->Fill( tk->numberOfValidHits() );

       double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
       d0 = -1.*tk->dxy(vtxPoint);
       dz = tk->dz(vtxPoint);
       d0sigma = sqrt(tk->d0Error()*tk->d0Error()+vxErr*vyErr);
       dzsigma = sqrt(tk->dzError()*tk->dzError()+vzErr*vzErr);
   
       pixeld0_->Fill( d0 );
       pixeldz_->Fill( dz );
       pixeldzerr_->Fill( fabs(dz/dzsigma) );
       pixeld0err_->Fill( fabs(d0/d0sigma) ); 

       pixeld0cent_[bin/4]->Fill( d0 );
       pixeldzcent_[bin/4]->Fill( dz );
       pixeldzerrcent_[bin/4]->Fill( fabs(dz/dzsigma) );
       pixeld0errcent_[bin/4]->Fill( fabs(d0/d0sigma) ); 

       pixelptcent_[bin/4]->Fill ( tk->pt() );

     }
  }

  if ( doGenerator_ )
  {

    Handle<reco::GenParticleCollection > gens;
    iEvent.getByLabel(genSrc_, gens);

    reco::GenParticleCollection::const_iterator q;
    for( q = gens->begin(); q != gens->end(); ++q )
    {
       if ( q->charge() != 0 && q->status() == 1 )
       {
         genptcent_[bin/4]->Fill( q->pt() );
         genetacent_[bin/4]->Fill( q->eta() );
         genphicent_[bin/4]->Fill( q->phi() );
       }

    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
PixelTrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelTrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelTrackAnalyzer);

