#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include <math.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/Vertex.h>

int main(int argc, char* argv[]) 
{

  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  fwlite::TFileService fs = fwlite::TFileService("pixel.root");
  TFileDirectory theDir = fs.mkdir("pixel");
  TH1F* pixelpt_  = theDir.make<TH1F>("pixelpt", "pixelpt",    200,  0.,10.);
  TH1F* pixelphi_ = theDir.make<TH1F>("pixelphi", "pixelphi", 200, -3.15, 3.15);
  TH1F* pixeleta_ = theDir.make<TH1F>("pixeleta", "pixeleta", 200, -2.5, 2.5);
  TH1F* pixeld0_ = theDir.make<TH1F>("pixeld0", "pixeld0", 200, -0.25, 0.25);
  TH1F* pixeld0err_ = theDir.make<TH1F>("pixeld0err", "pixeld0err", 200, 0, 5);
  TH1F* pixeldz_ = theDir.make<TH1F>("pixeldz", "pixeldz", 200, -0.25, 0.25);
  TH1F* pixeldzerr_ = theDir.make<TH1F>("pixeldzerr", "pixeldzerr", 200, 0, 5);

//  TFile* inFile = TFile::Open( "file:reco.root" );
//  TFile* inFile = TFile::Open( "file:tier0.root" );
//  TFile* inFile = TFile::Open( "file:/scratch/appelte1/a_few_more_pixel_tracks.root" );


  std::vector<std::string> files;

  files.push_back( "file:/scratch/appelte1/a_few_organ_pixel_tracks.root");

  fwlite::ChainEvent ev(files);

  unsigned int iEvent=0;
  for(ev.toBegin(); !ev.atEnd(); ++ev, ++iEvent){

    if ( iEvent >= 200 ) continue;

    std::cout << "Processing event " << iEvent << "\n" ;

    edm::EventBase const & event = ev;

    edm::Handle<std::vector<reco::Track> > tracks;
    edm::InputTag cmnLabel("hiLowPtPixelTracks");
    event.getByLabel(cmnLabel, tracks);

    edm::Handle<std::vector<reco::Vertex> > vertex;
    edm::InputTag vLabel("hiSelectedVertex");
    event.getByLabel(vLabel, vertex);


    math::XYZPoint vtxPoint(0.0,0.0,0.0);
    double vzErr =0.0, vxErr=0.0, vyErr=0.0;
  
    if(vertex->size()>0) {
      vtxPoint=vertex->begin()->position();
      vzErr=vertex->begin()->zError();
      vxErr=vertex->begin()->xError();
      vyErr=vertex->begin()->yError();
    }


    std::vector<reco::Track>::const_iterator tk = tracks->begin();
    for ( ; tk != tracks->end(); ++tk )
    {
      if ( fabs( tk->eta() ) <= 1.0 )
      {
         pixelpt_->Fill( tk->pt() );
         pixeleta_->Fill( tk->eta() );
         pixelphi_->Fill( tk->phi() );

         double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
         d0 = -1.*tk->dxy(vtxPoint);
         dz = tk->dz(vtxPoint);
         d0sigma = sqrt(tk->d0Error()*tk->d0Error()+vxErr*vyErr);
         dzsigma = sqrt(tk->dzError()*tk->dzError()+vzErr*vzErr);
   
         pixeld0_->Fill( d0 );
         pixeldz_->Fill( dz );
         pixeldzerr_->Fill( fabs(dz/dzsigma) );
         pixeld0err_->Fill( fabs(d0/d0sigma) ); 
      }


    }
  }


  return 0;

}


