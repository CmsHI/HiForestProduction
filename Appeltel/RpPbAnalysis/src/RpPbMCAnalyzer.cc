//
// Original Author:  Eric Appelt
//         Created:  Fri Nov 19 16:55:16 CST 2010
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


class RpPbMCAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RpPbMCAnalyzer(const edm::ParameterSet&);
      ~RpPbMCAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void printWeirdParticle( const reco::GenParticle &);
      void printRoguePion( const reco::GenParticle &, double vdist);

      // ----------member data ---------------------------

      std::map<std::string,TH1F*> spectra_;

      TH1F* events_;

      TH1F* vtxNtrk_;
      TH1F* vtxN_;
      TH1F* distPV_;

      int nevt_;

      edm::InputTag genSrc_;
      edm::InputTag vertexSrc_;
      double etaMax_;
      double rapidityShift_;
      double vzMax_;
      bool printDebug_;

};

RpPbMCAnalyzer::RpPbMCAnalyzer(const edm::ParameterSet& iConfig):
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
etaMax_(iConfig.getParameter<double>("etaMax")),
rapidityShift_(iConfig.getParameter<double>("rapidityShift")),
vzMax_(iConfig.getParameter<double>("vzMax")),
printDebug_(iConfig.getParameter<bool>("printDebug"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;

   std::vector<Double_t> ptBins; std::vector<Double_t> ptBins_finer;
   const Double_t small = 1e-3;
   Double_t pt; Double_t pt_finer;

   for(pt =   0.2  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt);
   for(pt =   1.2; pt <   2.4-small; pt +=  0.1 ) ptBins.push_back(pt); //12 bins
   for(pt =   2.4; pt <   7.2-small; pt +=  0.4 ) ptBins.push_back(pt); //12 bins
   for(pt =   7.2; pt <  16.8-small; pt +=  1.6 ) ptBins.push_back(pt);//it was 3.2
   for(pt =  16.8; pt <  30.0-small; pt +=  6.6 ) ptBins.push_back(pt);
   ptBins.push_back(30.0);

   //finer bins
   for(pt_finer =   0.2  ; pt_finer <   1.2-small; pt_finer +=  0.005) ptBins_finer.push_back(pt_finer);
   for(pt_finer =   1.2; pt_finer <   2.4-small; pt_finer +=  0.01 ) ptBins_finer.push_back(pt_finer);
   for(pt_finer =   2.4; pt_finer <   7.2-small; pt_finer +=  0.04 ) ptBins_finer.push_back(pt_finer);
   for(pt_finer =   7.2; pt_finer <  16.8-small; pt_finer +=  0.16 ) ptBins_finer.push_back(pt_finer);
   for(pt_finer =  16.8; pt_finer <  30.0-small; pt_finer +=  0.66 ) ptBins_finer.push_back(pt_finer);
   ptBins_finer.push_back(30.0);

    events_ = fs->make<TH1F>("events","",10,0,10);

    vtxNtrk_ = fs->make<TH1F>("vtxNtrk","Tracks per vertex",200,0,200);
    vtxN_ = fs->make<TH1F>("vtxN","Vertices per event",20,0,20);

    distPV_ = fs->make<TH1F>("distPV_","3d distance from best primary vertex",100,0,20);
 
    {    
      std::vector<std::string> vlabel{"chargedHadrons","chargedHadrons+","chargedHadrons-",
                                      "pions","pions+","pions-",
                                      "kaons","kaons+","kaons-",
                                      "protons","protons+","protons-"};
      for( auto label : vlabel )
        spectra_[label] = fs->make<TH1F>(label.c_str(),
            Form("%s, %4.2f < #eta < %4.2f;p_{T} [GeV/c];(1/N_{evt})d^{2}N/dp_{T}d#eta",
                 label.c_str(), -etaMax_+rapidityShift_, etaMax_+rapidityShift_),
            ptBins.size()-1, &ptBins[0]);
    }
    {
      std::vector<std::string> vlabel{"chargedHadrons_finer","chargedHadrons_finer+","chargedHadrons_finer-",
                                      "pions_finer","pions_finer+","pions_finer-",
                                      "kaons_finer","kaons_finer+","kaons_finer-",
                                      "protons_finer","protons_finer+","protons_finer-"};
      for( auto label : vlabel )
        spectra_[label] = fs->make<TH1F>(label.c_str(),
            Form("%s, %4.2f < #eta < %4.2f;p_{T} [GeV/c];(1/N_{evt})d^{2}N/dp_{T}d#eta", 
                 label.c_str(), -etaMax_+rapidityShift_, etaMax_+rapidityShift_),
            ptBins_finer.size()-1, &ptBins_finer[0]);
    }

    nevt_ = 0;
}

RpPbMCAnalyzer::~RpPbMCAnalyzer()
{

  std::cout << "Total passed events: " << nevt_ << "\n\n";

}

void
RpPbMCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   Handle<reco::GenParticleCollection> gcol;
   iEvent.getByLabel(genSrc_, gcol);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   double vz =0.0, vx=0.0, vy=0.0;

   if(vertex->size()>0) {
     vz=vertex->begin()->z();
     vx=vertex->begin()->x();
     vy=vertex->begin()->y();
     } else { return; }

   unsigned int ntrk = 0;
   int vcount = 1;
   int vp = 1;

   for( const auto & vi : * vertex )
   {
     if( ntrk < vi.tracksSize() )
     { 
        ntrk = vi.tracksSize();
        vx = vi.x();
        vy = vi.y();
        vz = vi.z();
        vp = vcount;
     }
     vtxNtrk_->Fill(ntrk);
     if ( printDebug_ )
       std::cout << "reco::Vertex " << vcount << "  ntrk= " << vi.tracksSize()
                 << " x= " << vi.x()
                 << " y= " << vi.y()
                 << " z= " << vi.z() << std::endl;
     vcount++;
   }

   vtxN_->Fill(vertex->size() );

   if ( printDebug_ )
   {
     std::cout << "Total reconstructed vertices: " << vertex->size() << " using vertex " << vp << std::endl; 
     std::cout << "Primary vertex x= " << vx << " y= " << vy 
               << " z= " << vz << " ntrk= " << ntrk << std::endl;
   }

   // skip displaced vertices

   if( fabs(vz) > vzMax_ ) return;

   nevt_++;

   events_->Fill(0.5); 

   for( const auto & particle : * gcol )
   {
     if (printDebug_) printWeirdParticle( particle );

     if ( particle.status() == 1 && 
          particle.charge() != 0 &&
          fabs(particle.eta() - rapidityShift_ ) <= etaMax_ )
     {

       bool isPos = (particle.charge())==1 ? true : false;

       spectra_["chargedHadrons"]->Fill( particle.pt() );
       spectra_["chargedHadrons_finer"]->Fill( particle.pt() );
       if(isPos) spectra_["chargedHadrons+"]->Fill( particle.pt() );
       if(isPos) spectra_["chargedHadrons_finer+"]->Fill( particle.pt() );
       if(!isPos) spectra_["chargedHadrons-"]->Fill( particle.pt() );
       if(!isPos) spectra_["chargedHadrons_finer-"]->Fill( particle.pt() );
 
       double vdist = sqrt( (particle.vx() - vx)* (particle.vx() - vx)
                           + (particle.vy() - vy)* (particle.vy() - vy)
                           + (particle.vz() - vz)* (particle.vz() - vz) );
      
       distPV_->Fill(vdist);
 
       if ( abs(particle.pdgId()) == 211 )
       {
         spectra_["pions"]->Fill( particle.pt() );
         spectra_["pions_finer"]->Fill( particle.pt() );
         if(isPos) spectra_["pions+"]->Fill( particle.pt() );
         if(isPos) spectra_["pions_finer+"]->Fill( particle.pt() );
         if(!isPos) spectra_["pions-"]->Fill( particle.pt() );
         if(!isPos) spectra_["pions_finer-"]->Fill( particle.pt() );
         // find rogue pions
         if ( vdist > 1. && printDebug_) printRoguePion(particle,vdist);
       } 
       if ( abs(particle.pdgId()) == 321 )
       {
         spectra_["kaons"]->Fill( particle.pt() );
         spectra_["kaons_finer"]->Fill( particle.pt() );
         if(isPos) spectra_["kaons+"]->Fill( particle.pt() );
         if(isPos) spectra_["kaons_finer+"]->Fill( particle.pt() );
         if(!isPos) spectra_["kaons-"]->Fill( particle.pt() );
         if(!isPos) spectra_["kaons_finer-"]->Fill( particle.pt() );
       }
       if ( particle.pdgId() == 2212 )
       {
         spectra_["protons"]->Fill( particle.pt() );
         spectra_["protons_finer"]->Fill( particle.pt() );
         if(isPos) spectra_["protons+"]->Fill( particle.pt() );
         if(isPos) spectra_["protons_finer+"]->Fill( particle.pt() );
         if(!isPos) spectra_["protons-"]->Fill( particle.pt() );
         if(!isPos) spectra_["protons_finer-"]->Fill( particle.pt() );
       }
     }
   }
}

void 
RpPbMCAnalyzer::printWeirdParticle( const reco::GenParticle & particle)
{
   if (particle.status() > 3 ||
         ( particle.status() == 1 &&
           (  ( abs(particle.pdgId()) > 400 &&  abs(particle.pdgId()) < 600 )
             || ( abs(particle.pdgId()) > 4000 && abs(particle.pdgId()) < 6000 )
           )
         )
      )
   {
       std::cout << "  Stable B/D status= " << particle.status()
                 << " pdgId= " << particle.pdgId()
                 << " pt= " << particle.pt()
                 << " eta= " << particle.eta()
                 << " phi= " << particle.phi()
                 << std::endl;
   }
}

void 
RpPbMCAnalyzer::printRoguePion( const reco::GenParticle & particle, double vdist)
{
  std::cout << "  Rogue pion px= " << particle.vx()
            << " py= " << particle.vy()
            << " pz= " << particle.vz()
            << " vdist= " << vdist
            << " pt= " << particle.pt()
            << " eta= " << particle.eta()
            << " phi= " << particle.phi()
            << std::endl;
  if( particle.numberOfMothers() > 0 )
    std::cout << "  Mother pdgId= " << particle.mother()->pdgId() << std::endl;
}

void
RpPbMCAnalyzer::beginJob()
{
}

void
RpPbMCAnalyzer::endJob()
{
  // properly scale spectra histograms
  for( auto elem : spectra_ )
  {
    auto & histo = elem.second;
    for(int i = 1; i <= histo->GetNbinsX(); i++) {
      float content = histo->GetBinContent(i);
      float error = histo->GetBinError(i);
      histo->SetBinContent(i,content/histo->GetBinWidth(i)/(float)nevt_/2./etaMax_);
      histo->SetBinError(i,error/histo->GetBinWidth(i)/(float)nevt_/2./etaMax_);
    }
  }  
}

DEFINE_FWK_MODULE(RpPbMCAnalyzer);
