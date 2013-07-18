//
// Original Author:  Eric Appelt
//         Created:  Fri Nov 19 16:55:16 CST 2010
//
//

#include <memory>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/Scalers/interface/LumiScalers.h>

class RpPbVertexAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RpPbVertexAnalyzer(const edm::ParameterSet&);
      ~RpPbVertexAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);

      // ----------member data ---------------------------

      std::map<std::string,TH1F*> evtPerf_;
      std::map<std::string,TH2F*> evtPerf2D_;
      std::map<std::string,TH1F*> vtxPerf_;
      std::map<std::string,TH2F*> vtxPerf2D_;
      std::map<std::string,TH3F*> vtxPerf3D_;

      TH1F* events_;
      TH1F* vertices_;

      int nevt_;
      int nvertex_;
      bool newLumi_;

      edm::InputTag vertexSrc_;
      CentralityProvider * centrality_;

 
};

RpPbVertexAnalyzer::RpPbVertexAnalyzer(const edm::ParameterSet& iConfig):
nevt_(0),
nvertex_(0),
newLumi_(true),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
   centrality_ = 0;
}

RpPbVertexAnalyzer::~RpPbVertexAnalyzer()
{
}

void
RpPbVertexAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const & b, edm::EventSetup const & e )
{
  newLumi_  = true;
}

void
RpPbVertexAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);
  
   std::vector<reco::Vertex> vsorted = *vertex;
   // sort the vertcies by number of tracks in descending order
   std::sort( vsorted.begin(), vsorted.end(), 
              []( reco::Vertex a, reco::Vertex b) 
   {
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   });

   nevt_++;
   events_->Fill(0.5); 
   evtPerf_["Nvtx"]->Fill(vsorted.size());

   int lumi = iEvent.getLuminosityBlock().luminosityBlock();
   evtPerf_["Lumi"]->Fill(lumi);
   evtPerf_["NvtxLumi"]->Fill(lumi,vsorted.size());

   if( newLumi_ ) 
   {
     Handle<std::vector<LumiScalers> > lcol;
     iEvent.getByLabel("scalersRawToDigi", lcol);
 
     if( lcol->size() > 0 )
       evtPerf_["puLumi"]->Fill(lumi,lcol->begin()->pileup());
     newLumi_ = false;
   } 
 
   if (!centrality_) centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup); 
   evtPerf_["centrality"]->Fill(centrality_->getBin());

   int vcount = 0; 
   for( const auto & vi :  vsorted )
   {
     vtxPerf_["Ntrk"]->Fill(vi.tracksSize());
     vtxPerf_["x"]->Fill(vi.x());
     vtxPerf_["y"]->Fill(vi.y());
     vtxPerf_["z"]->Fill(vi.z());
     vtxPerf_["xerr"]->Fill(vi.xError());
     vtxPerf_["yerr"]->Fill(vi.yError());
     vtxPerf_["zerr"]->Fill(vi.zError());
     vtxPerf_["chi2"]->Fill(vi.normalizedChi2());
     vtxPerf2D_["Ntrk2D"]->Fill(vcount,vi.tracksSize());
     vertices_->Fill(0.5);
     vcount++;
     nvertex_++;
   }

   for (unsigned int i =1; i<vsorted.size(); i++)
   {
     double dz = fabs( vsorted[i].z() - vsorted[0].z() );
     double dx = fabs( vsorted[i].x() - vsorted[0].x() );
     double dy = fabs( vsorted[i].y() - vsorted[0].y() );
   //  double dxErr = sqrt( vsorted[i].xError()*vsorted[i].xError() + vsorted[0].xError()*vsorted[0].xError() );
   //  double dyErr = sqrt( vsorted[i].yError()*vsorted[i].yError() + vsorted[0].yError()*vsorted[0].yError() );
   //  double dzErr = sqrt( vsorted[i].zError()*vsorted[i].zError() + vsorted[0].zError()*vsorted[0].zError() );
     double dxy  = sqrt ( dx*dx + dy*dy );
   //  double dxyErr = sqrt ( dxErr*dxErr + dyErr * dyErr );
     vtxPerf_["assocVtxDz"]->Fill(dz);
     vtxPerf_["assocVtxDxy"]->Fill(dxy);
     vtxPerf2D_["assocVtxDzNtrk"]->Fill(dz,vsorted[i].tracksSize() );
     vtxPerf2D_["assocVtxDxyNtrk"]->Fill(dxy,vsorted[i].tracksSize() );
     vtxPerf2D_["assocVtxDxyDz"]->Fill(dxy,dz);
     vtxPerf2D_["assocVtxChi2Ntrk"]->Fill(vsorted[i].normalizedChi2(),vsorted[i].tracksSize());
     vtxPerf2D_["assocVtxChi2Dxy"]->Fill(vsorted[i].normalizedChi2(),dxy); 
     vtxPerf2D_["assocVtxChi2Dz"]->Fill(vsorted[i].normalizedChi2(),dz); 
     vtxPerf2D_["vtxCorrZ"]->Fill( vsorted[0].z(), vsorted[i].z() );
     vtxPerf3D_["assocVtxDzNtrkNtrkLead"]->Fill(dz,vsorted[i].tracksSize(),vsorted[0].tracksSize()); 
   }

}


void
RpPbVertexAnalyzer::initHistos(const edm::Service<TFileService> & fs)
{

  events_ = fs->make<TH1F>("events","",1,0,1);
  vertices_ = fs->make<TH1F>("vertices","",1,0,1);

  evtPerf_["Ntrk"] = fs->make<TH1F>("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["Nvtx"] = fs->make<TH1F>("evtNvtx","Primary Vertices per event",10,0,10);
  evtPerf_["centrality"] = fs->make<TH1F>("centrality","Event centrality bin",100,0,100);

  evtPerf_["NvtxLumi"] = fs->make<TH1F>("evtNvtxLumi","Primary Vertices by Lumi",200,0,2000);
  evtPerf_["Lumi"] = fs->make<TH1F>("evtLumi","Events by Lumi",200,0,2000);
  evtPerf_["puLumi"] = fs->make<TH1F>("puLumi","Pileup by Lumi",2000,0,2000);

  vtxPerf_["Ntrk"] = fs->make<TH1F>("vtxNtrk","Tracks per vertex",50,0,200);
  vtxPerf_["x"] = fs->make<TH1F>("vtxX","Vertex x position",1000,-1,1);
  vtxPerf_["y"] = fs->make<TH1F>("vtxY","Vertex y position",1000,-1,1);
  vtxPerf_["z"] = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);
  vtxPerf_["xerr"] = fs->make<TH1F>("vtxXerr","Vertex x error",100,0.,0.05);
  vtxPerf_["yerr"] = fs->make<TH1F>("vtxYerr","Vertex y error",100,0.,0.05);
  vtxPerf_["zerr"] = fs->make<TH1F>("vtxZerr","Vertex z error",100,0.,0.05);
  vtxPerf_["chi2"] = fs->make<TH1F>("chi2","Vertex Normalized Chi2",200,0.,50.);
  vtxPerf_["assocVtxDz"] = fs->make<TH1F>("assocVtxDz","Z Distance from first PV; dz (cm)",500,0,50);
  vtxPerf_["assocVtxDxy"] = fs->make<TH1F>("assocVtxDxy","Rho Distance from first PV; dxy (cm)",800,0,4);

  vtxPerf2D_["Ntrk2D"] = fs->make<TH2F>("vtxNtrk2D","Tracks per vertex;vertex (sorted by Ntrk);Ntrk"
                                ,10,0,10,200,0,200);
  vtxPerf2D_["assocVtxDzNtrk"] = fs->make<TH2F>("assocVtxDzNtrk",
                                 "Z Distance from first PV vs Ntrk of assoc; dz (cm); Ntrk",
                                 500,0,50,200,0,200);
  vtxPerf2D_["assocVtxDxyNtrk"] = fs->make<TH2F>("assocVtxDxyNtrk",
                                 "Rho Distance from first PV vs Ntrk of assoc; dxy (cm); Ntrk",
                                 800,0,4,200,0,200);
  vtxPerf2D_["assocVtxChi2Ntrk"] = fs->make<TH2F>("assocVtxChi2Ntrk",
                                 "Normalized Chi2 vs Ntrk of assoc; chi2/ndof; Ntrk",
                                 200,0,50,200,0,200);
  vtxPerf2D_["assocVtxDxyDz"] = fs->make<TH2F>("assocVtxDxyDz",
                                 "Rho Distance from first PV vs Z distance; dxy (cm); dz (cm)",
                                 800,0,4,500,0,50);
  vtxPerf2D_["assocVtxChi2Dz"] = fs->make<TH2F>("assocVtxChi2Dz",
                                 "Normalized Chi2 vs Z distance; chi2/ndof; dz (cm)",
                                 200,0,50.,500,0,50);
  vtxPerf2D_["assocVtxChi2Dxy"] = fs->make<TH2F>("assocVtxChi2Dxy",
                                 "Normalized Chi2 vs Rho distance; chi2/ndof; dxy (cm)",
                                 200,0,50.,800,0,4);
  vtxPerf2D_["vtxCorrZ"] = fs->make<TH2F>("vtzCorrZ",
                                 "z position of first PV vs additional PVs; z_{trig} (cm); z_{assoc} (cm)",
                                 300,-30,30,300,-30,30);

  vtxPerf3D_["assocVtxDzNtrkNtrkLead"] = fs->make<TH3F>("assocVtxDzNtrkNtrkLead",
                                 "Z Distance from first PV vs Ntrk of assoc vs Ntrk of lead; dz (cm); Ntrk; NtrkLead",
                                 500,0,50,200,0,200,20,0,400);

}




void
RpPbVertexAnalyzer::beginJob()
{
}

void
RpPbVertexAnalyzer::endJob()
{
}

DEFINE_FWK_MODULE(RpPbVertexAnalyzer);
