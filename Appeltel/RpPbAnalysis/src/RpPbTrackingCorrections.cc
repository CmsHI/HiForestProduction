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
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"

class RpPbTrackingCorrections : public edm::EDAnalyzer {
   public:
      explicit RpPbTrackingCorrections(const edm::ParameterSet&);
      ~RpPbTrackingCorrections();

      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);
      bool passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex);
      double bestJetEt( const reco::Track & tr, 
                        const std::vector<const pat::Jet *> & jets );
      double bestJetEt( const TrackingParticle & tp, 
                        const std::vector<const pat::Jet *> & jets );
      void fillGenHistosWithTp( const TrackingParticle & tp, double occ, double w );
      void fillTrkPerfHistosWithTrk( const reco::Track & track, const reco::Vertex & vertex, double w);
      int pixelLayersHit( const TrackingParticle & tp );
      int getLayerId(const PSimHit & simHit);

      // ----------member data ---------------------------

      edm::ESHandle<TrackerGeometry> trackerGeo;

      TrackClassifier classifier_;

      std::map<std::string,TH3F*> trkCorr3D_;
      std::map<std::string,TH2F*> trkCorr2D_;
      std::map<std::string,TH3F*> trkPerf3D_;
      std::map<std::string,TH3F*> genHist3D_;
      TH3F * momRes_;
      TH1F * vtxZ_;
      TH1F * leadJetEt_;
      TF1 * vtxWeightFunc_;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      edm::InputTag tpFakSrc_;
      edm::InputTag tpEffSrc_;
      edm::InputTag associatorMap_;
      edm::InputTag jetSrc_;

      std::vector<double> ptBins_;
      std::vector<double> etaBins_;
      std::vector<double> occBins_;

      
      std::vector<double> vtxWeightParameters_;
      bool doVtxReweighting_;

      bool occByCentrality_;
      bool occByNPixelTrk_;
      bool occByJetEt_;
      bool occByLeadingJetEt_;

      double jetEtaMax_;
      double jetRadius_;

      CentralityProvider * centrality_;

      bool selectSpecies_;
      std::vector<int> pdgIdList_;
      bool selectHeavyFlavorDecays_;
      bool selectNonHeavyFlavorDecays_;

      bool applyVertexZCut_;
      double vertexZMax_;

      bool applyJetCuts_;
      bool invertJetCuts_;
      double jetEtMin_;
      double jetEtMax_;

      bool applyTrackCuts_;
      std::string qualityString_;
      double dxyErrMax_;
      double dzErrMax_;
      double ptErrMax_;

      bool fillGenHistos_;
      bool fillTrkPerfHistos_;
      bool doMomRes_;

      enum { BPix1=0, BPix2=1, BPix3=2,
	     FPix1_neg=3, FPix2_neg=4,
	     FPix1_pos=5, FPix2_pos=6, nLayers=7};

};

RpPbTrackingCorrections::RpPbTrackingCorrections(const edm::ParameterSet& iConfig):
classifier_(iConfig),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
tpFakSrc_(iConfig.getParameter<edm::InputTag>("tpFakSrc")),
tpEffSrc_(iConfig.getParameter<edm::InputTag>("tpEffSrc")),
associatorMap_(iConfig.getParameter<edm::InputTag>("associatorMap")),
jetSrc_(iConfig.getParameter<edm::InputTag>("jetSrc")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
occBins_(iConfig.getParameter<std::vector<double> >("occBins")),
vtxWeightParameters_(iConfig.getParameter<std::vector<double> >("vtxWeightParameters")),
doVtxReweighting_(iConfig.getParameter<bool>("doVtxReweighting")),
occByCentrality_(iConfig.getParameter<bool>("occByCentrality")),
occByNPixelTrk_(iConfig.getParameter<bool>("occByNPixelTrk")),
occByJetEt_(iConfig.getParameter<bool>("occByJetEt")),
occByLeadingJetEt_(iConfig.getParameter<bool>("occByLeadingJetEt")),
jetEtaMax_(iConfig.getParameter<double>("jetEtaMax")),
jetRadius_(iConfig.getParameter<double>("jetRadius")),
selectSpecies_(iConfig.getParameter<bool>("selectSpecies")),
pdgIdList_(iConfig.getParameter<std::vector<int> >("pdgIdList")),
selectHeavyFlavorDecays_(iConfig.getParameter<bool>("selectHeavyFlavorDecays")),
selectNonHeavyFlavorDecays_(iConfig.getParameter<bool>("selectNonHeavyFlavorDecays")),
applyVertexZCut_(iConfig.getParameter<bool>("applyVertexZCut")),
vertexZMax_(iConfig.getParameter<double>("vertexZMax")),
applyJetCuts_(iConfig.getParameter<bool>("applyJetCuts")),
invertJetCuts_(iConfig.getParameter<bool>("invertJetCuts")),
jetEtMin_(iConfig.getParameter<double>("jetEtMin")),
jetEtMax_(iConfig.getParameter<double>("jetEtMax")),
applyTrackCuts_(iConfig.getParameter<bool>("applyTrackCuts")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
dxyErrMax_(iConfig.getParameter<double>("dzErrMax")),
dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
ptErrMax_(iConfig.getParameter<double>("ptErrMax")),
fillGenHistos_(iConfig.getParameter<bool>("fillGenHistos")),
fillTrkPerfHistos_(iConfig.getParameter<bool>("fillTrkPerfHistos")),
doMomRes_(iConfig.getParameter<bool>("doMomRes"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
   centrality_ = 0;
   vtxWeightFunc_ = new TF1("vtxWeight","gaus(0)/gaus(3)",-50.,50.);
   // vtxWeightParameters should have size 6,
   // one really should throw an error if not
   if( (int)vtxWeightParameters_.size() == 6 )
   {
     for( unsigned int i=0;i<vtxWeightParameters_.size(); i++)
       vtxWeightFunc_->SetParameter(i,vtxWeightParameters_[i]);
   }
}

RpPbTrackingCorrections::~RpPbTrackingCorrections()
{
   delete vtxWeightFunc_;
}

void
RpPbTrackingCorrections::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Setup Classifier if it is being used
   if( selectHeavyFlavorDecays_ || selectNonHeavyFlavorDecays_)
     classifier_.newEvent(iEvent, iSetup);

   // obtain collections of simulated particles 
   edm::Handle<TrackingParticleCollection>  TPCollectionHeff, TPCollectionHfake;
   iEvent.getByLabel(tpEffSrc_,TPCollectionHeff);
   iEvent.getByLabel(tpFakSrc_,TPCollectionHfake);

   // obtain association map between tracks and simulated particles
   reco::RecoToSimCollection recSimColl;
   reco::SimToRecoCollection simRecColl;
   edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
   edm::Handle<reco::RecoToSimCollection > recotosimCollectionH;
   iEvent.getByLabel(associatorMap_,simtorecoCollectionH);
   simRecColl= *(simtorecoCollectionH.product());
   iEvent.getByLabel(associatorMap_,recotosimCollectionH);
   recSimColl= *(recotosimCollectionH.product());

   // obtain reconstructed tracks
   Handle<edm::View<reco::Track> > tcol;
   iEvent.getByLabel(trackSrc_, tcol);

   //obtain tracker geometry
   iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeo);

   // obtain jets, if we are considering them
   Handle<std::vector<pat::Jet> > jets;
   std::vector<const pat::Jet *> sortedJets;
   if( occByJetEt_  || applyJetCuts_ || occByLeadingJetEt_) 
   {
      iEvent.getByLabel(jetSrc_, jets);
      // take only accepted jets and sort them by Et
      for(unsigned it=0; it<jets->size(); ++it)
      {
        const pat::Jet* jet = &((*jets)[it]);
        if(fabs(jet->eta())>=jetEtaMax_) continue; 
        sortedJets.push_back(jet);
      }
      sortByPtRef(&sortedJets); 
   }
  
   // obtain primary vertices
   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);
  
   // sort the vertcies by number of tracks in descending order
   std::vector<reco::Vertex> vsorted = *vertex;
   std::sort( vsorted.begin(), vsorted.end(), RpPbTrackingCorrections::vtxSort );

   // obtain event centrality
   if (!centrality_) centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup); 

   // skip events with no PV, this should not happen
   if( vsorted.size() == 0) return;

   // skip events failing vertex cut
   if( applyVertexZCut_)
   {
     if( fabs(vsorted[0].z()) > vertexZMax_ ) return;
   }

   // determine vertex reweighting factor
   double w = 1.0;
   if ( doVtxReweighting_ )
     w *= vtxWeightFunc_->Eval(vsorted[0].z());

   vtxZ_->Fill(vsorted[0].z(),w);

   // determine occupancy variable for event
   double occ = 0.;  
   if( occByCentrality_) occ = (double) centrality_->getBin(); 
   if( occByNPixelTrk_) occ = centrality_->raw()->NpixelTracks();   
   double leadJetEt = 0.;
   if ( ! sortedJets.empty() ) leadJetEt = sortedJets[0]->pt();
   if( occByLeadingJetEt_ ) occ = leadJetEt;
   leadJetEt_->Fill( leadJetEt ,w );
   // ---------------------
   // loop through reco tracks to fill fake, reco, and secondary histograms
   // ---------------------

   for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){
    
     edm::RefToBase<reco::Track> track(tcol, i);
     reco::Track* tr=const_cast<reco::Track*>(track.get());
     // skip tracks that fail cuts, using vertex with most tracks as PV       
     if( !passesTrackCuts(*tr, vsorted[0]) ) continue;

     // if applying jet cuts, only take tracks within jet 
     // cone of a jet within Et range
     if ( applyJetCuts_ )
     {  
        bool cut = false;
        double jetEt = bestJetEt( *tr, sortedJets);
        if( jetEt < jetEtMin_ || jetEt > jetEtMax_ ) cut = true;
        if( invertJetCuts_ ) cut = !cut;
        if( cut ) continue;
     } 

     if( fillTrkPerfHistos_ ) fillTrkPerfHistosWithTrk( *tr, vsorted[0], w) ;

     // get Et of closest jet to track, or 0 if not in a cone, 
     // if jet occupancy is used
     if ( occByJetEt_ ) occ = bestJetEt( *tr, sortedJets);

     trkCorr3D_["hrec3D"]->Fill(tr->eta(), tr->pt(), occ,w);
     trkCorr2D_["hrec"]->Fill(tr->eta(), tr->pt(), w);

     // look for match to simulated particle, use first match if it exists
     std::vector<std::pair<TrackingParticleRef, double> > tp;
     const TrackingParticle *mtp=0;
     if(recSimColl.find(track) != recSimColl.end())
     {
       tp = recSimColl[track];
       mtp = tp.begin()->first.get();  
       if( mtp->status() < 0) trkCorr3D_["hsec3D"]->Fill(tr->eta(), tr->pt(), occ, w);     
       if( mtp->status() < 0) trkCorr2D_["hsec"]->Fill(tr->eta(), tr->pt(), w);     
     }
     else
     {
       trkCorr3D_["hfak3D"]->Fill(tr->eta(), tr->pt(), occ, w);
       trkCorr2D_["hfak"]->Fill(tr->eta(), tr->pt(), w);
     }

   }

   // ---------------------
   // loop through sim particles to fill matched, multiple,  and sim histograms 
   // ---------------------

   for(TrackingParticleCollection::size_type i=0; i<TPCollectionHeff->size(); i++) 
   {      
     TrackingParticleRef tpr(TPCollectionHeff, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
         
     if(tp->status() < 0 || tp->charge()==0) continue; //only charged primaries
 
     // select certain species for efficiency
     if( selectSpecies_ ) 
     { 
        bool accept = false;
        for( const int & id : pdgIdList_ )
        {
          if( tp->pdgId() == id ) accept = true; 
        }
        if(accept == false) continue;  
     }

     // select D and B if set
     if( selectHeavyFlavorDecays_ || selectNonHeavyFlavorDecays_)
     {
        bool fromHF = false;
        classifier_.evaluate(tpr);
        if( classifier_.is(TrackCategories::BWeakDecay) || 
            classifier_.is(TrackCategories::CWeakDecay)    ) fromHF = true;
        if( fromHF && selectNonHeavyFlavorDecays_ ) continue;
        if( !fromHF && selectHeavyFlavorDecays_ ) continue;
     }
    

     // if applying jet cuts, only take particles within jet 
     // cone of a jet within Et range
     if ( applyJetCuts_ )
     {
        bool cut = false;
        double jetEt = bestJetEt( *tp, sortedJets);
        if( jetEt < jetEtMin_ || jetEt > jetEtMax_ ) cut = true;
        if( invertJetCuts_ ) cut = !cut;
        if( cut ) continue;
     }

     // get Et of closest jet to track, or 0 if not in a cone, 
     // if jet occupancy is used 
     if ( occByJetEt_ ) occ = bestJetEt( *tp, sortedJets);

     trkCorr3D_["hsim3D"]->Fill(tp->eta(),tp->pt(), occ, w);
     trkCorr2D_["hsim"]->Fill(tp->eta(),tp->pt(), w);
     if( fillGenHistos_ ) fillGenHistosWithTp( *tp, occ, w );

     if( pixelLayersHit(*tp) >= 2) 
     {
       trkCorr3D_["hacc3D"]->Fill(tp->eta(),tp->pt(), occ, w);
       trkCorr2D_["hacc"]->Fill(tp->eta(),tp->pt(), w);
     }

     // find number of matched reco tracks that pass cuts
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
     size_t nrec=0;
     if(simRecColl.find(tpr) != simRecColl.end())
     {
       rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[tpr];
       for (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator rtit = rt.begin(); rtit != rt.end(); ++rtit)
       {
         const reco::Track* tmtr = rtit->first.get();
         if( ! passesTrackCuts(*tmtr, vsorted[0]) ) continue;
         nrec++;
         if( doMomRes_ ) momRes_->Fill( tp->eta(), tp->pt(), tmtr->pt(), w);
       }
     }
     if(nrec>0) trkCorr3D_["heff3D"]->Fill(tp->eta(),tp->pt(), occ, w);
     if(nrec>1) trkCorr3D_["hmul3D"]->Fill(tp->eta(),tp->pt(), occ, w);
     if(nrec>0) trkCorr2D_["heff"]->Fill(tp->eta(),tp->pt(), w);
     if(nrec>1) trkCorr2D_["hmul"]->Fill(tp->eta(),tp->pt(), w);
   }
}

bool
RpPbTrackingCorrections::passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex)
{
   if ( ! applyTrackCuts_ ) return true;

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   vtxPoint=vertex.position();
   vzErr=vertex.zError();
   vxErr=vertex.xError();
   vyErr=vertex.yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
 
   if(track.quality(reco::TrackBase::qualityByName(qualityString_)) != 1)
       return false;
   if(fabs(dxy/dxysigma) > dxyErrMax_) return false;
   if(fabs(dz/dzsigma) > dzErrMax_) return false;
   if(track.ptError() / track.pt() > ptErrMax_) return false;

   return true;
}


void
RpPbTrackingCorrections::initHistos(const edm::Service<TFileService> & fs)
{

  std::vector<std::string> hNames3D = { "hsim3D", "hrec3D", "hmul3D", "hfak3D",
                                        "heff3D", "hsec3D", "hacc3D" };

  for( auto name : hNames3D )
  {
     trkCorr3D_[name] = fs->make<TH3F>(name.c_str(),";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]); 
  }

  std::vector<std::string> hNames2D = { "hsim", "hrec", "hmul", "hfak",
                                        "heff", "hsec", "hacc" };

  for( auto name : hNames2D )
  {
     trkCorr2D_[name] = fs->make<TH2F>(name.c_str(),";#eta;p_{T}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0]);
  }

  std::vector<std::string> genNames3D = { "genPartH","genPartH+","genPartH-",
                                          "genPartPi","genPartPi+","genPartPi-",
                                          "genPartK","genPartK+","genPartK-",
                                          "genPartp","genPartp+","genPartp-",
                                          "genPartSigma","genPartXi","genPartOmega"
                                        };

  for( auto name : genNames3D )
  {
     genHist3D_[name] = fs->make<TH3F>(name.c_str(),";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]);
  }

  std::vector<double> dumBins;
  dumBins.clear(); 
  for( double i = 0.; i<36.; i += 1.) dumBins.push_back(i);
  trkPerf3D_["Nhit3D"] = fs->make<TH3F>("trkNhit3D", "Tracks by Number of Valid Hits;N hits",    
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();  
  for( double i = -3.15; i<3.151; i += 6.30/50.) dumBins.push_back(i);
  trkPerf3D_["phi3D"] = fs->make<TH3F>("trkPhi3D", "Track Azimuthal Distribution;#phi",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();    
  for( double i = 0.; i<6.01; i += 6./60.) dumBins.push_back(i);
  trkPerf3D_["chi23D"] = fs->make<TH3F>("trkChi23D", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = 0.0; i<0.201; i += 0.2/50.) dumBins.push_back(i);
  trkPerf3D_["pterr3D"] = fs->make<TH3F>("trkPterr3D", "Track p_{T} error;#delta p_{T} / p_{T}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = -8.; i<8.01; i += 16./50.) dumBins.push_back(i);
  trkPerf3D_["dxyErr3D"] = fs->make<TH3F>("trkDxyErr3D", "Transverse DCA Significance;dxy / #sigma_{dxy}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  trkPerf3D_["dzErr3D"] = fs->make<TH3F>("trkDzErr3D", "Longitudinal DCA Significance;dz / #sigma_{dz}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);

  vtxZ_ = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);
  leadJetEt_ = fs->make<TH1F>("leadJetEt","Vertex z position",1000,0,1000);

  std::vector<double> ptBinsFine;
  for( unsigned int bin = 0; bin<ptBins_.size()-1; bin++)
  {
    double bStart = ptBins_[bin];
    double bWid = ptBins_[bin+1] - ptBins_[bin];
    for( int i=0;i<5;i++)
      ptBinsFine.push_back( bStart + (double)i * bWid / 5. );
  }
  ptBinsFine.push_back(ptBins_[ptBins_.size()-1]);

  momRes_ = fs->make<TH3F>("momRes","momentum resolution sim vs reco",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBinsFine.size()-1, &ptBinsFine[0],
                           ptBinsFine.size()-1, &ptBinsFine[0]);

}

bool
RpPbTrackingCorrections::vtxSort( const reco::Vertex &  a, const reco::Vertex & b )
{
  if( a.tracksSize() != b.tracksSize() )
    return  a.tracksSize() > b.tracksSize() ? true : false ;
  else
    return  a.chi2() < b.chi2() ? true : false ;  
}

double
RpPbTrackingCorrections::bestJetEt( const reco::Track & tr,
                         const std::vector<const pat::Jet *> & jets )
{
   double drMin = 999.;
   double et = 0.;
   for( const auto * jet : jets )
   {
     double dr = deltaR(*jet,tr);
     if( dr < drMin && dr < jetRadius_ )
     {
       drMin = dr; et = jet->pt();
     }
   }
   return et;
}

double
RpPbTrackingCorrections::bestJetEt( const TrackingParticle & tp,
                         const std::vector<const pat::Jet *> & jets )
{
   double drMin = 999.;
   double et = 0.;
   for( const auto * jet : jets )
   {
     double dr = deltaR(*jet,tp);
     if( dr < drMin && dr < jetRadius_ )
     {
       drMin = dr; et = jet->pt();
     }
   }
   return et;
}

void
RpPbTrackingCorrections::beginJob()
{
}

void
RpPbTrackingCorrections::endJob()
{
}

void 
RpPbTrackingCorrections::fillGenHistosWithTp( const TrackingParticle & tp, double occ, double w )
{
  genHist3D_["genPartH"]->Fill( tp.eta(), tp.pt(), occ, w);
  if( tp.charge() == 1) genHist3D_["genPartH+"]->Fill( tp.eta(), tp.pt(), occ, w);
  if( tp.charge() == -1) genHist3D_["genPartH-"]->Fill( tp.eta(), tp.pt(), occ, w);

  if( fabs(tp.pdgId()) == 211 )
  {
    genHist3D_["genPartPi"]->Fill( tp.eta(), tp.pt(), occ, w);
    if( tp.charge() == 1) genHist3D_["genPartPi+"]->Fill( tp.eta(), tp.pt(), occ, w);
    if( tp.charge() == -1) genHist3D_["genPartPi-"]->Fill( tp.eta(), tp.pt(), occ, w);
  }
  if( fabs(tp.pdgId()) == 321 )
  {
    genHist3D_["genPartK"]->Fill( tp.eta(), tp.pt(), occ);
    if( tp.charge() == 1) genHist3D_["genPartK+"]->Fill( tp.eta(), tp.pt(), occ, w);
    if( tp.charge() == -1) genHist3D_["genPartK-"]->Fill( tp.eta(), tp.pt(), occ, w);
  }
  if( fabs(tp.pdgId()) == 2212 )
  {
    genHist3D_["genPartp"]->Fill( tp.eta(), tp.pt(), occ);
    if( tp.charge() == 1) genHist3D_["genPartp+"]->Fill( tp.eta(), tp.pt(), occ, w);
    if( tp.charge() == -1) genHist3D_["genPartp-"]->Fill( tp.eta(), tp.pt(), occ, w);
  }
  if( fabs(tp.pdgId()) == 3222 || fabs(tp.pdgId()) == 3112 )
    genHist3D_["genPartSigma"]->Fill( tp.eta(), tp.pt(), occ, w);
  if( fabs(tp.pdgId()) == 3312  )
    genHist3D_["genPartXi"]->Fill( tp.eta(), tp.pt(), occ, w);
  if( fabs(tp.pdgId()) == 3334  )
    genHist3D_["genPartOmega"]->Fill( tp.eta(), tp.pt(), occ, w);
}


int
RpPbTrackingCorrections::pixelLayersHit(const TrackingParticle & tp)
{
  std::vector<bool> f(nLayers, false);

  const std::vector<PSimHit> & simHits = tp.trackPSimHit(DetId::Tracker);
  
  for(std::vector<PSimHit>::const_iterator simHit = simHits.begin();
                                      simHit!= simHits.end(); simHit++)
  {
    int id = getLayerId(*simHit);

    if(id != -1)
      f[id] = true;
  }

  bool canBeTriplet =
    ( (f[BPix1] && f[BPix2]     && f[BPix3]) ||
      (f[BPix1] && f[BPix2]     && f[FPix1_pos]) ||
      (f[BPix1] && f[BPix2]     && f[FPix1_neg]) ||
      (f[BPix1] && f[FPix1_pos] && f[FPix2_pos]) ||
      (f[BPix1] && f[FPix1_neg] && f[FPix2_neg]) );
  if(canBeTriplet) return 3;

  bool canBePair =
    ( (f[BPix1] && f[BPix2]) ||
      (f[BPix1] && f[BPix3]) ||
      (f[BPix2] && f[BPix3]) ||
      (f[BPix1] && f[FPix1_pos]) ||
      (f[BPix1] && f[FPix1_neg]) ||
      (f[BPix1] && f[FPix2_pos]) ||
      (f[BPix1] && f[FPix2_neg]) ||
      (f[BPix2] && f[FPix1_pos]) ||
      (f[BPix2] && f[FPix1_neg]) ||
      (f[BPix2] && f[FPix2_pos]) ||
      (f[BPix2] && f[FPix2_neg]) ||
      (f[FPix2_neg] && f[FPix2_neg]) ||
      (f[FPix2_pos] && f[FPix2_pos]) );
  if(canBePair) return 2;

  bool canBeSingle =
    ( f[BPix1] || f[BPix2] || f[BPix3] ||
      f[FPix1_pos] || f[FPix1_neg] ||
      f[FPix2_pos] || f[FPix2_neg] );
  if(canBeSingle) return 1;


  return 0;
}

// ------------
int 
RpPbTrackingCorrections::getLayerId(const PSimHit & simHit)
{
  unsigned int id = simHit.detUnitId();

  if(trackerGeo->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelBarrel)
  {
    PXBDetId pid(id);
    return pid.layer() - 1; // 0, 1, 2
  }

  if(trackerGeo->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelEndcap)
  {
    PXFDetId pid(id);
    return BPix3 + ((pid.side()-1) << 1) + pid.disk(); // 3 -
  }

  // strip
  return -1;
}

void
RpPbTrackingCorrections::fillTrkPerfHistosWithTrk( const reco::Track & track, const reco::Vertex & vertex, double w )
{
   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   vtxPoint=vertex.position();
   vzErr=vertex.zError();
   vxErr=vertex.xError();
   vyErr=vertex.yError();
   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
   trkPerf3D_["Nhit3D"]->Fill(track.eta(), track.pt(), track.numberOfValidHits(), w);
   trkPerf3D_["phi3D"]->Fill(track.eta(), track.pt(), track.phi(),w );
   trkPerf3D_["dxyErr3D"]->Fill(track.eta(), track.pt(), dxy/dxysigma, w);
   trkPerf3D_["dzErr3D"]->Fill(track.eta(), track.pt(), dz/dzsigma, w);
   trkPerf3D_["chi23D"]->Fill(track.eta(), track.pt(), track.normalizedChi2(), w);
   trkPerf3D_["pterr3D"]->Fill(track.eta(), track.pt(), track.ptError() / track.pt(), w );
}

DEFINE_FWK_MODULE(RpPbTrackingCorrections);
