//
// Original Author:  Edward Wenger
//         Created:  Thu Apr 29 14:31:47 CEST 2010
// $Id: HiTrkEffAnalyzer.cc,v 1.18 2011/11/30 18:06:17 sungho Exp $
//

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "edwenger/HiTrkEffAnalyzer/interface/HiTrkEffAnalyzer.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <stdio.h>
#include <math.h>

//#define DEBUG

HiTrkEffAnalyzer::HiTrkEffAnalyzer(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
  jetTags_(iConfig.getUntrackedParameter<edm::InputTag>("jets")),
  label_tp_effic_(iConfig.getUntrackedParameter<edm::InputTag>("label_tp_effic")),
  label_tp_fake_(iConfig.getUntrackedParameter<edm::InputTag>("label_tp_fake")),
  associatorMap_(iConfig.getUntrackedParameter<edm::InputTag>("associatormap")),
  vtxTags_(iConfig.getUntrackedParameter<edm::InputTag>("vertices")),
  bsTags_(iConfig.getUntrackedParameter<edm::InputTag>("beamspot")),
  doAssociation_(iConfig.getUntrackedParameter<bool>("doAssociation",true)),
  hasSimInfo_(iConfig.getUntrackedParameter<bool>("hasSimInfo",false)),
  pixelMultMode_(iConfig.getUntrackedParameter<bool>("pixelMultMode",false)),
  useJetEtMode_(iConfig.getParameter<Int_t>("useJetEtMode")),
  trkAcceptedJet_(iConfig.getUntrackedParameter<bool>("trkAcceptedJet",false)),
  useSubLeadingJet_(iConfig.getUntrackedParameter<bool>("useSubLeadingJet",false)),
  jetTrkOnly_(iConfig.getUntrackedParameter<bool>("jetTrkOnly",false)),
  fiducialCut_(iConfig.getUntrackedParameter<bool>("fiducialCut",false)),
  useQaulityStr_(iConfig.getUntrackedParameter<bool>("useQaulityStr")),
  qualityString_(iConfig.getUntrackedParameter<std::string>("qualityString")),
  usePxlPair_(iConfig.getUntrackedParameter<bool>("usePxlPair",false)),
  centrality_(0)
{

  histograms = new HiTrkEffHistograms(iConfig);

}


// ------------ method called for each event  ------------
void
HiTrkEffAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::ESHandle<TrackAssociatorBase> theAssociator;
  const TrackAssociatorByHits * theAssociatorByHits;
  edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
  edm::Handle<reco::RecoToSimCollection > recotosimCollectionH;
  edm::Handle<TrackingParticleCollection>  TPCollectionHeff, TPCollectionHfake;
  reco::RecoToSimCollection recSimColl;
  reco::SimToRecoCollection simRecColl;


  // reco track, vertex, beamspot collections

  edm::Handle<edm::View<reco::Track> >  trackCollection;
  iEvent.getByLabel(trackTags_,trackCollection);

  iEvent.getByLabel(vtxTags_,vertexCollectionH);
  iEvent.getByLabel(bsTags_,beamSpotH);

  LogDebug("HiTrkEffAnalyzer") <<" number of rec tracks = "<<trackCollection->size()<<std::endl;

  // Centrality information ----------------------
  double pixelMult = 0.0;

  centrality_ = new CentralityProvider(iSetup);
  centrality_->newEvent(iEvent,iSetup);
  int cbin = centrality_->getBin();

  if(pixelMultMode_){
     pixelMult = centrality_->raw()->multiplicityPixel();
     pixelMult = pixelMult/100.; // scale it (120K -> 1200)
  } 

  // PAT jet, to get leading jet ET

  float jet_et = 0.0, jet_eta=-999.0, jet_phi=0.0;
  std::vector<const pat::Jet *> sortedJets;         // jets for event classfication

  if(useJetEtMode_>0){
     edm::Handle<std::vector<pat::Jet> > jets;
     iEvent.getByLabel(jetTags_, jets);
     

     for(unsigned it=0; it<jets->size(); ++it){
	const pat::Jet* jet = &((*jets)[it]);

	if(trkAcceptedJet_) { // fill the jet pull only when the jet axes are within trk acceptance
	   if(fabs(jet->eta())<2.) {
	      sortedJets.push_back(jet);
	      sortByEtRef (&sortedJets);
	   }
	}else{
	   sortedJets.push_back(jet);
	   sortByEtRef (&sortedJets);
	}
     }
     
     for(unsigned it=0; it<sortedJets.size(); ++it){
	if(useSubLeadingJet_){ // use sub-leading jet
	   if(sortedJets.size()>1) it++;
	   else break; // if not sub-leading jet, break
	}
	jet_et = sortedJets[it]->et();
	jet_eta = sortedJets[it]->eta();
	jet_phi= sortedJets[it]->phi();
	break;
     }
  }

  // sim track collections
  float occHandle = 0.0;

  if(pixelMultMode_) occHandle = (float) pixelMult;
  else occHandle = jet_et;


  if(hasSimInfo_) {

    // sim track collections
    
    iEvent.getByLabel(label_tp_effic_,TPCollectionHeff);
    iEvent.getByLabel(label_tp_fake_,TPCollectionHfake);

    LogDebug("HiTrkEffAnalyzer") <<" number of sim tracks (for eff) = "<<TPCollectionHeff->size()
				 <<" number of sim tracks (for fake) = "<<TPCollectionHfake->size()<<std::endl;
    
    
    // association map generated on-the-fly or read from file
    
    if(doAssociation_){
      iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
      theAssociatorByHits = (const TrackAssociatorByHits*) theAssociator.product();
      
      simRecColl= theAssociatorByHits->associateSimToReco(trackCollection,TPCollectionHeff,&iEvent);
      recSimColl= theAssociatorByHits->associateRecoToSim(trackCollection,TPCollectionHfake,&iEvent);
      LogDebug("HiTrkEffAnalyzer") <<" association is done on the fly! \n"
				   <<" sim for eff = ["<<label_tp_effic_<<" ] \n"
				   <<" sim for fake = ["<<label_tp_fake_<<" ] \n";
    }else{
      iEvent.getByLabel(associatorMap_,simtorecoCollectionH);
      simRecColl= *(simtorecoCollectionH.product());
      
      iEvent.getByLabel(associatorMap_,recotosimCollectionH);
      recSimColl= *(recotosimCollectionH.product());
      LogDebug("HiTrkEffAnalyzer") <<" association is from the following association map \n"
                                   <<" map = ["<<associatorMap_<<" ] \n";
    }
    
    // -------------------- SIM loop ----------------------------------------
    
    for(TrackingParticleCollection::size_type i=0; i<TPCollectionHeff->size(); i++) {
      
      TrackingParticleRef tpr(TPCollectionHeff, i);
      TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
      
      if(tp->status() < 0 || tp->charge()==0) continue; //only charged primaries

      double drs = deltaR(jet_eta,jet_phi,tp->eta(),tp->phi());
      if(jetTrkOnly_ && drs>0.8) continue; // if jetTrkOnly_, only those within dR = 0.8;
      
      std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
      const reco::Track* mtr=0;
      size_t nrec=0;

      
      if(simRecColl.find(tpr) != simRecColl.end()){
	rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[tpr];
	nrec=rt.size();   

	if(useQaulityStr_){ // if true, re-calculate nrec 
	   nrec = 0; // set it to 0
	   for (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator rtit = rt.begin(); rtit != rt.end(); ++rtit){
	      const reco::Track* tmtr = rtit->first.get();
	      if(tmtr->quality(reco::TrackBase::qualityByName(qualityString_))) nrec++;
	   }
	}

	if(nrec) mtr = rt.begin()->first.get();      
      }
      
      SimTrack_t s = setSimTrack(*tp, *mtr, nrec, occHandle, cbin, sortedJets);
      histograms->fillSimHistograms(s);  
      
#ifdef DEBUG
      if(nrec) edm::LogVerbatim("HiTrkEffAnalyzer") 
	<< "TrackingParticle #" << i << " with pt=" << tp->pt() 
	<< " associated with quality:" << rt.begin()->second;
#endif
    }
    

  } //end if(hasSimInfo_)
    
  // -------------------- RECO loop ----------------------------------------

  for(edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); ++i){
    
    edm::RefToBase<reco::Track> track(trackCollection, i);
    reco::Track* tr=const_cast<reco::Track*>(track.get());
    
    if(fiducialCut_ && hitDeadPXF(*tr)) continue; // if track hits the dead region, igonore it;
    if(useQaulityStr_ && !tr->quality(reco::TrackBase::qualityByName(qualityString_))) continue;

    std::vector<std::pair<TrackingParticleRef, double> > tp;
    const TrackingParticle *mtp=0;
    size_t nsim=0;

    double drr = deltaR(jet_eta,jet_phi,track->eta(),track->phi());
    if(jetTrkOnly_ && drr>0.8) continue; // if jetTrkOnly_, only those within dR = 0.8;  

    if(hasSimInfo_ && recSimColl.find(track) != recSimColl.end()){
      tp = recSimColl[track];
      nsim=tp.size();
      if(nsim) mtp = tp.begin()->first.get();       
    }

    RecTrack_t r = setRecTrack(*tr, *mtp, nsim, occHandle, cbin,sortedJets);
    histograms->fillRecHistograms(r); 

#ifdef DEBUG
    if(nsim) edm::LogVerbatim("HiTrkEffAnalyzer") 
      << "reco::Track #" << i << " with pt=" << track->pt() 
      << " associated with quality:" << tp.begin()->second;
#endif
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
HiTrkEffAnalyzer::beginJob()
{

  histograms->declareHistograms();

}

// ------------ method called at the beginning of each new run  ------------
void
HiTrkEffAnalyzer::beginRun(const edm::Run& r, const edm::EventSetup& es)
{

  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiTrkEffAnalyzer::endJob()
{

  histograms->writeHistograms();

}

// ------------
SimTrack_t 
HiTrkEffAnalyzer::setSimTrack(TrackingParticle& tp, const reco::Track& mtr, size_t nrec, float jet, int cent, std::vector<const pat::Jet *> & sortedJets)
{
   
  SimTrack_t s;
  s.ids = tp.pdgId();
  s.etas = tp.eta();
  s.pts = tp.pt();
  s.hits = tp.matchedHit();
  s.status = tp.status();
  std::pair<bool,bool> acc = isAccepted(tp);
  if(usePxlPair_) s.acc = acc.first || acc.second; // If pixel pair seeding is sued
  else s.acc = acc.first; // for HI tracking, only triplet should be taken 

#ifdef DEBUG
  edm::LogVerbatim("HiTrkEffAnalyzer")  
    << "primary simtrack pt=" << s.pts 
    << " eta=" << s.etas
    << " hits=" << s.hits
    << " pdgId=" << s.ids;
#endif

  // remove the association if the track hits the bed region in FPIX
  // nrec>0 since we don't need it for nrec=0 case 
  if(fiducialCut_ && nrec>0 && hitDeadPXF(*const_cast<reco::Track*>(&mtr))) 
     nrec=0;

  s.nrec = nrec;
  s.jetr = jet;

  // if do closest jet, equivalent to closestJet_ mode in spectra ana
  if (useJetEtMode_==2) {
    Float_t bestJetDR=99, dR=99;
    Int_t bestJetInd=-99;
    unsigned int maxnjet
       = (sortedJets.size()<2) ?  sortedJets.size() : 2;

    for (UInt_t j=0; j<maxnjet; ++j) { // leading, sub-leading only, if there's any
      if (sortedJets[j]->et()<40) continue; // fake jet meaningless
      dR=deltaR(*sortedJets[j],tp);
      if(dR<0.8 && dR<bestJetDR){ // dR>0.8, should not influence efficiency..                                                                             
	 bestJetDR=dR, bestJetInd=j;
      }
    }
    if (bestJetInd<0) {
      s.jetr=0; s.jetar=-99;
    } else {
      s.jetr = sortedJets[bestJetInd]->et();
      s.jetar = sortedJets[bestJetInd]->eta();
      if (sortedJets[bestJetInd]->genParton()) {
	s.jrflavor = sortedJets[bestJetInd]->genParton()->pdgId();
      }
    }
    s.jrdr = bestJetDR;
    s.jrind = bestJetInd;
  }
  
  if(nrec > 0) {
    double dxy=0.0, dz=0.0, d0err=0.0, dzerr=0.0;
    testVertex(*const_cast<reco::Track*>(&mtr),dxy,dz,d0err,dzerr);
    s.ptr = mtr.pt();
    s.d0 = dxy;
    s.dz = dz;
    s.d0err = d0err;
    s.dzerr = dzerr;
    s.pterr = mtr.ptError();
    s.hitr = mtr.numberOfValidHits();
    s.algo = mtr.algo();
  } else {
    s.ptr = 0.0;
    s.dz = 0.0;
    s.d0 = 0.0;
    s.pterr = 0.0;
    s.d0err = 0.0;
    s.dzerr = 0.0;
    s.hitr = 0;
    s.algo = 0;
  }

  s.cbin = cent;

  return s;

}

// ------------
RecTrack_t 
HiTrkEffAnalyzer::setRecTrack(reco::Track& tr, const TrackingParticle& tp, size_t nsim, float jet, int cent, std::vector<const pat::Jet *> & sortedJets)
{

  RecTrack_t r;
  r.charge = tr.charge();
  r.etar = tr.eta();
  r.ptr = tr.pt();
  r.phir = tr.phi();

  double dxy=0.0, dz=0.0, d0err=0.0, dzerr=0.0;
  if(!testVertex(tr,dxy,dz,d0err,dzerr)) 
    edm::LogWarning("HiTrkEffAnalyzer") << "used the beamspot as there is no reco::Vertex";
  r.d0 = dxy;
  r.dz = dz;
  r.d0err = d0err;
  r.dzerr = dzerr;
  
  r.pterr = tr.ptError();
  r.hitr = tr.numberOfValidHits();
  r.algo = tr.algo();

#ifdef DEBUG
  edm::LogVerbatim("HiTrkEffAnalyzer")  
    << "reco track pt=" << r.ptr
    << " eta=" << r.etar
    << " hits=" << r.hitr;
#endif

  r.nsim = nsim;
  r.jetr = jet;
  // if do closest jet
  if (useJetEtMode_==2) {
    Float_t bestJetDR=99, dR=99;
    Int_t bestJetInd=-99;
    unsigned int maxnjet
       = (sortedJets.size()<2) ?  sortedJets.size() : 2;
    for (UInt_t j=0; j<maxnjet; ++j) {
      if (sortedJets[j]->et()<40) continue; // fake jet meaningless
      dR=deltaR(*sortedJets[j],tr);
      if(dR<0.8 && dR<bestJetDR){ // dR>0.8, should not influence efficiency..
         bestJetDR=dR, bestJetInd=j;
      }
    }
    if (bestJetInd<0) {
      r.jetr=0; r.jetar=-99;
    } else {
      r.jetr = sortedJets[bestJetInd]->et();
      r.jetar = sortedJets[bestJetInd]->eta();
      if (sortedJets[bestJetInd]->genParton()) {
	r.jrflavor = abs(sortedJets[bestJetInd]->genParton()->pdgId());
      }
    }
    r.jrdr = bestJetDR;
    r.jrind = bestJetInd;
  }

  if(nsim>0) {
    r.status = tp.status();
    r.ids = tp.pdgId();
    
    int parentId=0;
    if(tp.parentVertex()->nSourceTracks() > 0) {
      parentId = (*(tp.parentVertex()->sourceTracks_begin()))->pdgId();
      edm::LogVerbatim("HiTrkEffAnalyzer") 
    	<< "pdg Id = " << tp.pdgId()
    	<< " parent Id = " << parentId;
    }
    r.parids = parentId;
    
    r.etas = tp.eta();
    r.pts = tp.pt();
  } else {
    r.status = 0;
    r.ids = 0;
    r.parids = 0;
    r.etas = 0.0;
    r.pts = 0.0;
  }
  
  r.cbin = cent;
  
  return r;

}

// ------------
bool
HiTrkEffAnalyzer::testVertex(reco::Track& tr, double &dxy, double &dz, double &d0err, double &dzerr)
{

  d0err = tr.d0Error();
  dzerr = tr.dzError();

  const reco::VertexCollection *vtxs = vertexCollectionH.product();

  if(vtxs->size() == 0 || vtxs->begin()->isFake()) {
    dxy = tr.dxy(beamSpotH->position());
    dz  = 0.0;
    d0err = sqrt ( (d0err*d0err) + (beamSpotH->BeamWidthX()*beamSpotH->BeamWidthY()) );
    dzerr = 0.0;
    return false;
  } else {
    dxy = tr.dxy(vtxs->begin()->position()); // FIX ME (use most populated vertex?)
    dz  = tr.dz(vtxs->begin()->position());  // FIX ME (or closest vertex to track?)
    d0err = sqrt ( (d0err*d0err) + (vtxs->begin()->xError()*vtxs->begin()->yError()) );
    dzerr = sqrt ( (dzerr*dzerr) + (vtxs->begin()->zError()*vtxs->begin()->zError()) );
    return true;
  }

}

// ------------
std::pair<bool,bool> 
HiTrkEffAnalyzer::isAccepted(TrackingParticle & tp)
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

  return std::pair<bool,bool>(canBeTriplet, canBePair);
}

// ------------
int 
HiTrkEffAnalyzer::getLayerId(const PSimHit & simHit)
{
  unsigned int id = simHit.detUnitId();

  if(theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelBarrel)
  {
    PXBDetId pid(id);
    return pid.layer() - 1; // 0, 1, 2
  }

  if(theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelEndcap)
  {
    PXFDetId pid(id);
    return BPix3 + ((pid.side()-1) << 1) + pid.disk(); // 3 -
  }

  // strip
  return -1;
}

// ---------------
bool
HiTrkEffAnalyzer::hitDeadPXF(reco::Track& tr){

   //-----------------------------------------------
   // For a given track, check whether this contains 
   // hits on the dead region in the forward pixel 
   //-----------------------------------------------

   bool hitDeadRegion = false;

   for(trackingRecHit_iterator recHit = tr.recHitsBegin();recHit!= tr.recHitsEnd(); recHit++){

      if((*recHit)->isValid()){

	 DetId detId = (*recHit)->geographicalId();
	 if(!theTracker->idToDet(detId)) continue;

	 Int_t diskLayerNum=0, bladeLayerNum=0, hcylLayerNum=0;
	 
	 unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());

	 if (subdetId == PixelSubdetector::PixelEndcap){
	    
	    PixelEndcapName pxfname(detId.rawId());
	    diskLayerNum = pxfname.diskName();
	    bladeLayerNum = pxfname.bladeName();
	    hcylLayerNum = pxfname.halfCylinder();
	    
	    // hard-coded now based on /UserCode/Appeltel/PixelFiducialRemover/pixelfiducialremover_cfg.py
	    if((bladeLayerNum==4 || bladeLayerNum==5 || bladeLayerNum==6) &&
	       (diskLayerNum==2) && (hcylLayerNum==4)) hitDeadRegion = true;
	 }
	 
      }// end of isValid
   }

   return hitDeadRegion;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HiTrkEffAnalyzer);
