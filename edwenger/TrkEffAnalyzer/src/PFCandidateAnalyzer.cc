#include "edwenger/TrkEffAnalyzer/interface/PFCandidateAnalyzer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"

//EAW
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "edwenger/VertexAnalyzer/interface/VertexComparator.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <Math/DistFunc.h>
#include "TMath.h"
//


#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"


using namespace std;
using namespace edm;
using namespace reco;


PFCandidateAnalyzer::PFCandidateAnalyzer(const edm::ParameterSet& iConfig) {
  


  inputTagPFCandidates_ = iConfig.getParameter<InputTag>("PFCandidates");
  inputTagVertices_ = iConfig.getParameter<InputTag>("Vertices");
  inputTagSimTracks_ = iConfig.getParameter<InputTag>("SimTracks");
  inputTagTracks_ = iConfig.getParameter<InputTag>("Tracks");

  verbose_ = iConfig.getUntrackedParameter<bool>("verbose",false);
  printBlocks_ = iConfig.getUntrackedParameter<bool>("printBlocks",false);
  thePtMin_ = iConfig.getUntrackedParameter<double>("ptMin",3.0);
  hasSimInfo_ = iConfig.getUntrackedParameter<bool>("hasSimInfo");

  minHits_ = iConfig.getUntrackedParameter<double>("minHits");
  maxPtErr_ = iConfig.getUntrackedParameter<double>("maxPtErr");
  maxD0_ = iConfig.getUntrackedParameter<double>("maxD0");
  maxDZ_ = iConfig.getUntrackedParameter<double>("maxDZ");
  maxD0Norm_ = iConfig.getUntrackedParameter<double>("maxD0Norm");
  maxDZNorm_ = iConfig.getUntrackedParameter<double>("maxDZNorm");
  pixelSeedOnly_ = iConfig.getUntrackedParameter<bool>("pixelSeedOnly");

  LogDebug("PFCandidateAnalyzer")
    <<" input collection : "<<inputTagPFCandidates_ ;



   
}



PFCandidateAnalyzer::~PFCandidateAnalyzer() { }



void 
PFCandidateAnalyzer::beginRun(const edm::Run& run, 
			      const edm::EventSetup & es) {}

void 
PFCandidateAnalyzer::beginJob() {

  nt = f->make<TNtuple>("nt","PF Testing","type:pt:tkptmax:tkptsum:eetmax:eetsum:hetmax:hetsum:nhits:relpterr:algo:nd0:ndz:fake");

}


void 
PFCandidateAnalyzer::analyze(const Event& iEvent, 
			     const EventSetup& iSetup) {
  
  LogDebug("PFCandidateAnalyzer")<<"START event: "<<iEvent.id().event()
			 <<" in run "<<iEvent.id().run()<<endl;

  
  // get vertices

  edm::Handle<reco::VertexCollection> vtxsH;
  iEvent.getByLabel(inputTagVertices_,vtxsH);
  reco::VertexCollection vtxs = *vtxsH;
  
  // do reco-to-sim association

  edm::Handle<TrackingParticleCollection>  TPCollectionHfake;
  edm::Handle<edm::View<reco::Track> >  trackCollection;
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  const TrackAssociatorByHits *theAssociatorByHits;
  reco::RecoToSimCollection recSimColl;
  
  if(hasSimInfo_) {
    iEvent.getByLabel(inputTagSimTracks_,TPCollectionHfake);
    iEvent.getByLabel(inputTagTracks_,trackCollection);
    iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
    theAssociatorByHits = (const TrackAssociatorByHits*) theAssociator.product();  
    recSimColl= theAssociatorByHits->associateRecoToSim(trackCollection,TPCollectionHfake,&iEvent);
  }

  // get PFCandidates

  Handle<PFCandidateCollection> pfCandidates;
  iEvent.getByLabel(inputTagPFCandidates_, pfCandidates);

  for( unsigned i=0; i<pfCandidates->size(); i++ ) {
    
    max_nhits=0.0;
    max_relpterr=0.0;
    max_algo=0.0;
    max_nd0=0.0;
    max_ndz=0.0;
    max_fake=0.0;

    const reco::PFCandidate& cand = (*pfCandidates)[i];
    
    cand_type = cand.particleId();
    cand_pt = cand.pt();

    //--------- EDIT: Apr 21, 2010 (EAW)
    if(cand_pt < thePtMin_) continue;
    
    if(!(cand_type == PFCandidate::h ||     //type1
	 cand_type == PFCandidate::e ||     //type2
	 cand_type == PFCandidate::mu ||     //type3
	 cand_type == PFCandidate::gamma || //type4
	 cand_type == PFCandidate::h0      //type5
	 )
       ) continue;

    if( verbose_ ) {
      cout<<cand<<endl;
      if (printBlocks_) printElementsInBlocks(cand);
    }   

    for(unsigned i=0; i<cand.elementsInBlocks().size(); i++) {

      PFBlockRef blockRef = cand.elementsInBlocks()[i].first;

      // get a copy of the link data
      PFBlock::LinkData linkData =  blockRef->linkData();
     
      if(i && verbose_) {
	cout << "WARNING: more than one block in this candidate" << endl; 
	continue;
      }

      cout << "\nCandidate Type #" << cand_type << endl;
      const edm::OwnVector< reco::PFBlockElement >& elements = (*blockRef).elements();
      
      max_trk=0.0, sum_trk=0.0, max_ecal=0.0, sum_ecal=0.0, max_hcal=0.0, sum_hcal=0.0;

      for(unsigned ie=0; ie<elements.size(); ie++) {

	PFBlockElement::Type type = elements[ie].type();
	if(type==PFBlockElement::TRACK) {
	  cout << "TRK:";
	  reco::TrackRef trackRef = elements[ie].trackRef();

	  //----- track quality selections

	  double nhits = trackRef->numberOfValidHits();
	  if(nhits<minHits_) continue;

	  double relpterr = trackRef->ptError()/trackRef->pt();
	  if(relpterr > maxPtErr_) continue;

	  double algo = trackRef->algo();
	  if(algo > 7 && pixelSeedOnly_) continue;

	  double d0 = trackRef->dxy(vtxs[0].position());
	  double dz = trackRef->dz(vtxs[0].position());

	  if(fabs(d0) > maxD0_) continue;
	  if(fabs(dz) > maxDZ_) continue;

	  double d0err = sqrt(trackRef->d0Error()*trackRef->d0Error()+vtxs[0].xError()*vtxs[0].yError());
	  double dzerr = sqrt(trackRef->dzError()*trackRef->dzError()+vtxs[0].zError()*vtxs[0].zError());

	  if(fabs(d0/d0err) > maxD0Norm_) continue;
	  if(fabs(dz/dzerr) > maxDZNorm_) continue;

	  bool fake=false;

	  if(hasSimInfo_)
	    if(recSimColl.find(edm::RefToBase<reco::Track>(trackRef)) == recSimColl.end())
	      fake=true;

	  //-----

	  double trkpt = trackRef->pt();
	  cout << "pt=" << trkpt << endl;
	  sum_trk+=trkpt;
	  if(trkpt>max_trk) {
	    max_trk=trkpt;
	    max_nhits=nhits;
	    max_relpterr=relpterr;
	    max_algo=algo;
	    max_nd0=d0/d0err;
	    max_ndz=dz/dzerr;
	    max_fake=fake;
	  }

	  // look for associated elements of all types
	  // all types of links are considered. 
	  // the elements are sorted by increasing distance
	  std::multimap<double, unsigned> ecalElems;
	  blockRef->associatedElements( ie,  linkData,
					ecalElems ,
					reco::PFBlockElement::ECAL,
					reco::PFBlock::LINKTEST_ALL );
	  
	  std::multimap<double, unsigned> hcalElems;
	  blockRef->associatedElements( ie,  linkData,
					hcalElems,
					reco::PFBlockElement::HCAL,
					reco::PFBlock::LINKTEST_ALL );

	  cout << "\t" << ecalElems.size() << " ecalElems" << endl;

	  std::multimap <double, unsigned> :: const_iterator  iee;
	  for(iee = ecalElems.begin(); iee != ecalElems.end(); ++iee ) {

	    unsigned index = iee->second;
	    double dist = iee->first * 1000;

	    cout << "\t\t ECAL element #" << index 
		 << " at distance x 1000 = " << dist << endl;
	  }

	  cout << "\t" << hcalElems.size() << " hcalElems" << endl;

	  std::multimap <double, unsigned> :: const_iterator  ihe;
	  for(ihe = hcalElems.begin(); ihe != hcalElems.end(); ++ihe ) {

	    unsigned index = ihe->second;
	    double dist = ihe->first * 1000;

	    cout << "\t\t HCAL element #" << index 
		 << " at distance x 1000 =  " << dist << endl;
	  }

	} 
	else if(type==PFBlockElement::ECAL) {
	  cout << "ECAL:";
	  PFClusterRef ecalRef = elements[ie].clusterRef();
	  double eet = ecalRef->energy()/cosh(ecalRef->eta());
	  cout << "et=" << eet << endl;
	  sum_ecal+=eet;
	  if(eet>max_ecal) max_ecal=eet;
	} 
	else if(type==PFBlockElement::HCAL) {
	  cout << "HCAL:";
	  PFClusterRef hcalRef = elements[ie].clusterRef();
	  double het = hcalRef->energy()/cosh(hcalRef->eta());
	  cout << "et=" << het << endl;
	  sum_hcal+=het;
	  if(het>max_hcal) max_hcal=het;
	}

      } //end loop over elements
 
    } //end loop over blocks

    nt->Fill(cand_type,cand_pt,max_trk,sum_trk,max_ecal,sum_ecal,max_hcal,sum_hcal,max_nhits,max_relpterr,max_algo,max_nd0,max_ndz,max_fake);

    //---------
  }
   
 
  LogDebug("PFCandidateAnalyzer")<<"STOP event: "<<iEvent.id().event()
			 <<" in run "<<iEvent.id().run()<<endl;
}




void PFCandidateAnalyzer::printElementsInBlocks(const PFCandidate& cand,
						ostream& out) const {
  if(!out) return;

  PFBlockRef firstRef;

  assert(!cand.elementsInBlocks().empty() );
  for(unsigned i=0; i<cand.elementsInBlocks().size(); i++) {
    PFBlockRef blockRef = cand.elementsInBlocks()[i].first;

    if(blockRef.isNull()) {
      cerr<<"ERROR! no block ref!";
      continue;
    }

    if(!i) {
      out<<(*blockRef);
      firstRef = blockRef;
    }
    else if( blockRef!=firstRef) {
      cerr<<"WARNING! This PFCandidate is not made from a single block"<<endl;
    }
 
    out<<"\t"<<cand.elementsInBlocks()[i].second<<endl;
  }
}


DEFINE_FWK_MODULE(PFCandidateAnalyzer);
