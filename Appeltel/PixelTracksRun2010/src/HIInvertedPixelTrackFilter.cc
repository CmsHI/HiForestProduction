#include "Appeltel/PixelTracksRun2010/interface/HIInvertedPixelTrackFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

using namespace std;
using namespace edm;

/*****************************************************************************/
HIInvertedPixelTrackFilter::HIInvertedPixelTrackFilter (const edm::ParameterSet& ps, const edm::EventSetup& es) :
ClusterShapeTrackFilter(ps,es),
theTIPMin( ps.getParameter<double>("tipMin") ),
theNSigmaTipMinTolerance( ps.getParameter<double>("nSigmaTipMinTolerance")),
theLIPMin( ps.getParameter<double>("lipMin") ),
theNSigmaLipMinTolerance( ps.getParameter<double>("nSigmaLipMinTolerance")),
theChi2Min( ps.getParameter<double>("chi2") ),
thePtMin( ps.getParameter<double>("ptMin") ),
useClusterShape( ps.getParameter<bool>("useClusterShape") ),
theVertexCollection( ps.getParameter<edm::InputTag>("VertexCollection")),
theVertices(0)
{ 
}

/*****************************************************************************/
HIInvertedPixelTrackFilter::~HIInvertedPixelTrackFilter()
{ }

/*****************************************************************************/
bool HIInvertedPixelTrackFilter::operator() (const reco::Track* track,const PixelTrackFilter::Hits & recHits) const
{

  if (!track) return false; 
  if (track->chi2() < theChi2Min || track->pt() < thePtMin) return false; 
  
  
  math::XYZPoint vtxPoint(0.0,0.0,0.0);
  double vzErr =0.0, vxErr=0.0, vyErr=0.0;
  
  if(theVertices->size()>0) {
    vtxPoint=theVertices->begin()->position();
    vzErr=theVertices->begin()->zError();
    vxErr=theVertices->begin()->xError();
    vyErr=theVertices->begin()->yError();
  } else {
    // THINK OF SOMETHING TO DO IF THERE IS NO VERTEX
  }
  
  double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
  d0 = -1.*track->dxy(vtxPoint);
  dz = track->dz(vtxPoint);
  d0sigma = sqrt(track->d0Error()*track->d0Error()+vxErr*vyErr);
  dzsigma = sqrt(track->dzError()*track->dzError()+vzErr*vzErr);
  
  if (fabs(d0)<theTIPMin) return false;
  if ((fabs(d0)/d0sigma)<theNSigmaTipMinTolerance) return false;
  if (fabs(dz)<theLIPMin) return false;
  if ((fabs(dz)/dzsigma)<theNSigmaLipMinTolerance) return false;
  
  bool ok = true;
  if(useClusterShape) ok = !ClusterShapeTrackFilter::operator() (track,recHits);
  
  return ok;
}

/*****************************************************************************/
void HIInvertedPixelTrackFilter::update(edm::Event& ev)
{
  
  // Get reco vertices
  edm::Handle<reco::VertexCollection> vc;
  ev.getByLabel(theVertexCollection, vc);
  theVertices = vc.product();
  
  if(theVertices->size()>0) {
    LogInfo("HeavyIonVertexing") 
      << "[HIPixelTrackFilter] Pixel track selection based on best vertex"
      << "\n   vz = " << theVertices->begin()->z()  
      << "\n   vz sigma = " << theVertices->begin()->zError();
  } else {
    LogError("HeavyIonVertexing") // this can be made a warning when operator() is fixed
      << "No vertex found in collection '" << theVertexCollection << "'";
  }

  return;
  
}

