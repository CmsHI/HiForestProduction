#ifndef PIXELTRACKLETANA_FINDER
#define PIXELTRACKLETANA_FINDER

#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MitHig/PixelTracklet/interface/TrackletCorrections.h"
#include "MitHig/PixelTracklet/interface/Tracklet.h"
#include "TNtuple.h"

class TrackletFinder{
public :

  TrackletFinder(
		 TrackletCorrections* corr, const TrackerGeometry* geo, bool count = false, bool checkSecondLayer = false
		 ): 
  corr_(corr), 
  geo_(geo), 
  countMode_(count), 
  checkSecondLayer_(checkSecondLayer) 
  {;}
 int getNHits1(){ return layer1_.size(); }
 int getNHits2(){ return layer2_.size(); }

 void sortLayers();
  std::vector<Tracklet> makeTracklets(bool invert = false);
  std::vector<Tracklet> cleanTracklets(std::vector<Tracklet> input);
  std::vector<Tracklet> getTracklets() {return cleanTracklets(makeTracklets());}
  void fillPixelEvent(TNtuple* nt1,TNtuple* nt2);
  void setEvent(const edm::Event & iEvent);
  void setVertex(math::XYZVector vtx) { vertex_ = vtx; }
  //  bool compareTracklet(Tracklet a,Tracklet b);
  
private:

  const edm::Event* event_;
  std::vector<const SiPixelRecHit*> layer1_;
  std::vector<const SiPixelRecHit*> layer2_;
  math::XYZVector vertex_;
  TrackletCorrections* corr_;  
  const TrackerGeometry* geo_;
  bool countMode_;
  bool checkSecondLayer_;

};

//bool compareTracklet(Tracklet a,Tracklet b) {    return fabs(a.dR2())<fabs(b.dR2()); }

#endif
