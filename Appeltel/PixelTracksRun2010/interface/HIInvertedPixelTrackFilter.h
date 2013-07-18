#ifndef _HIInvertedPixelTrackFilter_h_
#define _HIInvertedPixelTrackFilter_h_

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeTrackFilter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace edm { class ParameterSet; class EventSetup; class Event;}

class HIInvertedPixelTrackFilter : public ClusterShapeTrackFilter {
public:
	HIInvertedPixelTrackFilter(const edm::ParameterSet& ps, const edm::EventSetup& es);
	virtual ~HIInvertedPixelTrackFilter();
	virtual bool operator() (const reco::Track*, const PixelTrackFilter::Hits & hits) const;
	virtual void update(edm::Event& ev);
private:
	double theTIPMin, theNSigmaTipMinTolerance;
	double theLIPMin, theNSigmaLipMinTolerance;
	double theChi2Min, thePtMin;
	bool useClusterShape;
	edm::InputTag theVertexCollection; 	
	const reco::VertexCollection *theVertices;

};

#endif



