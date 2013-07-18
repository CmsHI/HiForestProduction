#include "NewVertexProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "../interface/FastPairwiseNearestNeighbor.h"
#include "../interface/OptimalTree.h"
#include "../interface/GaussianMixture.h"
#include "../interface/KMeansMethod.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "TMath.h"

#undef Debug

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
NewVertexProducer::NewVertexProducer
  (const edm::ParameterSet& ps) : theConfig(ps)
{
  produces<reco::VertexCollection>();

  thePtMin   = theConfig.getParameter<double>("PtMin");
  nSigma     = theConfig.getParameter<double>("nSigma");

  theMethod  = theConfig.getParameter<int>("Method");

  dMax       = theConfig.getParameter<double>("dMax");
  theProbMin = theConfig.getParameter<double>("ProbMin");

  theTrackMin = (unsigned int)theConfig.getParameter<int>("TrackMin");

  LogTrace("NewVertices")
            << " [NewVertices] method = " << theMethod;
}

/*****************************************************************************/
NewVertexProducer::~NewVertexProducer()
{ 
}

/*****************************************************************************/
double NewVertexProducer::getLambda(const NewVertexCollection & vertices)
{
  int n = 0;
  for(NewVertexCollection::const_iterator vertex = vertices.begin();
                                          vertex!= vertices.end(); vertex++)
    n += vertex->second.size();

  double chi2 = 0.;
  for(NewVertexCollection::const_iterator vertex = vertices.begin();
                                          vertex!= vertices.end(); vertex++)
  {
    int k = vertex->second.size();

    if(k > 0)
      chi2 += -2 * k * log(float(k)/n);
  }

  return chi2;
}

/*****************************************************************************/
void NewVertexProducer::numberOfVertices
  (const vector<pair<double,double> > & points,
   const vector<pair<TVectorD,TVectorD> > & clusters,
   NewVertexCollection & rec)
{
  rec.clear();

  for(int K = 1; K <= int(points.size()); K++)
  {
    // Get results from hierarchical clustering
    TVectorD mu(K); mu = clusters[K].first;
    TVectorD P(K) ; P  = clusters[K].second;

    double chi2 = 0.;

    if(theMethod == 1)
    {
      KMeansMethod theKMeansMethod;
      chi2 = theKMeansMethod.run(K,points, mu,P, rec);
    }

    if(theMethod == 2)
    {
      GaussianMixture theGaussianMixture;
      chi2 = theGaussianMixture.run(K,points, mu,P, rec);
    }

    // Stopping condition
    if(chi2 < getLambda(rec) ||
       TMath::Prob(chi2 - getLambda(rec), points.size()) > theProbMin)
      break;
  }
}


/*****************************************************************************/
void NewVertexProducer::produce
  (edm::Event& ev, const edm::EventSetup& es)
{
  // Get beamSpot
  edm::Handle<reco::BeamSpot>      beamSpotHandle;
  ev.getByLabel("offlineBeamSpot", beamSpotHandle);
  const reco::BeamSpot * theBeamSpot = beamSpotHandle.product();

  LogTrace("NewVertices")
    << " [NewVertices] beamSpot at " << theBeamSpot->position();

  // Get tracks
  edm::Handle<reco::TrackCollection> trackCollection;
  string trackCollectionName =
    theConfig.getParameter<string>("TrackCollection");
  ev.getByLabel(trackCollectionName, trackCollection);
  const reco::TrackCollection recTracks = *(trackCollection.product());

  // Select tracks and create points
  vector<pair<double, double> > points;
  reco::TrackRefVector trks;

#ifdef Debug
  ofstream file("z.dat");
#endif

  for(unsigned int i=0; i<recTracks.size(); i++)
    if(recTracks[i].pt() > thePtMin)
    {
      // Transverse impact 
      double dt = fabs(recTracks[i].dxy(theBeamSpot->position()));
      double st = sqrt(sqr(recTracks[i].dxyError()) +
                       sqr(theBeamSpot->BeamWidthX()));

      if(dt < nSigma * st)
      {
        // z wrt beam spot
        pair<double,double> point;
        point.first  = recTracks[i].dz(theBeamSpot->position()) +
                                       theBeamSpot->position().z();

        // sigma_z^2 taking into account the beam width
        point.second =
          sqr(recTracks[i].dzError()) + 
          sqr(theBeamSpot->BeamWidthX()) *
          sqr(recTracks[i].pz() / recTracks[i].pt());

#ifdef Debug
        file << " " << recTracks[i].dz(theBeamSpot->position()) +
                                       theBeamSpot->position().z()
             << " " << recTracks[i].dzError() << endl;
#endif

        points.push_back(point);

        // Just technical
        trks.push_back( reco::TrackRef(trackCollection, i) );
      }
      else
         LogTrace("NewVertices")
            << " [NewVertices] not seleted track #" << i
            << " | dt/st = " << dt << " " << st << endl;
    }

#ifdef Debug
  file.close();
#endif

  LogTrace("NewVertices")
            << " [NewVertices] selected tracks: "
            << points.size() << " (out of " << recTracks.size()
            << ")"; 

  auto_ptr<reco::VertexCollection> vertices(new reco::VertexCollection);

  if(points.size() > 0)
  {
    // Initialize clusters, at most maxVertices
    unsigned int maxVertices = 100;
    vector<pair<TVectorD,TVectorD> > clusters;
    for(unsigned int i = 0; i <= maxVertices; i++)
    {
      TVectorD mu(i);
      TVectorD  P(i);

      clusters.push_back(pair<TVectorD,TVectorD>(mu,P));
    }

    //////////////////
    // fPNN
    unsigned int nOptimal;
    vector<vector<int> > lists;

    FastPairwiseNearestNeighbor thePairGroupMethod(dMax);
    thePairGroupMethod.run(points, clusters, nOptimal, lists, maxVertices);

    NewVertexCollection rec;

    //////////////////
    // OptimalTree 
    if(theMethod == 0)
    {
      int K = nOptimal;
      OptimalTree theOptimalTree;
      TVectorD mu(K); mu = clusters[K].first;
      TVectorD P(K) ; P  = clusters[K].second;

      theOptimalTree.run(K,points, mu,P, lists, rec);
    }
    else
      numberOfVertices(points, clusters, rec);


    // Look at vertices
    for(NewVertexCollection::const_iterator vertex = rec.begin();
                                            vertex!= rec.end(); vertex++)
    if(vertex->second.size() >= 1)
    {
      LogTrace("NewVertices")
        << " [NewVertices]  vertex with "
        << vertex->second.size() << " tracks, at "
        << vertex->first.first << " +/- "
        << vertex->first.second
        << " cm";

      // Position
      double pos  = vertex->first.first;
      double err2 = vertex->first.second;

      // Error
      reco::Vertex::Error err;
      err(2,2) = err2; 

      // Use transverse coordinates of the offline beam spot  
      reco::Vertex ver(reco::Vertex::Point(theBeamSpot->position().x(),
                                           theBeamSpot->position().y(),
                                           pos), err, 0, 1, 1);

      int sum_q = 0;
      reco::TrackBase::Vector sum_p(0.,0.,0.);

      for(vector<int>::const_iterator i = vertex->second.begin();
                                      i < vertex->second.end(); i++)
      {
        ver.add(reco::TrackBaseRef(trks[*i]));

        sum_q += trks[*i]->charge();
        sum_p += trks[*i]->momentum();
      }

      // Is it a looper?
      if(vertex->second.size() == 2 && sum_q == 0 && sum_p.R() < 0.050)
        LogTrace("NewVertices")
          << " [NewVertices]  looper, sum_p = " << sum_p.R() * 1e+3 << " MeV/c";
      else 
      {
        // enough tracks?
        if(vertex->second.size() >= theTrackMin ||
                recTracks.size() <  theTrackMin)
          vertices->push_back(ver); // store
      }
    }

    LogTrace("NewVertices")
      << " [NewVertices] found vertices: " << vertices->size();
  }
  ev.put(vertices);
}

