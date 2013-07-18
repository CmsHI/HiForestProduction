#ifndef NewVertexProducer_H
#define NewVertexProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TVectorD.h"

namespace edm { class Event; class EventSetup; }

// (position,error2), list of track indexes
typedef std::pair<std::pair<double,double>, std::vector<int> > NewVertex;

typedef std::vector<NewVertex> NewVertexCollection;

class NewVertexProducer : public edm::EDProducer
{
public:
  explicit NewVertexProducer(const edm::ParameterSet& ps);
  ~NewVertexProducer();
  virtual void produce(edm::Event& ev, const edm::EventSetup& es);
 
private:
  double getLambda(const NewVertexCollection & vertices);

  void numberOfVertices
    (const std::vector<std::pair<double,double> > & points,
     const std::vector<std::pair<TVectorD,TVectorD> > & clusters,
     NewVertexCollection & rec);


  edm::ParameterSet theConfig;

  double thePtMin;
  double nSigma;

  int    theMethod;

  double dMax; 
  double theProbMin;

  unsigned int theTrackMin;
};
#endif
