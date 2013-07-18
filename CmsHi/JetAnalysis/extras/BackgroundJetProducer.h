#ifndef RecoJets_JetProducers_plugins_BackgroundJetProducer_h
#define RecoJets_JetProducers_plugins_BackgroundJetProducer_h

#include "RecoJets/JetProducers/plugins/VirtualJetProducer.h"

class BackgroundJetProducer : public VirtualJetProducer
{

public:
  //
  // construction/destruction
  //
  explicit BackgroundJetProducer(const edm::ParameterSet& iConfig);
  virtual ~BackgroundJetProducer();

protected:

  //
  // member functions
  //

  virtual void runAlgorithm( edm::Event& iEvent, const edm::EventSetup& iSetup );
  virtual void output(  edm::Event & iEvent, edm::EventSetup const& iSetup );
  template< typename T >
  void writeBkgJets( edm::Event & iEvent, edm::EventSetup const& iSetup );

 private:

  // trackjet clustering parameters
  bool useOnlyVertexTracks_;
  bool useOnlyOnePV_;
  float dzTrVtxMax_;

  int nFill_;
  float etaMax_;
  bool avoidNegative_;

  const CaloGeometry *geo;
};


#endif
