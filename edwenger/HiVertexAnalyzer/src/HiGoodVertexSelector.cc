#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// define your producer name
typedef SingleObjectSelector<
  reco::VertexCollection,
  StringCutObjectSelector<reco::Vertex> 
  > HiGoodVertexSelection;

// declare the module as plugin
DEFINE_FWK_MODULE( HiGoodVertexSelection );
