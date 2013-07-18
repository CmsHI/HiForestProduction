#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace reco
{
   namespace modules
   {
      
      // define your producer name
      typedef ObjectSelector<
	      SortCollectionSelector<
	      reco::VertexCollection,
	      VertexHigherPtSquared
	      > 
	      > SumPtOrdVertexSelector;
      
      // declare the module as plugin
      DEFINE_FWK_MODULE( SumPtOrdVertexSelector );
   }
}
