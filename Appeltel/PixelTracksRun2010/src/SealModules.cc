#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "Appeltel/PixelTracksRun2010/interface/HiTrackListMerger.h"


using cms::HiTrackListMerger;



DEFINE_FWK_MODULE(HiTrackListMerger);

// Inverted Pixel track filter
#include "Appeltel/PixelTracksRun2010/interface/HIInvertedPixelTrackFilter.h"
#include "RecoPixelVertexing/PixelTrackFitting/interface/PixelTrackFilterFactory.h"
DEFINE_EDM_PLUGIN(PixelTrackFilterWithESFactory, HIInvertedPixelTrackFilter, "HIInvertedPixelTrackFilter");


