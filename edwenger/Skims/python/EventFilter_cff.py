import FWCore.ParameterSet.Config as cms

from edwenger.Skims.eventSelection_cff import *
from edwenger.Skims.hfCoincFilter_cff import *
from edwenger.Skims.Analyzers_cff import *

eventFilter = cms.Sequence(preTrgTest *
                           minBiasBscFilter *    ## BSC OR, !BSCHALO
                           postTrgTest *
                           hfCoincFilter *       ## E=3 GeV tower threshold
                           purityFractionFilter)

eventFilterLoose = cms.Sequence(preTrgTest *
                                minBiasBscFilter *    ## BSC OR, !BSCHALO
                                postTrgTest *
                                hfEitherOfTwoFilter *       ## E=3 GeV tower threshold
                                purityFractionFilter)


eventFilterMC = cms.Sequence(preTrgTest *
                             minBiasBscFilterMC *    ## BSC OR, !BSCHALO
                             postTrgTest *
                             hfCoincFilter *         ## E=3 GeV tower threshold
                             purityFractionFilter)
