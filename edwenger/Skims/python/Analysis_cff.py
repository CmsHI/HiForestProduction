import FWCore.ParameterSet.Config as cms

from edwenger.Skims.patAna_cff import *   # PAT jet sequence
from edwenger.Skims.eventSelection_cff import *
from edwenger.Skims.Analyzers_cff import *

analysisSeq = cms.Sequence(patAnaSequence *
			   cleanedPatJets *
                           selectedVertex *
                           (postExtraTrkVtxSel + postExtraTrkVtxSel_jet60 + postExtraTrkVtxSel_jet120) * # trk,vtx,jet 
                           (looseTrackAna + looseTrackAna_STD) *
                           loosetrkEffAnalyzer *
                           primaryVertexFilter *
                           (trackAna + trackAna_STD) *
                           trkEffAnalyzer)

analysisGEN = cms.Sequence(preTrackAna)
