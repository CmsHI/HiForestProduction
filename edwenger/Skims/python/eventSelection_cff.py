import FWCore.ParameterSet.Config as cms

# ================ Physics Declared Bit ================
physDeclFilter = cms.EDFilter("PhysDecl",
    applyfilter = cms.untracked.bool(False),
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
    )


# ================ L1 Filters =====================
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
from L1Trigger.Skimmer.l1Filter_cfi import l1Filter

bptxAnd = hltLevel1GTSeed.clone(
    L1TechTriggerSeeding = cms.bool(True),
    L1SeedsLogicalExpression = cms.string('0')
    )

bscOr = hltLevel1GTSeed.clone(
    L1TechTriggerSeeding = cms.bool(True),
    L1SeedsLogicalExpression = cms.string('34 OR 40 OR 41')
    )

bscNoHalo = hltLevel1GTSeed.clone(
    L1TechTriggerSeeding = cms.bool(True),
    L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')
    )

bscOrBptxOr = l1Filter.clone(
    algorithms = ["L1_BscMinBiasOR_BptxPlusORMinus"]
    )


# ================ HLT Filters =====================
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel

hltMinBias = hltHighLevel.clone(
    HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus'),
    andOr = cms.bool(True)
    )

hltJets = hltHighLevel.clone(
    HLTPaths = cms.vstring('HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U'),
    andOr = cms.bool(True)
    )


# ================ Monster Event Rejection ==============
purityFractionFilter = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )


# ================ Vertex Filter ========================
primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(15),	
    maxd0 = cms.double(2)	
    )

# filter on events without vertex and produce vertex collection
# containing only the vertex with the most tracks
from edwenger.VertexAnalyzer.SelectedVertex_cfi import *

# select only vertices passing a string cut and sort them
from edwenger.VertexAnalyzer.GoodVertices_cff import *

minBiasBscFilter = cms.Sequence(physDeclFilter * bscOrBptxOr * hltMinBias * bptxAnd * bscOr * bscNoHalo)
minBiasBscFilterMC = cms.Sequence(bscOrBptxOr * bscOr * bscNoHalo)


# =========== Jet Eta and Phi selector ===================
from SpectraAna.JetEtaPhiFilter.jetetaphifilter_cfi import *
jetcut = jetfilter.clone(minEta = cms.untracked.double(-1.5),
                         maxEta = cms.untracked.double(0.0),
                         minPhi = cms.untracked.double(0.2),
                         maxPhi = cms.untracked.double(0.8)
                         )
jetEvtFilter = cms.Sequence(jetcut)

