import FWCore.ParameterSet.Config as cms

# this is to enable or disable features as necessary for running over MC
def enableSIM(process):
    process.preTrgAna.isGEN=True
    process.postTrgAna.isGEN=True
    process.postEvtSelAna.isGEN=True
    process.postVtxAna.isGEN=True
    process.postTrkVtxAna.isGEN=True
    process.trackAna.isGEN=True
    process.looseTrackAna.isGEN=True
    process.rootpleProducer.OnlyRECO=False
    process.eventFilter.remove(process.physDeclFilter) # always false in MC
    process.eventFilter.remove(process.bptxAnd)        # always false in MC
    process.trkEffAnalyzer.hasSimInfo=True
    process.loosetrkEffAnalyzer.hasSimInfo=True
    process.extraReco *= process.trackingParticleRecoTrackAsssociation
    process.preTrgTest *= process.preTrackAna
    process.postVtxAna.vtxWeight=cms.untracked.bool(True)
    process.postTrkVtxAna.vtxWeight=cms.untracked.bool(True)
    return process

def enableSIMnoTP(process):
    process.preTrgAna.isGEN=True
    process.postTrgAna.isGEN=True
    process.postEvtSelAna.isGEN=True
    process.postVtxAna.isGEN=True
    process.postTrkVtxAna.isGEN=True
    process.trackAna.isGEN=True
    process.looseTrackAna.isGEN=True
    process.rootpleProducer.OnlyRECO=False
    process.eventFilter.remove(process.physDeclFilter) # always false in MC
    process.eventFilter.remove(process.bptxAnd)        # always false in MC
    process.trkEffAnalyzer.hasSimInfo=False
    process.loosetrkEffAnalyzer.hasSimInfo=False
    #process.extraReco *= process.trackingParticleRecoTrackAsssociation
    process.preTrgTest *= process.preTrackAna
    process.postVtxAna.vtxWeight=cms.untracked.bool(True)
    process.postTrkVtxAna.vtxWeight=cms.untracked.bool(True)
    return process

def removeTPAssociation(process):
    #process.trkEffAnalyzer.hasSimInfo=False #overwriteen above
    #process.loosetrkEffAnalyzer.hasSimInfo=False #overwritten above
    process.extraReco.remove(process.trackRefit)
    return process

def removeDijetAna(process):
    #process.analysisSeq.remove(process.djcalo)
    return process

def enable41X(process):
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.postEvtSelAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.postVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.postTrkVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.preTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.trackAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.looseTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.trackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.looseTrackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    process.refitTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','HLT')
    #process.djcalo.hltsrc=cms.InputTag('TriggerResults','','HLT')
    process.eventFilter.remove(process.hltMinBias)     #
    return process
                                                     
    
# this is for Summer 09 samples where the HLT has been re-run during Spring 10 production
def enableREDIGI(process):
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.postEvtSelAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.postVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.postTrkVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.preTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.trackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.looseTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.trackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.looseTrackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.refitTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI')
    process.djcalo.hltsrc=cms.InputTag('TriggerResults','','REDIGI')
    process.eventFilter.remove(process.hltMinBias)     # 
    return process


def enableREDIGI2(process):
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.postEvtSelAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.postVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.postTrkVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.preTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.trackAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.looseTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.trackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.looseTrackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.refitTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','ANASKIM')
    process.djcalo.hltsrc=cms.InputTag('TriggerResults','','ANASKIM')
    process.eventFilter.remove(process.hltMinBias)     #
    return process

def enableREDIGI3(process):
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.postEvtSelAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.postVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.postTrkVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.preTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.trackAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.looseTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.trackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.looseTrackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.refitTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','RECO')
    process.djcalo.hltsrc=cms.InputTag('TriggerResults','','RECO')
    process.eventFilter.remove(process.hltMinBias)     #
    return process


def enableREDIGI36X(process):
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.postEvtSelAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.postVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.postTrkVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.preTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.trackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.looseTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.trackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.looseTrackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.refitTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36X')
    process.djcalo.hltsrc=cms.InputTag('TriggerResults','','REDIGI36X')
    process.eventFilter.remove(process.hltMinBias)     #
    return process

def enableREDIGI36(process):
    process.preTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.postTrgAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.postEvtSelAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.postVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.postTrkVtxAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.preTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.trackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.looseTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.trackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.looseTrackAna_STD.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.refitTrackAna.triglabel=cms.untracked.InputTag('TriggerResults','','REDIGI36')
    process.djcalo.hltsrc=cms.InputTag('TriggerResults','','REDIGI36')
    process.eventFilter.remove(process.hltMinBias)     #
    return process
                                                    

# this is to replace the MB trigger for MC
#def gettriglist():
#    list = {

def replaceMinBiasHLTPath(process):
    list = ['HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U']
    print "hlt list for event selection analyzer = ", list
    process.evtselanalyzer.trignames=cms.untracked.vstring(list)
    process.preTrgAna.trignames=cms.untracked.vstring(list)
    process.postTrgAna.trignames=cms.untracked.vstring(list)
    process.postEvtSelAna.trignames=cms.untracked.vstring(list)
    process.postVtxAna.trignames=cms.untracked.vstring(list)
    process.postTrkVtxAna.trignames=cms.untracked.vstring(list)
    return process

def replaceMinBiasHLTPath41X(process):
    list = ['HLT_L1_BptxXOR_BscMinBiasOR','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U_v3','HLT_Jet30U_v3','HLT_Jet50U_v3']
    print "hlt list for event selection analyzer = ", list
    process.evtselanalyzer.trignames=cms.untracked.vstring(list)
    process.preTrgAna.trignames=cms.untracked.vstring(list)
    process.postTrgAna.trignames=cms.untracked.vstring(list)
    process.postEvtSelAna.trignames=cms.untracked.vstring(list)
    process.postVtxAna.trignames=cms.untracked.vstring(list)
    process.postTrkVtxAna.trignames=cms.untracked.vstring(list)
    return process


# this is to update hlt names for 900 GeV sample
def run900GeVmode(process):
    process.preTrgAna.trignames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.postTrgAna.trignames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.postEvtSelAna.trignames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.postVtxAna.trignames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.postTrkVtxAna.trignames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.preTrackAna.trignames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_MinBiasPixel_SingleTrack','HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.trackAna.hltNames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus', 'HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.trackAna.mode900GeV=cms.untracked.bool(True)
    process.looseTrackAna.hltNames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus', 'HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.looseTrackAna.mode900GeV=cms.untracked.bool(True)
    process.trackAna_STD.hltNames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus', 'HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.trackAna_STD.mode900GeV=cms.untracked.bool(True)
    process.looseTrackAna_STD.hltNames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus', 'HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.looseTrackAna_STD.mode900GeV=cms.untracked.bool(True)
    process.refitTrackAna.hltNames=cms.untracked.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus', 'HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U')
    process.refitTrackAna.mode900GeV=cms.untracked.bool(True)
    process.trkEffAnalyzer.mode900GeV=cms.bool(True)
    process.loosetrkEffAnalyzer.mode900GeV=cms.bool(True)
    return process


def run2760GeVmode(process):
    print 'bscOr and bscOrBptxOr are removed'
    process.minBiasBscFilter.remove(process.bscOr)
    process.minBiasBscFilter.remove(process.bscOrBptxOr)
    #process.bscOrBptxOr.algorithms=['L1_BscMinBiasOR_BptxPlusANDMinus']
    process.hltMinBias.HLTPaths=cms.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1')
    process.preTrgAna.trignames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_ZeroBiasPixel_SingleTrack_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    process.postTrgAna.trignames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_ZeroBiasPixel_SingleTrack_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    process.postEvtSelAna.trignames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_ZeroBiasPixel_SingleTrack_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    process.postVtxAna.trignames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_ZeroBiasPixel_SingleTrack_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    process.postTrkVtxAna.trignames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_ZeroBiasPixel_SingleTrack_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    process.preTrackAna.trignames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_ZeroBiasPixel_SingleTrack_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    process.trackAna.hltNames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    #process.trackAna.mode900GeV=cms.untracked.bool(True)
    process.looseTrackAna.hltNames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    #process.looseTrackAna.mode900GeV=cms.untracked.bool(True)
    process.trackAna_STD.hltNames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    #process.trackAna_STD.mode900GeV=cms.untracked.bool(True)
    process.looseTrackAna_STD.hltNames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    #process.looseTrackAna_STD.mode900GeV=cms.untracked.bool(True)
    process.refitTrackAna.hltNames=cms.untracked.vstring('HLT_L1BscMinBiasORBptxPlusANDMinus_v1','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1')
    #process.refitTrackAna.mode900GeV=cms.untracked.bool(True)
    #process.trkEffAnalyzer.mode900GeV=cms.bool(True)
    #process.loosetrkEffAnalyzer.mode900GeV=cms.bool(True)
    return process

def run2760GeVmodeMC(process):
    print 'bscOr and bscOrBptxOr are removed'
    process.eventFilter.remove(process.bscOr)
    process.eventFilter.remove(process.bscOrBptxOr)
    #process.minBiasBscFilter.remove(process.bscOr)
    #process.minBiasBscFilter.remove(process.bscOrBptxOr)
    return process

# this is for choosing a right global tag and regid name for an input sample
def setGlobTagAndRedigi(process,inputFileType):
    if inputFileType=='MinBias':
        print "Notice: MinBias with no REDIGI"
        process.GlobalTag.globaltag = 'START3X_V26A::All'
    if inputFileType=='MinBias_REDIGI35':
        print "Notice: MinBias with REDIGI35"
        process.GlobalTag.globaltag = 'START3X_V26A::All'
        process = enableREDIGI(process)
    if inputFileType=='Herwig_REDIGI35':
        print "Notice: Herwig with REDIGI35"
        process.GlobalTag.globaltag = 'START3X_V26A::All'
        process = enableREDIGI(process)
    if inputFileType=='MinBias900GeV':
        print "Notice: MinBias900GeV with no REDIGI"
        process.GlobalTag.globaltag = 'START36_V10::All'
    if inputFileType=='MinBias_REDIGI36':
        print "Notice: MinBias with REDIGI36"
        process.GlobalTag.globaltag = 'START36_V10::All'
        process = enableREDIGI36(process)
    if inputFileType=='QCDPtX_REDIGI36X':
        print "Notice: QCDPtX with REDIGI36X"
        process.GlobalTag.globaltag = 'START36_V9::All'
        process = enableREDIGI36X(process)
    if inputFileType=='MC_41X':
        print "Notice: MB and QCDPtX in 41X"
        process.GlobalTag.globaltag = 'START311_V2A::All'
        process = enable41X(process)
    print ("Notice: GlobalTag = " + str(process.GlobalTag.globaltag))
    return process


# this is for running of AOD skims that already have extra reco/pat objects
def enableAOD(process):
    process.extraReco.remove(process.extraVertex)
    process.analysisSeq.remove(process.patAnaSequence)
    process.extraReco.remove(process.trackRefit)
    #process.analysisSeq.remove(process.refitTrackAna)
    return process


### this is for selecting jet events instead of prescaled MB
def enableHLTJet(process,hltname='HLT_Jet15U'):
    process.minBiasBscFilter.remove(process.bscOrBptxOr)
    process.minBiasBscFilter.remove(process.hltMinBias)
    process.hltJets.HLTPaths = [hltname]
    process.minBiasBscFilter = cms.Sequence(process.hltJets*process.minBiasBscFilter)
    # check
    print "Skim jet HLT path: ", process.hltJets.HLTPaths
    return process

### this is for make a single hit HF event selection for the AGR (1charged trk) evt selection
def enableEitherHFEvtSel(process):
    process.eventFilter.replace(process.hfCoincFilter,process.hfEitherOfTwoFilter)
    return process

### this is for GEN analyzer for Pt_hat: 0 ~ 15
def enableMinPtHatCut(process):
    process.preTrackAna.pthatCut=cms.untracked.double(15.0)
    return process

def enableMinPtHatCutAuto(process,genTypePtHatRange):
    if genTypePtHatRange=='NSD_0_to_5':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(5.0)
    if genTypePtHatRange=='NSD_0_to_10':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(10.0)
    if genTypePtHatRange=='NSD_0_to_15':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(15.0)
    if genTypePtHatRange=='NSD_0_to_20':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(20.0)
    if genTypePtHatRange=='ENSD_0_to_5':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(5.0)
    if genTypePtHatRange=='ENSD_0_to_10':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(10.0)
    if genTypePtHatRange=='ENSD_0_to_15':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(15.0)
    if genTypePtHatRange=='ENSD_0_to_20':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(20.0)
    if genTypePtHatRange=='ENSD_0_to_30':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(30.0)
    if genTypePtHatRange=='MB_0_to_5':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(5.0)
    if genTypePtHatRange=='MB_0_to_10':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(10.0)
    if genTypePtHatRange=='MB_0_to_15':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(15.0)
    if genTypePtHatRange=='MB_0_to_20':
        print "Min pT_hat cut enable"
        process.preTrackAna.pthatCut=cms.untracked.double(20.0)
    return process

## to generate and analyze pythia 8 events
def enablePythia8Mode(process):
    print "Analyzer in Pythia8 mode"
    process.preTrackAna.pythia6=cms.untracked.bool(False)
    return process
    

def enable900GeVGENMode(process,genSqrts):
    if genSqrts < 910 :
        print "Sqrts is less then 910 GeV"
        process.preTrackAna.mode900GeV=cms.untracked.bool(True)
    return process                    

## this is to change vertex collection to sortedSumPtOrdVertices
def runWithsortedSumPtOrdVertices(process):
    print "sortedSumPtOrdVertices is used for final spectra!!"
    process.postTrkVtxAna.vtxlabel=cms.untracked.InputTag("sortedSumPtOrdVertices")
    process.postTrkVtxSel.vtxlabel=cms.untracked.InputTag("sortedSumPtOrdVertices")
    process.selectTracks.vertices = cms.InputTag("sortedSumPtOrdVertices")
    process.trackAna.vsrc=cms.untracked.InputTag("sortedSumPtOrdVertices")
    process.refitTrackAna.vsrc=cms.untracked.InputTag("sortedSumPtOrdVertices")
    process.trackAna_STD.vsrc=cms.untracked.InputTag("sortedSumPtOrdVertices")
    process.trkEffAnalyzer.vertices=cms.untracked.InputTag("sortedSumPtOrdVertices")
    return process

## this is to change the number of max vertices in event
def setMaxNumberVertices(process,maxnum):
    print "Max number of vertices in each event (GoodVertices, SumPtOrdVertices) set to", maxnum
    process.sortedGoodVertices.maxNumber=cms.uint32(maxnum)
    process.sortedSumPtOrdVertices.maxNumber=cms.uint32(maxnum)
    print "double-check (GoodVertices): ", process.sortedGoodVertices.maxNumber
    print "double-check (SumPtOrdVertices): ", process.sortedSumPtOrdVertices.maxNumber
    return process

### this is for taking the correct event selection efficiency
def getEvtSelEff(evtSelType):
    eff = {
      # NSD with STD
      'STD_NSD_TrkVtx':[30,0,9.0384e-05,0.00223004,0.0304791,0.217088,0.460085,0.673905,0.830166,0.92802,0.978678,1.00044,1.00497,1.00382,1.00041,0.994874,0.99174,0.990136,0.989347,0.989595,0.989071,0.990078,0.990065,0.990801,0.991051,0.991619,0.99187,0.992327,0.992703,0.99281,0.992914],
      'STD_NSD_PixVtx':[30,0.0266389,0.352064,0.55773,0.709436,0.821248,0.89875,0.956152,0.995745,1.01817,1.02657,1.02549,1.01831,1.01144,1.00505,0.998037,0.994116,0.992235,0.991112,0.991103,0.990337,0.991335,0.991312,0.991832,0.992135,0.992577,0.992819,0.993234,0.993617,0.993562,0.993773],
      # Inel with AGR Single HF EvtSel
      'AGR_Inel_EitherHF_TrkVtx':[30,0.30361,0.765493,0.915144,0.968228,0.985074,0.990709,0.992796,0.993328,0.993417,0.993582,0.993913,0.993383,0.99379,0.993978,0.993733,0.993953,0.99387,0.993817,0.993735,0.993942,0.995149,0.994863,0.993719,0.993205,0.992893,0.993857,0.993767,0.99413,0.997728,0.995064],
      'AGR_Inel_EitherHF_PixVtx':[30,0.566018,0.920571,0.963983,0.981874,0.989162,0.992304,0.993871,0.994226,0.994387,0.99452,0.994759,0.994392,0.994667,0.995005,0.99485,0.994987,0.994639,0.994662,0.994711,0.994989,0.995553,0.995464,0.994314,0.993933,0.994078,0.994233,0.99439,0.994913,0.997981,0.995393],
      # 900GeV NSD with STD
      'STD_NSD_900GeV_TrkVtx': [42,0,9.60646e-05,0.00262191,0.0337782,0.228173,0.463418,0.661808,0.798992,0.88422,0.929808,0.954356,0.964662,0.968794,0.970072,0.970321,0.969551,0.970247,0.971157,0.972393,0.97493,0.974667,0.976915,0.977959,0.979166,0.979742,0.98169,0.982051,0.982586,0.983429,0.984924,0.984952,0.98544,0.98525,0.986517,0.986371,0.986602,0.988167,0.988759,0.989259,0.98922,0.989398,0.9895],
      'STD_NSD_900GeV_PixVtx': [30,0.051323,0.459641,0.654303,0.779912,0.862804,0.918859,0.955648,0.976691,0.987668,0.992159,0.993694,0.991707,0.989005,0.987004,0.985711,0.984101,0.984,0.984631,0.985698,0.987195,0.987634,0.989591,0.990561,0.991783,0.992561,0.993581,0.993937,0.994328,0.994914,0.995531],
      # 2760GeV NSD with STD
      'STD_NSD_2760GeV_TrkVtx': [44,0,9.21022e-05,0.00129748,0.0235032,0.179341,0.40027,0.614643,0.76444,0.864228,0.923655,0.951388,0.96699,0.973472,0.973865,0.969131,0.969939,0.969537,0.968224,0.967431,0.970173,0.968282,0.972264,0.974348,0.973122,0.972457,0.974793,0.979822,0.980527,0.98035,0.979858,0.97984,0.983024,0.984624,0.984346,0.986562,0.986131,0.989628,0.989913,0.98857,0.990146,0.991374,0.990627,0.989836,0.991397],
      'STD_NSD_2760GeV_PixVtx' : [36,0.0346931,0.494957,0.686681,0.82791,0.907836,0.956462,0.993745,1.0006,1.00152,1.00214,0.993881,0.990084,0.988229,0.981616,0.975313,0.973781,0.973778,0.971727,0.970195,0.973937,0.972184,0.975222,0.977703,0.976333,0.975967,0.978696,0.982506,0.984059,0.983734,0.983647,0.982092,0.985906,0.987879,0.987245,0.989494,0.98893]
      }
    return eff[evtSelType]

def updateEvtSelEff(ana,evtSelType):
  # Set the eff sel mode
  if evtSelType.find('STD') >= 0:
    ana.ptMin = 0.0
    ana.evtEffCorrType = 0
  if evtSelType.find('AGR') >= 0:
    ana.ptMin = 0.5
    ana.evtEffCorrType = 2
  # Set the eff sel eff numbers
  ana.evtSelEffv = getEvtSelEff(evtSelType)

def enableDJetAna(process,mode="MC",output="FF"):
  from Saved.DiJetAna.customise_cfi import enablePp
  from Saved.DiJetAna.customise_cfi import djOutputThreshold
  enablePp(process)
  for m in [process.djcalo,process.djcalo_tower,process.djcalo_genp,process.djgen]:
    m.vtxsrc = "sortedGoodVertices"
    m.verbosity = 1
    m.nearJetPtMin = djOutputThreshold[output]
    m.refjetsrc = "cleanedPatJets"
    if mode=="MC":
      m.hltNames[-1] = "HLT_Photon20_L1R"
  for m in [process.djcalo,process.djcalo_tower,process.djcalo_genp]:
    m.jetsrc = "cleanedPatJets"
    m.anaJetType = 2
  process.djcalo.trksrc = "selectTracks"
  process.djcalo.anaTrkType = 2
  if mode=="Data":
    from Saved.DiJetAna.customise_cfi import enableData
    enableData(process)

def enableOpenHlt(process,seq):
	process.load("HLTrigger.HLTanalyzers.HLTAnalyser_cfi")
	process.hltanalysis.RunParameters.Debug = cms.bool(False)
	process.hltanalysis.RunParameters.UseTFileService = cms.untracked.bool(True)
	process.hltanalysis.RunParameters.Monte = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoMC = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoAlCa = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoTracks = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoVertex = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoJets = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoMuons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoL1Muons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoL2Muons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoL3Muons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoOfflineMuons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoQuarkonias = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoPhotons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoSuperClusters = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoElectrons = cms.untracked.bool(False)
	process.hltanalysis.RunParameters.DoBJets = cms.untracked.bool(False)
	process.hltanalysis.l1GtReadoutRecord = "gtDigis"
	# add to seq
	seq*=process.hltanalysis
	print "Added openhlt"
