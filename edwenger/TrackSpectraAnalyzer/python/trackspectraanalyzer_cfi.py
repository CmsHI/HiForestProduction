import FWCore.ParameterSet.Config as cms

trackAna = cms.EDAnalyzer('TrackSpectraAnalyzer',
                          src = cms.untracked.InputTag("generalTracks"),
                          vsrc = cms.untracked.InputTag("offlinePrimaryVertices"),
                          jsrc = cms.untracked.InputTag("selectedPatJets"),
                          gsrc = cms.untracked.InputTag("genParticles"),
                          gjsrc = cms.untracked.InputTag("ak5GenJets"),
                          src_evtCorr = cms.untracked.InputTag("generalTracks"),
                          setQualityBit = cms.untracked.bool(True),
                          isGEN = cms.untracked.bool(False),
                          pureGENmode = cms.untracked.bool(False),
                          nsdOnly = cms.untracked.bool(False), # look at only NSD GEN
                          pythia6 = cms.untracked.bool(True), # false => pythia8
                          mode900GeV = cms.untracked.bool(False),
                          pthatCut = cms.untracked.double(0.0),
                          doJet = cms.untracked.bool(True),
                          histOnly = cms.untracked.bool(True),
                          includeExtra = cms.untracked.bool(False),
                          etaMax = cms.untracked.double(2.5),
                          ptMin = cms.untracked.double(0.0),
                          applyEvtEffCorr = cms.untracked.bool(False),
                          evtEffCorrType = cms.untracked.int32(0), # 0 = STD, 1 = SPEC, 2 = AGR
                          efit_type = cms.untracked.int32(0), # 0 = 1 pixel,  1 = 4 tracks
                          evtMultCut = cms.untracked.int32(3), # skip events with mult <= evtMultCut
                          evtSelEffv = cms.untracked.vdouble(19,0.0119454,0.212963,0.376093,0.450549,0.535354,0.607345,0.684058,0.723837,0.805643,0.83391,
                                                             0.880702,0.93007,0.959839,0.972656,0.966942,0.965812,0.982609,0.985646,0.978495),
                          hltNames = cms.untracked.vstring('HLT_L1Tech_BSC_minBias_OR', # the 1st 5 paths are hardcoded in nt, please add more HLT paths after these
                                                           'HLT_L1Jet6U',
                                                           'HLT_Jet15U',
                                                           'HLT_Jet30U',
                                                           'HLT_Jet50U'),
                          neededTrigSpectra = cms.untracked.vint32(0,0,1,0,1), # 1st, 3rd, in the hltNames list 
                          triglabel = cms.untracked.InputTag('TriggerResults','','HLT'),
                          trkAcceptedJet = cms.untracked.bool(True), # # limit jets to |eta|<1.5
                          )

