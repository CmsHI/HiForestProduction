// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// for jet,track,vertices correlation 
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"


// ROOT includes
#include <TH2.h>
#include <TMath.h>

class VertexAnalyzer : public edm::EDAnalyzer {
   public:
      explicit VertexAnalyzer(const edm::ParameterSet&);
      ~VertexAnalyzer(){};

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      typedef math::XYZPoint Point;

      // ----------member data ---------------------------

      edm::InputTag vtxlabel_;
      edm::InputTag trklabel_;
      edm::InputTag jetlabel_;

      double dzcut_;
      double minJetEt_;
      unsigned int nthleadingMax_;
      int nTracksBins_;
      bool jetTrkVerticesCorr_;

      edm::Service<TFileService> f;

      TH1D *hVtxSize;
      TH1D *hVtxTracks;
      TH1D *hVtxZ;
      TH2D *hMultVtxTracks;
      TH2D *hMultVtxZ;
      TH1D *hAllVtxZ;
      
      TH1F *hLeadingTrkPVdZ;
      TH1F *hSLeadingTrkPVdZ;
      TH1F *hSSLeadingTrkPVdZ;
      TH1F *hLeadingTrkSVdZ;
      TH1F *hSLeadingTrkSVdZ;
      TH1F *hSSLeadingTrkSVdZ;

      TH1F *hLeadingTrkPVdZ_tightest;
      TH1F *hLeadingTrkSVdZ_tightest;


      TH1F *hMostCorrTrkPVdZ;
      TH1F *hMostCorrTrkSVdZ;
      
      TH1F *hMostCorrTrkPVdZ_Acc;
      TH1F *hMostCorrTrkSVdZ_Acc;

      TH1F *hMostCorrTrkPVdZ_nonAcc;
      TH1F *hMostCorrTrkSVdZ_nonAcc;

      TH1F *hTrkIndex;
      TH1F *dRofMostCorrTrk;
      TH1F *dRofMostCorrTrk_Acc;
      TH1F *dRofMostCorrTrk_nonAcc;
      TH2F *dEtadPhi_leading;

      TH2F *hLeadingAndSLeadingTrkPVdZ;
      TH2F *hLeadingAndSSLeadingTrkPVdZ;
      TH2F *hSLeadingAndSSLeadingTrkPVdZ;

      TH2F *hLeadingAndSLeadingTrkSVdZ;
      TH2F *hLeadingAndSSLeadingTrkSVdZ;
      TH2F *hSLeadingAndSSLeadingTrkSVdZ;

      TH1F *hLeadingTrkPVdZ_narrowEta; // for jet |eta|<2.0! 
      TH1F *hLeadingTrkSVdZ_narrowEta;

      TH1F *hLeadingJetEtaFromPV;
      TH1F *hLeadingJetEtaFromSV;
      TH1F *hSLeadingJetEtaFromPV;
      TH1F *hSLeadingJetEtaFromSV;
      TH1F *hSSLeadingJetEtaFromPV;
      TH1F *hSSLeadingJetEtaFromSV;

      TH1F *hLeadingJetEtFromPV;
      TH1F *hLeadingJetEtFromSV;
      TH1F *hSLeadingJetEtFromPV;
      TH1F *hSLeadingJetEtFromSV;
      TH1F *hSSLeadingJetEtFromPV;
      TH1F *hSSLeadingJetEtFromSV;

      TH1F *hLeadingdRFromPV;
      TH1F *hLeadingdRFromSV;
      TH1F *hSLeadingdRFromPV;
      TH1F *hSLeadingdRFromSV;
      TH1F *hSSLeadingdRFromPV;
      TH1F *hSSLeadingdRFromSV;

      // track-jet
      TH1F *hLeadingdR_eitherJets;
      TH1F *hSLeadingdR_eitherJets;
      TH1F *hSSLeadingdR_eitherJets;

      TH1F *hLeadingdR_leadingJets;
      TH1F *hSLeadingdR_leadingJets;
      TH1F *hSSLeadingdR_leadingJets;

      TH1F *hLeadingdR_sleadingJets;
      TH1F *hSLeadingdR_sleadingJets;
      TH1F *hSSLeadingdR_sleadingJets;

      TH2F *hLeadingPVdZdR;
      TH2F *hSLeadingPVdZdR;
      TH2F *hSSLeadingPVdZdR;
      TH2F *hLeadingSVdZdR;
      TH2F *hSLeadingSVdZdR;
      TH2F *hSSLeadingSVdZdR;

      TH1F *hTrkPtFromAV;
      TH1F *hTrkPtFromPV;
      TH1F *hTrkPtFromSV;

      TH2F *hVtxZTrkEta;
      TH1F *hLeadMonoJetEta;
      TH2F *hLeadnSLeadJetEta;

      TH2F *hJetEtTrkPtSumPV;
      TH2F *hJetEtTrkPtSumSV;

      std::vector<double> ptBins, jetBins, cptBins;

};
