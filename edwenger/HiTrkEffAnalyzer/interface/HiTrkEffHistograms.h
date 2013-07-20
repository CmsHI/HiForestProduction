#ifndef HiTrkEffAnalyzer_HiTrkEffHistograms
#define HiTrkEffAnalyzer_HiTrkEffHistograms

#include <vector>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

namespace edm { class ParameterSet ; }
class TFile;
class TTree;
class TH1F;
class TH2F;
class TH3F;

// define matched track structs

typedef struct
{
  Int_t ids;
  Float_t etas;
  Float_t pts;
  Int_t hits;
  Int_t status;
  Int_t acc;
  Int_t nrec;
  Float_t ptr;
  Float_t dz;
  Float_t d0;
  Float_t pterr;
  Float_t d0err;
  Float_t dzerr;
  Int_t hitr;
  Int_t algo;
  Float_t jetr;
  Int_t cbin;
  Float_t jetar;
  Float_t jrdr;
  Int_t jrind;
  Int_t jrflavor;
} SimTrack_t;

typedef struct
{
  Int_t charge;
  Float_t etar;
  Float_t ptr;
  Float_t phir;
  Float_t dz;
  Float_t d0;
  Float_t pterr;
  Float_t d0err;
  Float_t dzerr;
  Int_t hitr; 
  Int_t algo;
  Int_t nsim;
  Int_t status;
  Int_t ids;
  Int_t parids;
  Float_t etas;
  Float_t pts;
  Float_t jetr;
  Int_t cbin;
  Float_t jetar;
  Float_t jrdr;
  Int_t jrind;
  Int_t jrflavor;
} RecTrack_t;

// define class to hold histograms and track trees

class HiTrkEffHistograms
{

 public:
   HiTrkEffHistograms(const edm::ParameterSet& pset);
   ~HiTrkEffHistograms();

  void declareHistograms();
  void fillSimHistograms(const SimTrack_t & s);
  void fillRecHistograms(const RecTrack_t & r);
  void writeHistograms();

 private: 
   std::vector<TTree*> trackTrees;
   SimTrack_t simTrackValues;
   RecTrack_t recTrackValues;
   bool fillHistograms;
   bool fillNtuples;
   Int_t ptBinScheme;
   bool lowPtMode;
   Double_t trkPtMin;
   edm::Service<TFileService> f;

   std::vector<Double_t> etaBins, ptBins, jetBins;

   std::vector<int32_t> neededCentBins;


   // SimTrack
   TH2F* hsim;
   TH2F* hacc;
   TH2F* heff;
   TH2F* hmul;

   TH3F* hsim3D;
   TH3F* heff3D;
   TH3F* hmul3D;
   TH3F* hresStoR3D;
   TH3F* hresStoR3D_etaS;
   TH3F* hresStoR3D_etaL;
   TH3F* hresStoR3D_v2;

   TH1F* hsimSPt1D; // sim dn/dpt with sim pt for mom. res. closure
   TH1F* hsimRPt1D; // sim dn/dpt with rec pt for mom. res. closure

   // RecTrack
   TH2F* hrec;
   TH2F* hfak;
   TH2F* hsec;

   TH3F* hrec3D;
   TH3F* hfak3D;
   TH3F* hsec3D;
   
   TH2F* hetaphi;


   // vector of histograms
   std::vector<TH3F*> vhsim3D;
   std::vector<TH3F*> vheff3D;
   std::vector<TH3F*> vhmul3D;

   std::vector<TH3F*> vhrec3D;
   std::vector<TH3F*> vhfak3D;
   std::vector<TH3F*> vhsec3D;
   
   std::vector<TH3F*> vhresStoR3D;





};

#endif
