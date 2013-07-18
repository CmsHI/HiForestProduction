#include <vector>
#include "Math/Vector3D.h"

using namespace std;

class RecoHit {
   public:

   RecoHit(double _eta,double _phi,double _r) 
   { 
      eta = _eta;
      phi = _phi;
      r = _r;
   }; 
   ~RecoHit(){};
   
      double eta;
      double phi;
      double r;
};

class SelectionCriteria {
 public:

  double drCut   ;       // to remove double hit
  double dPhiCut ;       // to remove double hit
  double vzCut   ;       // vertex cut

  bool verbose_ ;
  bool useDeltaPhi_;
  bool checkSecondLayer_;
};

class Parameters {
 public:

  float eta1[1000],phi1[1000],r1[1000],eta2[1000],phi2[1000],r2[1000],vz[100];
  int nhits1,nhits2,mult,nv;
};

bool compareEta(RecoHit a,RecoHit b) { return a.eta<b.eta;}
bool comparePhi(RecoHit a,RecoHit b) { return a.phi<b.phi;}

vector<RecoHit> removeDoubleHits(Parameters par, SelectionCriteria cuts,Int_t layer)
{
  vector<RecoHit> hits;
  vector<RecoHit> cleanedHits;

  if (layer == 1) {
    for(int ihit = 0; ihit < par.nhits1; ++ihit){
      RecoHit tmp(par.eta1[ihit],par.phi1[ihit],par.r1[ihit]);
      hits.push_back(tmp);
    }
  } else {
    for(int ihit = 0; ihit < par.nhits2; ++ihit){
      RecoHit tmp(par.eta2[ihit],par.phi2[ihit],par.r2[ihit]);
      hits.push_back(tmp);
    }
  }
  sort (hits.begin(),hits.end(),comparePhi);
    
  for(int ihit = 0; ihit < (int)hits.size(); ++ihit) {
    double dr=0;
    double dphi=10;
    if (ihit !=0) {
      dphi = fabs(hits[ihit-1].phi - hits[ihit].phi);
      dr   = fabs(hits[ihit-1].r - hits[ihit].r);
    }
      
    if (dr>cuts.drCut && dphi < cuts.dPhiCut) continue;
      
    // recalculate eta and phi
    double x = hits[ihit].r*cos(hits[ihit].phi);
    double y = hits[ihit].r*sin(hits[ihit].phi);
    double z = hits[ihit].r/tan(atan(exp(-hits[ihit].eta))*2);
    ROOT::Math::XYZVector tmpVector(x,y,z-par.vz[1]);
    RecoHit tmpHit(tmpVector.eta(),tmpVector.phi(),tmpVector.rho());
    cleanedHits.push_back(tmpHit);      
  }
  return cleanedHits;
}


