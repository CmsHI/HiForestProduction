#include <vector>
#include <iostream>
#include <algorithm>
#include "Math/Vector3D.h"
#include <TRandom.h>

#define MAXHITS 1000
#define PI 3.141592653589

using namespace std;

class RecoHit {
   public:

   RecoHit(double _eta,double _phi,double _r) 
   { 
      eta = _eta;
      phi = _phi;
      r = _r;

      geneta=-9999;
      genphi=-9999;
      genpt=-9999;
      genid=-9999;
   }; 
   RecoHit(double _eta,double _phi,double _r,double _l) 
   { 
      eta = _eta;
      phi = _phi;
      r = _r;
      layer = _l;
      
      geneta=-9999;
      genphi=-9999;
      genpt=-9999;
      genid=-9999;
   }; 
   
   ~RecoHit(){};
   
      double eta;
      double phi;
      double r;
      double layer;
      
      double geneta;
      double genphi;
      double genpt;
      double genid;
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

  float eta1[MAXHITS],phi1[MAXHITS],r1[MAXHITS],eta2[MAXHITS],phi2[MAXHITS],r2[MAXHITS],eta[10*MAXHITS],phi[10*MAXHITS],pt[10*MAXHITS],vz[100];
  int nhits1,nhits2,mult,nv,npart;

  int gp1[MAXHITS],gp2[MAXHITS], pdg[MAXHITS], chg[MAXHITS];
  int id1[MAXHITS],id2[MAXHITS];

};

class TrackletData {
 public:

  float eta1[1000],phi1[1000],eta2[1000],phi2[1000],vz[100];
  float deta[1000],dphi[1000];
  float eta[1000],phi[1000],chg[1000],pdg[1000],nhad[12];
  int nTracklet,nhit1,nhit2,mult,nv,npart;
};

bool compareEta(RecoHit a,RecoHit b) { return a.eta<b.eta;}
bool comparePhi(RecoHit a,RecoHit b) { return a.phi<b.phi;}
bool compareAbsEta(RecoHit a,RecoHit b) { return fabs(a.eta)<fabs(b.eta);}

vector<RecoHit> removeDoubleHits(Parameters par, SelectionCriteria cuts,Int_t layer)
{
  vector<RecoHit> hits;
  vector<RecoHit> cleanedHits;

  if (layer == 1) {
    for(int ihit = 0; ihit < par.nhits1; ++ihit){
      RecoHit tmp(par.eta1[ihit],par.phi1[ihit],par.r1[ihit]);

      if(par.gp1[ihit] != -9999){
	tmp.genid = par.gp1[ihit];
	tmp.geneta = par.eta[par.gp1[ihit]];
	tmp.genphi = par.phi[par.gp1[ihit]];
	tmp.genpt = par.pt[par.gp1[ihit]];
      }

      hits.push_back(tmp);
    }
  } else {
    for(int ihit = 0; ihit < par.nhits2; ++ihit){
      RecoHit tmp(par.eta2[ihit],par.phi2[ihit],par.r2[ihit]);

      if(par.gp2[ihit] != -9999){
	tmp.geneta = par.eta[par.gp2[ihit]];
	tmp.genphi = par.phi[par.gp2[ihit]];
	tmp.genpt = par.pt[par.gp2[ihit]];
	tmp.genid = par.gp2[ihit];
      }

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
    //ROOT::Math::XYZVector tmpVector(x,y,z);

    RecoHit tmpHit(tmpVector.eta(),tmpVector.phi(),tmpVector.rho());
    tmpHit.geneta = hits[ihit].geneta;
    tmpHit.genphi = hits[ihit].genphi;
    tmpHit.genpt = hits[ihit].genpt;
    tmpHit.genid = hits[ihit].genid;

    cleanedHits.push_back(tmpHit);      
  }
  return cleanedHits;
}

RecoHit RandomHit(double etaMin, double etaMax, double phiMin, double phiMax)
{
   double eta = etaMin + (etaMax-etaMin)*gRandom->Rndm();
   double phi = phiMin + (phiMax-phiMin)*gRandom->Rndm();
   RecoHit myRandomHit(eta,phi,0);
   return myRandomHit;
}
