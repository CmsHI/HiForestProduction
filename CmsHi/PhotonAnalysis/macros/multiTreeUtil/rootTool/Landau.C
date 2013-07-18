#include <TF1.h>
#include <TMath.h>
Double_t LandauFunction(Double_t *x, Double_t *par) {
  Double_t xx = x[0];
  return par[0]*TMath::Landau(xx,par[1],par[2]);
}

