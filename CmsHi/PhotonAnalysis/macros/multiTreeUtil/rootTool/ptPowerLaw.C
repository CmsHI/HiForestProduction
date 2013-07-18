#include <TF1.h>
#include <TMath.h>
Double_t ptPowerLaw(Double_t *x, Double_t *par) {

   double dndy = par[0];
   double n    = par[1];
   double T    = par[2];
   
   double pi = 3.14159265358979;
   double m = 0.13957018;  // pion mass
   
   double a = (n-1)*(n-2)/n/T/(n*T+(n-2)*m);
   double b = x[0]*pow((1+x[0]/n/T),-n);
   
   return a*b;

}

