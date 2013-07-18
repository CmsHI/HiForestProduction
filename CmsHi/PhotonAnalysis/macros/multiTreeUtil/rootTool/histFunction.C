#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>

class histFunction
{ 
   public:
      histFunction(TH1F *h);
      ~histFunction(){ delete hist;};

      int GetBinNumber(double c);
      Double_t evaluate(Double_t *x, Double_t *par);
      TH1F *hist;
      Double_t lowEdge;
      Double_t highEdge;
      Double_t nbin;
};

histFunction::histFunction(TH1F *h)
{
   hist=(TH1F*)h->Clone();
   nbin=h->GetNbinsX();
   lowEdge=h->GetBinLowEdge(1);
   highEdge=h->GetBinLowEdge(nbin+1);
   hist->SetName("hHistFunction");
   hist->Scale(1./hist->GetEntries());
}

int histFunction::GetBinNumber(double c)
{

   double delta = c - lowEdge;
   double n = delta/(highEdge-lowEdge)*nbin;
   return (int)(n+1);
}

Double_t histFunction::evaluate(Double_t *x, Double_t *par) {
     
  Double_t xx = x[0];
  int nbin=hist->FindBin(xx);  //
  nbin=GetBinNumber(xx);
  
  return par[0]*hist->GetBinContent(nbin);
}

/* example: 

   histFunction *myfun = new histFunction(hist);
   TF1 *test = new TF1("landau",myfun,&histFunction::evaluate,0,0.7,3);

*/
