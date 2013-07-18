#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>

// x: observable, y: truth

class histFunction2D
{ 
   public:
      histFunction2D(TH2F *h);
      ~histFunction2D();

      int GetBinNumber(Double_t c, Double_t lowEdge, Double_t highEdge, Double_t nbin);
      Double_t evaluate(Double_t *x, Double_t *par);
      TH2F *hist;
      Double_t lowEdgeX;
      Double_t highEdgeX;
      Double_t lowEdgeY;
      Double_t highEdgeY;
      Double_t nbinX;
      Double_t nbinY;
};

histFunction2D::histFunction2D(TH2F *h)
{
   hist  = (TH2F*)h->Clone();
   nbinX = hist->GetNbinsX();
   nbinY = hist->GetNbinsY();
   lowEdgeX  = hist->GetXaxis()->GetBinLowEdge(1);
   highEdgeX = hist->GetXaxis()->GetBinLowEdge(nbinX+1);
   lowEdgeY  = hist->GetYaxis()->GetBinLowEdge(1);
   highEdgeY = hist->GetYaxis()->GetBinLowEdge(nbinY+1);
   
   // normalize
   for (int i=1;i<=nbinX;i++) {
      Double_t total=0;
      for (int j=1;j<=nbinY;j++) {
         total+=hist->GetBinContent(i,j);
      }
      for (int j=1;j<=nbinY;j++) {
         if (total!=0) hist->SetBinContent(i,j,hist->GetBinContent(i,j)/total);
      }
   }
}

int histFunction2D::GetBinNumber(Double_t c, Double_t lowEdge, Double_t highEdge, Double_t nbin)
{
   Double_t delta = c - lowEdge;
   Double_t n = delta/(highEdge-lowEdge)*nbin;
   return (int)(n+1);
}

Double_t histFunction2D::evaluate(Double_t *x, Double_t *par) {
     
  Double_t xx = x[0];
  int binY=GetBinNumber(xx,lowEdgeY,highEdgeY,nbinY);
  
  Double_t total=0;
  for (int j=0;j<nbinX;j++) {  
     total+=fabs(par[j])*hist->GetBinContent(j+1,binY); 
  }
  
  // add a constant term 
  total+=par[(int)nbinX];
  return total;
}

/* example: 

   histFunction2D *myfun = new histFunction2D(h);
   
   TF1 *test = new TF1("histFun",myfun,&histFunction2D::evaluate,0,maxDphi,nPtBin+1);
*/
