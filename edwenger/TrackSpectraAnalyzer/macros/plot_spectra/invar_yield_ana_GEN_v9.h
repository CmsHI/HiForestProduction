//------------------------------------------------------------------- 
// June 03 2010, Andre S. Yoon 
// - header file for invar_yield_ana_v9.C
// - efficiency correction is done for tracking
// - efficiency correction is done for event selection 
// - return various histograms/tgraph in struct
//--------------------------------------------------------------------


//using namespace std;
//#include <utility>


TGraphErrors* TgraphIt(TH1D* hist){

   TGraphErrors *tg;
   int nbins = hist->GetNbinsX();

   const int nlines = nbins;

   float pt[nlines], xsec[nlines];
   float pterr[nlines], xsecerr[nlines];

   for(int i = 0; i<nbins; i++ ){
      pt[i] = hist->GetBinCenter(i+1);
      xsec[i] = hist->GetBinContent(i+1);
      xsecerr[i] = hist->GetBinError(i+1);
      pterr[i] = 0;
   }
   
   tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
   return tg;
}

TH1D* RebinIt(TH1D* hist, bool REBIN){

   // ########################################### prepare variable bin ######                                                                                             
   Double_t dBin = hist->GetBinWidth(1);//assume hist has constan bin width                                                                                               

   Int_t NumBins = hist->GetNbinsX();
   const Int_t Nlines = NumBins;

   Int_t binTemp = 0;
   Int_t nbins = 0;
   Double_t binsTemp[Nlines+1];
   Int_t totBins = Nlines;

   for(Int_t i = 0; i < totBins ; i++){
      binsTemp[i] = binTemp;
      binTemp++;
   }
   binsTemp[totBins] = totBins;

   Int_t bin = 0;
   Double_t bins[Nlines+1];
   Int_t nbins = 0;
   Double_t binWidth = hist->GetBinWidth(1);
   cout<<"bin width of original histogram "<<binWidth<<endl;

   if(REBIN){
      while (bin < totBins) {
         bins[nbins] = binWidth*binsTemp[bin];
         //cout<<"bins[nbins] = "<<bins[nbins]<<endl;                                                                                                                     
         nbins++;
	 // MC
	 /*
         if (bin < 2)          bin += 1;
         else if(bin < 8)      bin += 2;
         else if (bin < 12)    bin += 4;
         else if (bin < 30)    bin += 6;
         else if (bin < 40)    bin += 8;
         else if (bin < 70)    bin += 10;
         else if (bin < 100)   bin += 15;
         else if (bin < 200)   bin += 20;
         else                  bin += 30;
	 */
	 // DATA

	 if (bin < 3)          bin += 1;
	 else if(bin < 8)      bin += 2;
	 else if (bin < 12)    bin += 4;
	 else if (bin < 30)    bin += 6;
	 else if (bin < 40)    bin += 8;
	 else if (bin < 70)    bin += 25;
	 else if (bin < 100)   bin += 40;
	 else if (bin < 200)   bin += 55;
	 else                  bin += 100;

      }
      bins[nbins] = binWidth*binTemp[totBins];
      //cout<<"bins[nbins] = "<<bins[nbins]<<endl;                                                                                                                        
   }

   cout<<"number of bins after rebinned = "<<nbins<<endl;

   //################################################# Rebin ###################                                                                                          

   Char_t hname[500];
   sprintf(hname,"pre_hRInvX");
   hist->Rebin(nbins,hname,bins);

   pre_hRInvX->SetNameTitle(hname,hname);
   cout<<"[RebinIt][check] name of the rebinned histogram  "<<pre_hRInvX->GetName()<<endl;

   const Int_t RNlines = nbins;

   Float_t ptR[RNlines], xsecR[RNlines];
   Float_t err_ptR[RNlines], err_xsecR[RNlines];


   for(Int_t j = 0; j<nbins; j++ ){

      // In order to scale the content properly (due to rebinning)                                                                                                       \
                                                                                                                                                                          
      Float_t dbin = pre_hRInvX->GetBinWidth(j+1);
      Float_t ratio = dbin/dBin;
      Float_t content = pre_hRInvX->GetBinContent(j+1);
      Float_t err = pre_hRInvX->GetBinError(j+1);

      pre_hRInvX->SetBinContent(j+1,(content/ratio));
      pre_hRInvX->SetBinError(j+1,(err/ratio));

      ptR[j] = pre_hRInvX->GetBinCenter(j+1);
      xsecR[j] = (content/ratio);
      err_ptR[j] = 0;
      err_xsecR[j] = (err/ratio);

   }

   //data.hRInvX = (TH1D*) pre_hRInvX->Clone("hRInvX");                                                                                                                   
   return pre_hRInvX;

}


