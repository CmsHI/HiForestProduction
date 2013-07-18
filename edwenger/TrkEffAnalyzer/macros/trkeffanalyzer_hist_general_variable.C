////

using namespace std;

struct trkeffanalyzer_hist_data {
   TGraphAsymmErrors *gEfficiency;
   TGraphAsymmErrors *gFakerate;
   TGraph2DErrors  *gEff2D;
   TGraph2DErrors  *gFR2D;
   TH2D  *sim_eta_pt;
   //TH2D  *gEff2D;
   //TH2D  *gFR2D;


};

trkeffanalyzer_hist_data trkeffanalyzer_hist_graph(Char_t *cFile, Bool_t REBIN = true, Bool_t RANGE1 = true, Int_t variable, 
						       Double_t ieta, Double_t feta, Double_t imult, Double_t fmult, Double_t minpt, Double_t maxpt)
{

   cout<<"Beginning.."<<endl;
   cout<<"File to be processed is "<<cFile<<endl;
   TFile *tracking = new TFile(cFile,"read"); 

   TDirectoryFile *dSub = tracking->GetDirectory("trkEffAnalyzer");


   dSub->cd();
   TIter next(dSub->GetListOfKeys());                 
   TKey *key;

   Char_t name[100];
   Char_t hist1[100], hist2[100], hist3[100], hist4[100];

   TH1D *hDEff,*hDenEff;
   TH1D *hNEff,*hNumEff;
   TH1D *hDFake,*hDenFake;
   TH1D *hNFake,*hNumFake,*hNumFakeS;

   sprintf(hist1,"heff3D");
   sprintf(hist2,"hsim3D");
   sprintf(hist3,"hfak3D");
   sprintf(hist4,"hrec3D");

   while ((key=(TKey*)next())) {
      sprintf(name,"%s",key->GetName());
      if(strcmp((key->GetName()),(hist1))==0){
         cout<<"Your numerator for Eff "<<name<<endl;
         TH3F *hPNEff3D = (TH3F*) dSub->Get(name);
      }
      if(strcmp((key->GetName()),(hist2))==0){
         cout<<"Your denominator for Eff "<<name<<endl;
         TH3F *hPDEff3D = (TH3F*) dSub->Get(name);
      }
      if(strcmp((key->GetName()),(hist3))==0){
         cout<<"Your numerator for FR "<<name<<endl;
         TH3F *hPNFake3D = (TH3F*) dSub->Get(name);
      }
      if(strcmp((key->GetName()),(hist4))==0){
         cout<<"Your denominator for FR "<<name<<endl;
         TH3F *hPDFake3D = (TH3F*) dSub->Get(name);
      }
   }


   
   // x-eta, y-pt, z-mult 

   // --------------  variable for eta projection
   Int_t numBinsX = hPNEff3D->GetNbinsX();
   Double_t etaMin = hPNEff3D->GetXaxis()->GetBinLowEdge(1);
   Double_t etaMax = hPNEff3D->GetXaxis()->GetBinUpEdge(numBinsX);

   cout<<"Number of bins in eta (x-axis) :"<<numBinsX<<endl;
   cout<<"[in this histogram] Min eta: "<<etaMin<<" && Max eta: "<<etaMax<<endl;
   cout<<"[range for eff/fake] |Min eta|: "<<ieta<<" && |Max eta|: "<<feta<<endl;

   Int_t binMinP = hPNEff3D->GetXaxis()->FindBin(ieta);
   Int_t binMaxP = hPNEff3D->GetXaxis()->FindBin(feta);

   Int_t binMinN = hPNEff3D->GetXaxis()->FindBin(-1.0*ieta);
   Int_t binMaxN = hPNEff3D->GetXaxis()->FindBin(-1.0*feta);

   // --------------  variable for pt projection 
   Int_t numBinsY = hPNEff3D->GetNbinsY();
   Double_t ptMin = hPNEff3D->GetYaxis()->GetBinLowEdge(1);
   Double_t ptMax = hPNEff3D->GetYaxis()->GetBinUpEdge(numBinsY);

   cout<<"Number of bins in pt (y-axis) :"<<numBinsY<<endl;
   cout<<"[in this histogram] Min pt: "<<ptMin<<" && Max pt: "<<ptMax<<endl;
   cout<<"[range for eff/fake] Min pt: "<<minpt<<" && Min pt: "<<maxpt<<endl;

   Int_t binYMin = hPNEff3D->GetYaxis()->FindBin(minpt);
   Int_t binYMax = hPNEff3D->GetYaxis()->FindBin(maxpt);

   // --------------  variable for jet projection
   Int_t numBinsZ = hPNEff3D->GetNbinsZ();
   Double_t jetMin = hPNEff3D->GetZaxis()->GetBinLowEdge(1);
   Double_t jetMax = hPNEff3D->GetZaxis()->GetBinUpEdge(numBinsZ);

   cout<<"Number of bins in jet (z-axis) :"<<numBinsZ<<endl;
   cout<<"[in this histogram] Min jet: "<<jetMin<<" && Max jet: "<<jetMax<<endl;
   cout<<"[range for eff/fake] Min jet: "<<imult<<" && Min jet: "<<fmult<<endl;

   Int_t binZMin = hPNEff3D->GetZaxis()->FindBin(imult);
   Int_t binZMax = hPNEff3D->GetZaxis()->FindBin(fmult);

   if(variable==1){
      hNEff_pos = (TH1D*) hPNEff3D->ProjectionX("hNEff_pos",binYMin,binYMax,binZMin,binZMax,"e");
      hDEff_pos = (TH1D*) hPDEff3D->ProjectionX("hDEff_pos",binYMin,binYMax,binZMin,binZMax,"e");
      hNFake_pos = (TH1D*) hPNFake3D->ProjectionX("hNFake_pos",binYMin,binYMax,binZMin,binZMax,"e");
      hDFake_pos = (TH1D*) hPDFake3D->ProjectionX("hDFake_pos",binYMin,binYMax,binZMin,binZMax,"e");
   }else if(variable==2){
      hNEff_pos = (TH1D*) hPNEff3D->ProjectionY("hNEff_pos",binMaxN,binMaxP,binZMin,binZMax,"e");
      hDEff_pos = (TH1D*) hPDEff3D->ProjectionY("hDEff_pos",binMaxN,binMaxP,binZMin,binZMax,"e");
      hNFake_pos = (TH1D*) hPNFake3D->ProjectionY("hNFake_pos",binMaxN,binMaxP,binZMin,binZMax,"e");
      hDFake_pos = (TH1D*) hPDFake3D->ProjectionY("hDFake_pos",binMaxN,binMaxP,binZMin,binZMax,"e");
   }else if(variable==3){
      hNEff_pos = (TH1D*) hPNEff3D->ProjectionZ("hNEff_pos",binMaxN,binMaxP,binYMin,binYMax,"e");
      hDEff_pos = (TH1D*) hPDEff3D->ProjectionZ("hDEff_pos",binMaxN,binMaxP,binYMin,binYMax,"e");
      hNFake_pos = (TH1D*) hPNFake3D->ProjectionZ("hNFake_pos",binMaxN,binMaxP,binYMin,binYMax,"e");
      hDFake_pos = (TH1D*) hPDFake3D->ProjectionZ("hDFake_pos",binMaxN,binMaxP,binYMin,binYMax,"e");
   }

   hNEff   = (TH1D*) hNEff_pos->Clone("hNEff");
   hDEff   = (TH1D*) hDEff_pos->Clone("hDEff");
   hNFake  = (TH1D*) hNFake_pos->Clone("hNFake");
   hDFake  = (TH1D*) hDFake_pos->Clone("hDFake");

   
   //Rebinning machine..
   Int_t bin = 1;
   Double_t bins[5000];
   
   //Float_t fPtBins[]={3.,4.0,10.,25.,55.,200.};                                                  
   //Int_t iPtBins=6;                                                                                               
   Int_t nbins = 0;
   Int_t totBins = hNEff->GetNbinsX();
   cout<<"Number of bins "<<totBins<<endl;
   
   TAxis *xaxis = hNEff->GetXaxis();
   
   while (bin < totBins) {
     bins[nbins] = xaxis->GetBinLowEdge(bin);
     nbins++;
     if(REBIN){
       if(RANGE1){

	  //bin += 4;
	  /*
	  if      (bin < 4)     bin += 2;
	  else if (bin < 8)    bin += 3;
          else if (bin < 16)    bin += 4;
          else if (bin < 26)    bin += 6;
          else if (bin < 35)    bin += 18;
          //else if (bin < 55)    bin += 40;
          else                  bin += 80;
	  */


	  /*
	  if (bin < 5)          bin += 2;
	  else if (bin < 11)    bin += 3;
	  else if (bin < 17)    bin += 5;
	  else if (bin < 60)    bin += 7;
	  else if (bin < 100)   bin += 12;
	  else                  bin += 18;
	  */


	  if (bin < 2)          bin += 1;
	  else if(bin < 8)      bin += 2;
	  else if (bin < 12)    bin += 4;
	  else if (bin < 30)    bin += 6;
	  else if (bin < 40)    bin += 8;
	  else if (bin < 70)    bin += 25;
	  else if (bin < 100)   bin += 40;
	  else if (bin < 200)   bin += 50;
	  else                  bin += 65;

	  //bin += 8;

	  /*
	  if (bin < 2)          bin += 1;
	  else if(bin < 8)      bin += 2;
	  else if (bin < 12)    bin += 4;
	  else if (bin < 30)    bin += 6;
	  else if (bin < 40)    bin += 8;
	  else if (bin < 70)    bin += 25;
	  else if (bin < 100)   bin += 40;
	  else if (bin < 200)   bin += 55;
	  else                  bin += 100;
	  */

	  /*
	  if      (bin < 5)     bin += 2;
          else if (bin < 10)    bin += 4;
          else if (bin < 20)    bin += 6;
          else if (bin < 25)    bin += 8;
          else if (bin < 40)    bin += 12;
          else                  bin += 20;
	  */

	  /*
	  if      (bin < 5)     bin += 2;
          else if (bin < 10)    bin += 3;
          else if (bin < 20)    bin += 4;
	  else if (bin < 25)    bin += 10;
          else if (bin < 30)    bin += 20;
          else if (bin < 40)    bin += 20;
          else if (bin < 100 )  bin += 20;
          else if (bin < 250 )  bin += 50;
          else if (bin < 550 )  bin += 60;
          else                  bin += 70;
	  */
       }else{
	  /*
	 if (bin < 30)         bin += 30;
         else if (bin < 40)    bin += 30;
         else if (bin < 100 )  bin += 60;
         else if (bin < 250 )  bin += 150;//150
         else if (bin < 550 )  bin += 300;//300
         else                  bin += 400;//400
	  */
	  /*
	  if (bin < 30)         bin += 10;
	  else if (bin < 40)    bin += 15;
	  else if (bin < 100 )  bin += 25;
	  else if (bin < 250 )  bin += 35;
	  else if (bin < 550 )  bin += 45;
	  else                  bin += 55;
	  */
	  /*
	  if      (bin < 5)     bin += 5;
	  else if (bin < 10)    bin += 10;
	  else if (bin < 20)    bin += 15;
	  else if (bin < 30)    bin += 15;
	  else if (bin < 40)    bin += 15;
	  else if (bin < 100 )  bin += 20;
	  else if (bin < 250 )  bin += 50;
	  else if (bin < 550 )  bin += 60;
	  else                  bin += 70;
	  */

	  bin += 4;

       }
     }else{
	bin += 4;
	/*
       if (bin < 30)         bin += 1;
       else if (bin < 40)    bin += 1;
       else if (bin < 100 )  bin += 1;
       else if (bin < 250 )  bin += 1;
       else if (bin < 550 )  bin += 1;
       else                  bin += 1;
	*/
     }
     //hardcoded
     //cout<<"bin*binWidth is "<<bin*0.5<<endl;
   }
   bins[nbins] = xaxis->GetBinUpEdge(totBins);
   
   hNumEff = new TH1D("hNumEff","",nbins,bins);
   hDenEff = new TH1D("hDenEff","",nbins,bins);
   hNumFake = new TH1D("hNumFake","",nbins,bins);
   hNumFakeS = new TH1D("hNumFakeS","",nbins,bins);
   hDenFake = new TH1D("hDenFake","",nbins,bins);

   for (bin=1;bin<totBins;bin++){
     hNumEff->Fill(xaxis->GetBinCenter(bin),(hNEff->GetBinContent(bin)));
     hDenEff->Fill(xaxis->GetBinCenter(bin),(hDEff->GetBinContent(bin)));
     hNumFake->Fill(xaxis->GetBinCenter(bin),(hNFake->GetBinContent(bin)));
     hNumFakeS->Fill(xaxis->GetBinCenter(bin),(hDFake->GetBinContent(bin)));
     hDenFake->Fill(xaxis->GetBinCenter(bin),(hDFake->GetBinContent(bin)));
   }

   hNumEff->Sumw2();
   hDenEff->Sumw2();
   hNumFake->Sumw2();
   hDenFake->Sumw2();

   //hNumFakeS->Add(hNumFake,-1); // (hden-hnum)/hden => 1-(hnum/hden) => fake rate!
   //hNumFakeS->Sumw2();
   
   //hNum->Rebin(10);
   //hDen->Rebin(10);


   // calling                                                                                                           
   trkeffanalyzer_hist_data data;

   data.gEfficiency = new TGraphAsymmErrors();
   data.gFakerate =  new TGraphAsymmErrors();


   data.gEfficiency->BayesDivide(hNumEff,hDenEff);
   Double_t dYhigh=data.gEfficiency->GetErrorYhigh(0);
   Double_t dYlow=data.gEfficiency->GetErrorYlow(0);

   cout <<dYhigh<<" "<<dYlow<<" "<<endl;;
   data.gEfficiency->SetPointError(0,0.0,0.0,dYlow,dYhigh);

   data.gFakerate->BayesDivide(hNumFake,hDenFake);
   Double_t dYhighFR=data.gFakerate->GetErrorYhigh(0);
   Double_t dYlowFR=data.gFakerate->GetErrorYlow(0);

   cout <<dYhighFR<<" "<<dYlowFR<<" "<<endl;
   data.gFakerate->SetPointError(0,0.0,0.0,dYlowFR,dYhighFR);

   return (data);

}
