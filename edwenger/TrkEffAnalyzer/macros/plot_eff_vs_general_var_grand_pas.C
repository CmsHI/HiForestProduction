//

#include "commonStyle.h"
Int_t color[12]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,1,1,1};

void plot_eff_vs_general_var_grand_pas(Bool_t rebin=false, Bool_t range=false){

   //-------------------------------------------------------------------------------------
   // variable: 1 - eta, 2 - pt, 3 - jet 
   //-------------------------------------------------------------------------------------
   

   int variable = 2;
   rebin = false;
   range = true;
   
   bool MBonly = false;
   bool effonly = true;
   bool fakonly = false;
   bool lowerhalf = false;
   int logX = 1;

   double minpt = 5.0;
   double maxpt = 200;

   double minjet = 0;
   double maxjet = 0;
   
   double minjet2 = 130;
   double maxjet2 = 170;



   double maxeta = 1.0;

   Char_t VARIABLE[100];


   if(variable==1){
      TH1F *dum = new TH1F("dum","",100,-2.5,2.5);
      dum->GetXaxis()->SetTitle("#eta");
      sprintf(VARIABLE,"ETA_%1.1f",maxeta);
      if(lowerhalf) dum->SetMaximum(0.5);
   }else if(variable==2){
      TH1F *dum = new TH1F("dum","",100,0.5,320);
      //if(range) dum->SetAxisRange(0,160);
      dum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      dum->GetXaxis()->CenterTitle();
      dum->GetYaxis()->SetDecimals();
      sprintf(VARIABLE,"PT_minpt_%1.1f_%1.1f_%1.1f",minpt,minjet,maxjet);
   }else if(variable==3){
      TH1F *dum = new TH1F("dum","",100,0,1200);
      dum->GetXaxis()->SetTitle("Corrected leading jet E_{T} (GeV)");
      sprintf(VARIABLE,"JETET_minpt_%1.1f",minpt);
   }   


   dum->SetMinimum(-0.05);
   dum->SetMaximum(1.0);
   dum->SetStats(0);
   //dum->GetYaxis()->SetTitle("Fractional efficiency");
   dum->GetYaxis()->SetTitle("A #times #epsilon_{TR}");
   dum->GetXaxis()->SetTitleSize(0.05);
   dum->GetYaxis()->SetTitleSize(0.05);
   dum->GetYaxis()->SetTitleOffset(1.5);
   dum->GetYaxis()->SetNdivisions(911);
   
   //gROOT->LoadMacro("trkeffanalyzer_hist_jet.C");      
   gROOT->LoadMacro("trkeffanalyzer_hist_general_variable.C");
   
   Char_t fileName1[100];
   Char_t fileName2[100];
   Char_t fileName3[100];
   Char_t fileName4[100];
   Char_t fileName5[100];
   Char_t fileName6[100];

   
   Char_t fileString1[200];
   Char_t fileString2[200];
   Char_t fileString3[200];
   Char_t fileString4[200];
   Char_t fileString5[200];
   Char_t fileString6[200];


   //sprintf(fileName1,"TrkHistMC_QCD_10M.root");
   //sprintf(fileName1,"TrkHistMC_july09v2_qcdMBD6TdJuly09V3_div2.root");
   //sprintf(fileString1,"/Users/andreyoon/Research/ana/spectra/pp_spectra/7TeV/root_files/%s",fileName1);

   sprintf(fileName1,"TrkHistGEN_oct06_qcdMBD6TV1.root");
   sprintf(fileString1,"/Users/andreyoon/Research/pp7TeV/root_files_postApp/mc/%s",fileName1);
   cout<<"input file 1 is "<<fileString1<<endl;
   
   //sprintf(fileName2,"TrkHistMC_QCD_Pt15.root");
   //sprintf(fileName2,"TrkHistMC_july09v2_qcdPt15dJuly09_div2.root");
   //sprintf(fileString2,"/Users/andreyoon/Research/ana/spectra/pp_spectra/7TeV/root_files/%s",fileName2);

   sprintf(fileName2,"TrkHistGEN_oct01_qcdPt15V1.root");
   sprintf(fileString2,"/Users/andreyoon/Research/pp7TeV/root_files_postApp/mc/%s",fileName2);
   cout<<"input file 2 is "<<fileString2<<endl;
   
   //sprintf(fileName3,"TrkHistMC_QCD_Pt30v2.root");
   //sprintf(fileName3,"TrkHistMC_july09v2_qcdPt30dJuly09V2_div2.root");
   //sprintf(fileString3,"/Users/andreyoon/Research/ana/spectra/pp_spectra/7TeV/root_files/%s",fileName3);

   sprintf(fileName3,"TrkHistGEN_oct01_qcdPt30V2.root");
   sprintf(fileString3,"/Users/andreyoon/Research/pp7TeV/root_files_postApp/mc/%s",fileName3);
   cout<<"input file 3 is "<<fileString3<<endl;
   
   //sprintf(fileName4,"TrkHistMC_QCD_Pt80v2.root");
   //sprintf(fileName4,"TrkHistMC_july09v2_qcdPt80dJuly09V2_div2.root");
   //sprintf(fileString4,"/Users/andreyoon/Research/ana/spectra/pp_spectra/7TeV/root_files/%s",fileName4);

   sprintf(fileName4,"TrkHistGEN_oct01_qcdPt80V2.root");
   sprintf(fileString4,"/Users/andreyoon/Research/pp7TeV/root_files_postApp/mc/%s",fileName4);
   cout<<"input file 4 is "<<fileString4<<endl;
   
   //sprintf(fileName5,"TrkHistMC_QCD_Pt170_V2_1st.root");  
   //sprintf(fileName5,"TrkHistMC_july09v2_qcdPt170dJuly09_div2.root");
   //sprintf(fileString5,"/Users/andreyoon/Research/ana/spectra/pp_spectra/7TeV/root_files/%s",fileName5);

   sprintf(fileName5,"TrkHistGEN_oct01_qcdPt170V1.root");
   sprintf(fileString5,"/Users/andreyoon/Research/pp7TeV/root_files_postApp/mc/%s",fileName5);
   cout<<"input file 5 is "<<fileString5<<endl;

   //sprintf(fileName6,"TrkHistMC_QCD_Pt470_V2_1st.root"); 
   //sprintf(fileName6,"TrkHistMC_july09v2_dec14_qcdPt300dJuly09V2.root");
   //sprintf(fileString6,"/Users/andreyoon/Research/ana/spectra/pp_spectra/7TeV/root_files/%s",fileName6);

   sprintf(fileName6,"TrkHistGEN_oct01_qcdPt300V1.root");
   sprintf(fileString6,"/Users/andreyoon/Research/pp7TeV/root_files_postApp/mc/%s",fileName6);
   cout<<"input file 6 is "<<fileString6<<endl;

   char jetrange1[100], jetrange2[100], jetrange3[100], jetrange4[100], jetrange5[100], jetrange6[100];
   char jetrange7[100], jetrange8[100], jetrange9[100], jetrange10[100], jetrange11[100], jetrange12[100];

   // 0-20  MB

   minjet = 0, maxjet = 20;
   sprintf(jetrange1," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   trkeffanalyzer_hist_data d1 = trkeffanalyzer_hist_graph(fileString1,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff1 = d1.gEfficiency;
   TGraphAsymmErrors *FR1  = d1.gFakerate;

   TGraphAsymmErrors *Eff1Skm = removeLastPoint(Eff1,1);
   TGraphAsymmErrors *FR1Skm = removeLastPoint(FR1,1);

   // 20-40 MB, Pt15
   minjet = 20, maxjet =40;
   sprintf(jetrange2," %1.0f<E_{T}<%1.0f",minjet,maxjet);
   
   trkeffanalyzer_hist_data d2 = trkeffanalyzer_hist_graph(fileString1,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff2 = d2.gEfficiency;
   TGraphAsymmErrors *FR2  = d2.gFakerate;

   TGraphAsymmErrors *Eff2Skm = removeLastPoint(Eff2,1);
   TGraphAsymmErrors *FR2Skm = removeLastPoint(FR2,1);

   //TGraphAsymmErrors *Eff2Skm = removeLastPoint(Eff2,0);
   //TGraphAsymmErrors *FR2Skm = removeLastPoint(FR2,0);

   trkeffanalyzer_hist_data d3 = trkeffanalyzer_hist_graph(fileString2,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff3 = d3.gEfficiency;
   TGraphAsymmErrors *FR3  = d3.gFakerate;

   TGraphAsymmErrors *Eff3Skm = removeLastPoint(Eff3,1);
   TGraphAsymmErrors *FR13km = removeLastPoint(FR3,1);

   //TGraphAsymmErrors *Eff3Skm = removeLastPoint(Eff3,0);
   //TGraphAsymmErrors *FR13km = removeLastPoint(FR3,0);

   // 40-60 Pt15, Pt30
   minjet = 40, maxjet =60;
   sprintf(jetrange3," %1.0f<E_{T}<%1.0f ",minjet,maxjet);

   trkeffanalyzer_hist_data d4 = trkeffanalyzer_hist_graph(fileString2,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff4 = d4.gEfficiency;
   TGraphAsymmErrors *FR4  = d4.gFakerate;

   TGraphAsymmErrors *Eff4Skm = removeLastPoint(Eff4,2);
   TGraphAsymmErrors *FR4Skm = removeLastPoint(FR4,2);

   //TGraphAsymmErrors *Eff4Skm = removeLastPoint(Eff4,0);
   //TGraphAsymmErrors *FR4Skm = removeLastPoint(FR4,0);

   trkeffanalyzer_hist_data d5 = trkeffanalyzer_hist_graph(fileString3,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff5 = d5.gEfficiency;
   TGraphAsymmErrors *FR5 = d5.gFakerate;
 
   TGraphAsymmErrors *Eff5Skm = removeLastPoint(Eff5,2);
   TGraphAsymmErrors *FR5Skm = removeLastPoint(FR5,2);

   //TGraphAsymmErrors *Eff5Skm = removeLastPoint(Eff5,0);
   //TGraphAsymmErrors *FR5Skm = removeLastPoint(FR5,0);

   // 60-80 Pt30 
   minjet = 60, maxjet =80;
   sprintf(jetrange4," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   trkeffanalyzer_hist_data d6 = trkeffanalyzer_hist_graph(fileString3,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff6 = d6.gEfficiency;
   TGraphAsymmErrors *FR6 = d6.gFakerate;

   TGraphAsymmErrors *Eff6Skm = removeLastPoint(Eff6,3);
   TGraphAsymmErrors *FR6Skm = removeLastPoint(FR6,3);


   // 80-100 Pt80 -> Pt30
   minjet = 80, maxjet =100;
   sprintf(jetrange5," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   //trkeffanalyzer_hist_data d7 = trkeffanalyzer_hist_graph(fileString4,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);
   trkeffanalyzer_hist_data d7 = trkeffanalyzer_hist_graph(fileString3,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt); 
   
   TGraphAsymmErrors *Eff7 = d7.gEfficiency;
   TGraphAsymmErrors *FR7 = d7.gFakerate;
   
   TGraphAsymmErrors *Eff7Skm = removeLastPoint(Eff7,2);
   TGraphAsymmErrors *FR7Skm = removeLastPoint(FR7,2);

   //TGraphAsymmErrors *Eff7Skm = removeLastPoint(Eff7,0);
   //TGraphAsymmErrors *FR7Skm = removeLastPoint(FR7,0);

 
   // 100-120 Pt80
   minjet = 100, maxjet =120;
   sprintf(jetrange6," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   trkeffanalyzer_hist_data d8 = trkeffanalyzer_hist_graph(fileString4,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff8 = d8.gEfficiency;
   TGraphAsymmErrors *FR8 = d8.gFakerate;

   TGraphAsymmErrors *Eff8Skm = removeLastPoint(Eff8,3);
   TGraphAsymmErrors *FR8Skm = removeLastPoint(FR8,3);

   //TGraphAsymmErrors *Eff8Skm = removeLastPoint(Eff8,0);
   //TGraphAsymmErrors *FR8Skm = removeLastPoint(FR8,0);


   // 120-140
   
   // 140-160

   // 160-180

   // 180-200 Pt80, Pt170
   minjet = 180, maxjet =200;
   sprintf(jetrange7," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   trkeffanalyzer_hist_data d9 = trkeffanalyzer_hist_graph(fileString4,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff9 = d9.gEfficiency;
   TGraphAsymmErrors *FR9 = d9.gFakerate;
   
   TGraphAsymmErrors *Eff9Skm = removeLastPoint(Eff9,1);
   TGraphAsymmErrors *FR9Skm = removeLastPoint(FR9,1);

   //TGraphAsymmErrors *Eff9Skm = removeLastPoint(Eff9,0);
   //TGraphAsymmErrors *FR9Skm = removeLastPoint(FR9,0);


   trkeffanalyzer_hist_data d10 = trkeffanalyzer_hist_graph(fileString5,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff10 = d10.gEfficiency;
   TGraphAsymmErrors *FR10 = d10.gFakerate;
   
   TGraphAsymmErrors *Eff10Skm = removeLastPoint(Eff10,2);
   TGraphAsymmErrors *FR10Skm = removeLastPoint(FR10,2);

   //TGraphAsymmErrors *Eff10Skm = removeLastPoint(Eff10,0);
   //TGraphAsymmErrors *FR10Skm = removeLastPoint(FR10,0);

   //250 -300 
   minjet = 250, maxjet =300;
   sprintf(jetrange8," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   trkeffanalyzer_hist_data d11 = trkeffanalyzer_hist_graph(fileString4,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff11 = d11.gEfficiency;
   TGraphAsymmErrors *FR11 = d11.gFakerate;

   TGraphAsymmErrors *Eff11Skm = removeLastPoint(Eff11,1);
   TGraphAsymmErrors *FR11Skm = removeLastPoint(FR11,1);

   //TGraphAsymmErrors *Eff11Skm = removeLastPoint(Eff11,0);
   //TGraphAsymmErrors *FR11Skm = removeLastPoint(FR11,0);
   

   //500 - 2400
   minjet = 500, maxjet =2400;
   sprintf(jetrange9," %1.0f<E_{T}<%1.0f",minjet,maxjet);

   trkeffanalyzer_hist_data d12 = trkeffanalyzer_hist_graph(fileString6,rebin,range,variable,0.0,maxeta,minjet,maxjet,minpt,maxpt);

   TGraphAsymmErrors *Eff12 = d12.gEfficiency;
   TGraphAsymmErrors *FR12 = d12.gFakerate;

   TGraphAsymmErrors *Eff12Skm = removeLastPoint(Eff12,0);
   TGraphAsymmErrors *FR12Skm = removeLastPoint(FR12,0);


   

   Char_t xTitle[100],yTitle[100];
   sprintf(xTitle,"p_{T} [GeV/c]");
   sprintf(yTitle,"Ratio");
   
   TCanvas *call = new TCanvas("call","call",500,550);
   call->SetGridx(),call->SetGridy();
   //if(logX) call->SetLogx();
   call->cd();
   dum->Draw();
   
   float binsize = 1.3;
   float binsize2 = binsize*0.86;
   

   int color=0;
   color = 15;
   int dopt=4;

   th1Style1(Eff1Skm,color,20,binsize,color,1,1,dopt);       
   th1Style1(FR1Skm,color,21,binsize,color,1,1,dopt);

   color = 16;
   th1Style1(Eff2Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR2Skm,color,21,binsize,color,1,1,dopt);

   //th1Style1(Eff3Skm,color,24,binsize,color,1,1,dopt);
   //th1Style1(FR3Skm,color,25,binsize,color,1,1,dopt);

   color = 17; 
   th1Style1(Eff4Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR4Skm,color,21,binsize,color,1,1,dopt);

   //th1Style1(Eff5,color,24,binsize,color,1,1,1);
   //th1Style1(FR5,color,25,binsize,color,1,1,1);

   color = 18;
   th1Style1(Eff6Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR6Skm,color,21,binsize,color,1,1,dopt);

   color = 19;
   th1Style1(Eff7Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR7Skm,color,21,binsize,color,1,1,dopt);

   color = 20;
   th1Style1(Eff8Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR8Skm,color,21,binsize,color,1,1,dopt);

   color = 21;
   //th1Style1(Eff9,color,20,binsize,color,1,1,1);
   //th1Style1(FR9,color,21,binsize,color,1,1,1);

   th1Style1(Eff10Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR10Skm,color,21,binsize,color,1,1,dopt);

   color = 22;
   //th1Style1(Eff11Skm,color,20,binsize,color,1,1,dopt);
   //th1Style1(FR11Skm,color,21,binsize,color,1,1,dopt);

   th1Style1(Eff12Skm,color,20,binsize,color,1,1,dopt);
   th1Style1(FR12Skm,color,21,binsize,color,1,1,dopt);

   

   /*
   th1Style1(FR3,15,4,binsize,15,1,1,1);
   TGraphAsymmErrors* FR3_dum = FR3->Clone("FR3_dum");
   th1Style1(FR3_dum,10,20,binsize2,15,1,1,1);

   if(!MBonly){
   th1Style1(Eff4,16,20,binsize,16,1,1,1);
   
   th1Style1(FR4,16,4,binsize,16,1,1,1);
   TGraphAsymmErrors* FR4_dum = FR4->Clone("FR4_dum");
   th1Style1(FR4_dum,10,20,binsize2,16,1,1,1);
   
   
   th1Style1(Eff5,17,20,binsize,17,1,1,1);
   th1Style1(FR5,17,4,binsize,17,1,1,1);
   TGraphAsymmErrors* FR5_dum = FR5->Clone("FR5_dum");
   th1Style1(FR5_dum,10,20,binsize2,17,1,1,1);
   
   
   th1Style1(Eff6,15,24,binsize,15,1,1,1);
   th1Style1(FR6,18,4,binsize,18,1,1,1);
   TGraphAsymmErrors* FR6_dum = FR6->Clone("FR6_dum");
   th1Style1(FR6_dum,10,20,binsize2,18,1,1,1);

   th1Style1(Eff7,16,24,binsize,16,1,1,1);
   th1Style1(FR7,19,4,binsize,19,1,1,1);
   TGraphAsymmErrors* FR7_dum = FR7->Clone("FR7_dum");
   th1Style1(FR7_dum,10,20,binsize2,19,1,1,1);

   th1Style1(Eff8,17,24,binsize,17,1,1,1);
   th1Style1(FR8,20,4,binsize,20,1,1,1);
   TGraphAsymmErrors* FR8_dum = FR7->Clone("FR8_dum");
   th1Style1(FR8_dum,10,20,binsize2,19,1,1,1);
   }
   */


   if(!MBonly){
      //TLegend *leg = new TLegend(0.54,0.4,0.89,0.62);
      //TLegend *leg = new TLegend(0.54,0.36,0.89,0.62);      
      //TLegend *leg = new TLegend(0.19,0.37,0.93,0.63); 
      //TLegend *leg = new TLegend(0.19,0.4,0.65,0.63);     
     //TLegend *leg = new TLegend(0.23,0.3,0.69,0.52);  
     TLegend *leg = new TLegend(0.23,0.25,0.69,0.47);   
      leg->SetBorderSize(1);
      leg->SetFillColor(kWhite);
      //leg->SetFillStyle(0);       
      leg->SetTextSize(0.031);
      leg->SetMargin(0.55);
      //leg->SetMargin(0.5);   
      /*
      char jetrange1[100];
      sprintf(jetrange1," %1.0f<E_{T}<%1.0f",minjet,maxjet); 
      char jetrange2[100];
      sprintf(jetrange2," %1.0f<E_{T}<%1.0f",minjet2,maxjet2);
      char header[100];
      sprintf(header,"    p_{T}^{hat}>15   p_{T}^{hat}>30   p_{T}^{hat}>80");
      leg->SetHeader(header); 
      */
      leg->SetNColumns(2);
      leg->SetColumnSeparation(0.000);

      leg->AddEntry(Eff1Skm,jetrange1,"lp");
      leg->AddEntry(dum,"MB   ","");

      leg->AddEntry(Eff2Skm,jetrange2,"lp");
      leg->AddEntry(dum,"MB   ","");

      leg->AddEntry(Eff4Skm,jetrange3,"lp");
      leg->AddEntry(dum,"Pt15  ","");

      leg->AddEntry(Eff6Skm,jetrange4,"lp");
      leg->AddEntry(dum,"Pt30  ","");

      leg->AddEntry(Eff7Skm,jetrange5,"lp");
      leg->AddEntry(dum,"Pt30  ","");

      leg->AddEntry(Eff8Skm,jetrange6,"lp");
      leg->AddEntry(dum,"Pt80  ","");

      leg->AddEntry(Eff10Skm,jetrange7,"lp");
      leg->AddEntry(dum,"Pt170  ","");
      
      leg->AddEntry(Eff12Skm,jetrange9,"lp");
      leg->AddEntry(dum,"Pt300  ","");

   }else{
     TLegend *leg = new TLegend(0.54,0.4,0.89,0.62);
     leg->SetBorderSize(1);
     leg->SetFillColor(kWhite);
     leg->SetTextSize(0.031);
     leg->SetMargin(0.35);
     //leg->SetHeader("  0 < Jet E_{T}< 20 ");
     leg->AddEntry(Eff3,"PYTHIA 7 TeV","lp");
   }
   leg->Draw();
   

   TLatex *tex = new TLatex(0.19,0.20,"Fake rate");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   //tex->Draw();


   
     Char_t outName1[100];
     //sprintf(outName1,"Eff_FR_GRAND_%s",VARIABLE);
     sprintf(outName1,"Eff_FR_GRAND_wider_v2");
     Char_t outName2[100];
     //sprintf(outName2,"Eff_FR_GRAND_logx_%s",VARIABLE);
     sprintf(outName2,"Eff_FR_GRAND_logx_wider_v2");


     tex->Draw();
     dum->SetMaximum(0.91);
     //dum->SetMinimum(0.65);
     dum->SetMinimum(0.50);
     printCanvases(call,outName1,0);

     dum->SetMinimum(-0.05);
     dum->SetMaximum(1.0);
     printCanvases(call,outName2,1);

}


TGraphAsymmErrors* ratio_tgraphs(TGraphAsymmErrors *n, TGraphAsymmErrors *d, double scale){

   cout<<"[ratio_tgraphs]"<<endl;
 
   if(!n || !d) return 0;
   int N = d->GetN();

   cout<<"number of bins "<<N<<endl;

   TGraphAsymmErrors *ratio = new TGraphAsymmErrors(N);
   char name[300];
   sprintf(name,"ratio_%s_over_%s",n->GetName(),d->GetName());
   
   ratio->SetName(name);

   for(int i=0;i<N;i++)
      {
	 double nx, ny;
	 double nex, ney;
	 
	 n->GetPoint(i,nx,ny);
	 nex = n->GetErrorX(i);
	 ney = n->GetErrorY(i);
	 
	 double dx, dy;
	 double dex, dey;
	 
	 d->GetPoint(i,dx,dy);
	 //cout<<" dy "<<dy<<endl;
	 dex = d->GetErrorX(i);
	 dey = d->GetErrorY(i);
	 
	 double rx, ry;
	 double rex, rey;
	 
	 //cout<<"nx :"<<nx<<"ny :"<<ny<<endl;
	 //cout<<"dx :"<<dx<<"dy :"<<dy<<endl;

	 if(dy==0) ry = 999;
	 else ry = (ny/dy)*scale;
	 //rex = rey = 0;
	 rex = 0;
	 //rey = sqrt((ney*ney)/(ny*ny)+(dey*dey)/(dy*dy));
	 rey = sqrt((ney*ney)+(dey*dey))*scale;          

	 ratio->SetPoint(i,dx,ry);
	 ratio->SetPointError(i,rex,rex,rey,rey);
   }

   return ratio;
}



TH1D* GetDummyHist(Float_t min, Float_t max,Char_t *xttl,Char_t *yttl) {

   TH1D *dum;
   dum = new TH1D("dum","",100,0.0,80);

   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(0);

   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetYaxis()->SetDecimals();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();

   return dum;

}


TGraphAsymmErrors *removeLastPoint(TGraphAsymmErrors *tgi, int nskip){
  TGraphAsymmErrors *tgo=0;

  int nbin = tgi->GetN();
  cout<<"nbin = "<<nbin<<endl;
  const int nlines = nbin-nskip;
  //const int nlines = nbin;
  double pt[nlines], pterrL[nlines], pterrH[nlines];
  double eff[nlines], efferrL[nlines], efferrH[nlines]; 

  for(int i=0;i<nlines;i++){
    tgi->GetPoint(i,pt[i],eff[i]);
    efferrL[i]=tgi->GetErrorYlow(i);
    efferrH[i]=tgi->GetErrorYhigh(i);
    pterrL[i]=pterrH[i]=0;
  }
  tgo = new TGraphAsymmErrors(nlines,pt,eff,pterrL,pterrH,efferrL,efferrH);
  return tgo;
}
