// =========================================================================
// Multi Tree Utility
// Author: Yen-Jie Lee (MIT)
//       A tool for multiTree handling
// =========================================================================
#include <TTree.h>
#include <TCut.h>
#include <TFile.h>
#include <vector.h>
#include <TH1.h>

class multiTreeUtil
{
   public:
      multiTreeUtil() {
         nTrees=0;
	 for (int i=0;i<100;i++)
	 {
	    color[i]=i;
	 }
	 color[0]=2;
	 color[1]=4;
	 color[2]=kGreen+2;
	 color[3]=6;
	 color[4]=7;
	 color[5]=kGreen+4;
	 color[6]=kBlue-2;
         color[7]=kAzure+3;
         color[8]=kOrange+3;
	 color[9]=kSpring-7;
	 color[10]=kPink+4;

	 /*
	   handsomeTH1(tmp1[1],2);
	   handsomeTH1(tmp1[2],4);
	   handsomeTH1(tmp1[3],kGreen+2);
	   handsomeTH1(tmp1[4],6);
	   handsomeTH1(tmp1[5],7);
	   handsomeTH1(tmp1[6],kGreen+4);
	   handsomeTH1(tmp1[7],kBlue-2);
	   handsomeTH1(tmp1[8],kAzure+3);
	   handsomeTH1(tmp1[9],kOrange+3);
	   handsomeTH1(tmp1[10],kSpring-7);
	   handsomeTH1(tmp1[11],kPink+4);
	 */
	 
      }
      ~multiTreeUtil() {
	 
      }
      void addFile(char *filename, char *treeName, TCut cut, Float_t scaleFactor);
      void Draw(TH1D* h,char* expression,char* cut = "" , char * opt="");
      void Draw2(TH1D* h,char* expression,char* cut = "", char * opt="");
      void Print(Option_t* option = "");
      TTree *getTree(int i) { return trees_[i]; };
      
 private:
      int nTrees;
      vector<TFile*>  tfiles_;                // TFile 
      vector<TTree*>  trees_;                 // Trees
      vector<char*>   tfileNames_;            // File names
      vector<TCut>    tcuts_;
      vector<Float_t> scaleFactors_;
      int color[100];

      void cleanHist(TH1* h);

};

// =========================================================================
// Add file
// =========================================================================
void multiTreeUtil::addFile(char *filename, char *treeName, TCut cut, Float_t scaleFactor)
{
   TFile *f = new TFile(filename);
   TTree *t = (TTree*) f->FindObjectAny(treeName);
   
   if (t!=0) {
      nTrees++;
      tfileNames_.push_back(filename);
      tfiles_.push_back(f);
      trees_.push_back(t);
      tcuts_.push_back(cut);
      scaleFactors_.push_back(scaleFactor);
   }
}

// =========================================================================
// Draw
// =========================================================================
void multiTreeUtil::Draw(TH1D *h, char *expression, char *cut, char *opt)
{

   h->Reset();
   
   for (int i=0;i<trees_.size();i++)
   {
      TH1D *htmp = (TH1D*)h->Clone();
      htmp->SetName(Form("my_htmp_%d",i));
      trees_[i]->Draw(Form("%s>>my_htmp_%d",expression,i),cut&&tcuts_[i]);
      htmp->Sumw2();
      htmp->Scale(scaleFactors_[i]);
      h->Add(htmp);
      delete htmp;
   }
   
   h->Draw(opt);
}

// =========================================================================
// Draw components
// =========================================================================
void multiTreeUtil::Draw2(TH1D *h, char *expression, char *cut, char *opt)
{
   TH1D *hComponent[100];
   
   h->Reset();
   
   for (int i=0;i<trees_.size();i++)
      {
      char *hName = Form("%s_Draw_%d",h->GetName(),i);
      TH1D *htmp= (TH1D*)h->Clone();
      htmp->SetName(hName);
      trees_[i]->Draw(Form("%s>>%s",expression,hName),cut&&tcuts_[i]);
      htmp->Sumw2();
      htmp->Scale(scaleFactors_[i]);
      h->Add(htmp);
      hComponent[i]=htmp;
      }
   
   h->Draw(opt);
   
   for (int i=0;i<trees_.size();i++)
      {
	 hComponent[i]->SetLineColor(color[i]);
         hComponent[i]->SetMarkerColor(color[i]);
	 hComponent[i]->Draw(Form("%s same",opt));
      }
}

// =========================================================================
// Print variables from the first tree
// =========================================================================
void multiTreeUtil::Print(Option_t* option)
{
   trees_[0]->Print(option);   
}

//--------------------------------PRIVATE-----------------------------------

// =========================================================================
// Clear
// =========================================================================
void multiTreeUtil::cleanHist(TH1* h)
{
   for (int i=0;i<=h->GetNbinsX()+1;i++) {
      h->SetBinContent(i,0);
      h->SetBinError(i,0);
   }
}
