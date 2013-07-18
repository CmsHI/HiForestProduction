#include "multiTreeUtil.h"
#include <TCanvas.h>

void testMultiTree()
{
   multiTreeUtil myTrees;
   myTrees.addFile("mc/ntuple-HiMPA_photon15GeV.root","Analysis",TCut("ptHat<30"),1);
   myTrees.addFile("mc/ntuple-HiMPA_photon30GeV.root","Analysis",TCut("ptHat<50"),896./9517.);
   myTrees.addFile("mc/ntuple-HiMPA_photon50GeV.root","Analysis",TCut(""),896./9517.*1521./11551.);

   // Draw Leading Et
   TCanvas *c = new TCanvas("c","",600,600);
   TH1D *h = new TH1D("h","",100,10,110);
   myTrees.Draw2(h,"et[0]","");
   h->SetXTitle("E_{T} (GeV)");

   // Draw r9
   TCanvas *c2 = new TCanvas("c2","",600,600);
   TH1D *h2 = new TH1D("h2","",100,0,0.02);
   myTrees.Draw2(h2,"sigmaIetaIeta[0]","abs(eta[0])<1.479");
   h2->SetXTitle("#sigma_{i#eta i#eta}");

   // Draw ptHat
   TCanvas *c3 = new TCanvas("c3","",600,600);
   c3->SetLogy();
   TH1D *h3 = new TH1D("h3","",100,10,110);
   myTrees.Draw(h3,"ptHat","");
   h3->SetXTitle("#hat{p_{T}} (GeV/c)");

   // Print the first tree
   myTrees.Print();

   // Get the first tree and draw something
   TCanvas *c4 = new TCanvas("c4","",600,600);
   myTrees.getTree(0)->Draw("et");
}
