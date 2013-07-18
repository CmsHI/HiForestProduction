#ifndef GETEFFANDFAKE_H
#define GETEFFANDFAKE_H


#include <utility>
#include <iostream>

#include <TROOT.h>
#include <TStyle.h>
#include "TError.h"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TF1.h"

#include "TFile.h"
#include "TMath.h"

using namespace std;

//----------------------------------------------------------------------------------          
TFile *loadFile(TFile *file, char *cfile);
TFile *loadFile(TFile *file, char *cfile, char *option);
//----------------------------------------------------------------------------------          



TFile *loadFile(TFile *file, char *cfile, char *option){

   if(file!=0) delete file;
   file = new TFile(cfile,option);
   if(file->IsOpen()==false) Error("[FileLoading]","could not open hist file [%s].",cfile);
   else cout<<"[FileLoading] input hist file:"<<cfile<<" with an option "<<option<<endl;
   return file;

}

TFile *loadFile(TFile *file, char *cfile){

   //receive a pointer and return a pointer!                                                                                                                                                                  
   if(file!=0) delete file;
   file = new TFile(cfile);
   if(file->IsOpen()==false) Error("[FileLoading]","could not open hist file [%s].",cfile);
   else cout<<"[FileLoading] input hist file:"<<cfile<<endl;
   return file;

}

#endif
