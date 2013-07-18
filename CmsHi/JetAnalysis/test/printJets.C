

void printJets(const char* file = "jectest_sqlite.root"){

  TFile* inf = new TFile(file);

  TTree* t = (TTree*)inf->Get("ak3PFJetAnalyzer/t");
  TCut realjets("refpt > 0");

  t->Scan("refpt:rawpt:jtpt:corrpt",realjets,"",100);


}

