#include <TFile.h>
#include <TTree.h>
#include <iostream>

void count_one(const char* fname)
{
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie())
  {
    std::cout << "NEVENTS -1" << std::endl;
    if (f) { f->Close(); delete f; }
    return;
  }

  TTree* T = (TTree*) f->Get("T");
  if (!T)
  {
    std::cout << "NEVENTS -1" << std::endl;
    f->Close();
    delete f;
    return;
  }

  std::cout << "NEVENTS " << T->GetEntries() << std::endl;

  f->Close();
  delete f;
}