/*
   This code takes the three FD MC swapfiles and divides them into 7 trees
   One tree for each PDG/swapfile combination (i.e Nutau from tauswap, Nue from nueswap, etc.)
   One tree for NC events
   Additionally adds a branch to each tree with the approximate oscillation probability and
      appropriate POT scaling.
*/


#include <iostream>
#include <cstdlib>

using namespace TMath;

const double DELMSQ_31 = 2.515e-3; //In eV^2
const double LOSC = 1300.; //In km

const double THETA_23 = 0.859;
const double THETA_13 = 0.150;

TString nonswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nonswap.root";
TString nueswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nueswap.root";
TString tauswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_tauswap.root";

double POTperYear = 1.1e21;
double nonswapPOT = 1.62824e24;
double nueswapPOT = 1.64546e24;
double tauswapPOT = 5.18551e24;

double OscWeight(const TString filename, const double Ev, const int nuPDG) {
  //First check flavor, then check which file it came from. Apply appropriate oscillation probability
  if (nuPDG == 12 || nuPDG == -12) {
    if (filename == nonswap)
      return 1. - Power(Sin(2 * THETA_13) * Sin((Pi() * DELMSQ_31 * LOSC) / (2.48 * Ev)), 2.); //Apply nue disappearance
    else if (filename == nueswap)
      return Power(Sin(THETA_23) * Sin(2 * THETA_13) * Sin((Pi() * DELMSQ_31 * LOSC) / (2.48 * Ev)), 2.); //Apply numu->nue
    else if (filename == tauswap) {
      std::cerr << "This should not happen: (nu_e in tauswap) \n";
      return 0.;
    } else {
      std::cerr << "This should not happen: (File name not recognized) \n";
        return 0.;
    }
  }
  else if (nuPDG == 14 || nuPDG == -14) {
    if (filename == nonswap) {
      return 1. - Power(Sin(2 * THETA_23) * Sin((Pi() * DELMSQ_31 * LOSC) / (2.48 * Ev)), 2.); //Apply numu disappearance
    }
    else if (filename == nueswap) {
      std::cerr << "This should not happen: (nu_mu in nueswap) \n";
      return 0.;
    }
    else if (filename == tauswap) {
      return Power(Sin(2 * THETA_13) * Sin(THETA_23) * Sin((Pi() * DELMSQ_31 * LOSC) / (2.48 * Ev)), 2.); //Apply nue->numu
    }
    else {
        std::cerr << "This should not happen: (File name not recognized) \n";
        return 0.;
    }
  }
  else if (nuPDG == 16 || nuPDG == -16) {
    if (filename == nonswap) {
      std:cerr << "This should not happen (nu_tau in nonswap) \n";
      return 0.;
    }
    else if (filename == nueswap) {
      return Power(Sin(2 * THETA_13) * Cos(THETA_23) * Sin((Pi() * DELMSQ_31 * LOSC) / (2.48 * Ev)), 2.); //Apply nue->nutau
    }
    else if (filename == tauswap) {
      return Power(Power(Cos(THETA_13), 2.) * Sin(2 * THETA_23) * Sin((Pi() * DELMSQ_31 * LOSC) / (2.48 * Ev)), 2.); //Apply numu->nutau
    }
    else {
        std::cerr << "This should not happen: (File name not recognized) \n";
        return 0.;
    }
  }
  else {
    std::cerr << "This should not happen: (nuPDG not recognized) \n";
    return 0.;
  }

}

void SplitMCTree(void) {

   //This will be the tree we train with
   TChain* ch1 = new TChain("cafTree");
   ch1->Add(nonswap);
   ch1->Add(nueswap);
   ch1->Add(tauswap);

   //This will be the tree we test on
   TChain* ch2 = new TChain("caf");
   ch2->Add(nonswap);
   ch2->Add(nueswap);
   ch2->Add(tauswap);

   //Use this file strictly to pull neutral current events
   TFile* NCFile = new TFile(nonswap, "OPEN");
   TTree* NCInputTree = (TTree*)NCFile->Get("cafTree");
   TTree* TestNCInputTree = (TTree*)NCFile->Get("caf");

   //Charged current events will be read from the chain

   TFile* outputFile = new TFile("/storage/shared/ncchambe/FDMonteCarlo/FHCSplitMC.root", "RECREATE");

   //Need 6 trees for CC: taus from tauswap, taus from nueswap, etc...
   TTree* NuTauFromTauswap = ch1->CloneTree(0);
   TTree* NuTauFromNueswap = ch1->CloneTree(0);
   TTree* NueFromNonswap = ch1->CloneTree(0);
   TTree* NueFromNueswap = ch1->CloneTree(0);
   TTree* NumuFromNonswap = ch1->CloneTree(0);
   TTree* NumuFromTauswap = ch1->CloneTree(0);
   TTree* NCTree = NCInputTree->CloneTree(0);

   TTree* TestNuTauFromTauswap = ch2->CloneTree(0);
   TTree* TestNuTauFromNueswap = ch2->CloneTree(0);
   TTree* TestNueFromNonswap = ch2->CloneTree(0);
   TTree* TestNueFromNueswap = ch2->CloneTree(0);
   TTree* TestNumuFromNonswap = ch2->CloneTree(0);
   TTree* TestNumuFromTauswap = ch2->CloneTree(0);
   TTree* TestNCTree = TestNCInputTree->CloneTree(0);


   int nuPDG = 0.;
   int isCC = 0.;
   double Ev = 0.;
   double POTScaledOscweight = 0.;
   ch1->SetBranchStatus("*", true);
   ch1->SetBranchAddress("nuPDG", &nuPDG);
   ch1->SetBranchAddress("isCC", &isCC);
   ch1->SetBranchAddress("Ev", &Ev);

   ch2->SetBranchStatus("*", true);
   ch2->SetBranchAddress("nuPDG", &nuPDG);
   ch2->SetBranchAddress("isCC", &isCC);
   ch2->SetBranchAddress("Ev", &Ev);

   NuTauFromTauswap->SetBranchStatus("*", true);
   NuTauFromNueswap->SetBranchStatus("*", true);
   NueFromNonswap->SetBranchStatus("*", true);
   NueFromNueswap->SetBranchStatus("*", true);
   NumuFromNonswap->SetBranchStatus("*", true);
   NumuFromTauswap->SetBranchStatus("*", true);
   NCTree->SetBranchStatus("*", true);
   NCTree->SetBranchAddress("isCC", &isCC);

   TestNuTauFromTauswap->SetBranchStatus("*", true);
   TestNuTauFromNueswap->SetBranchStatus("*", true);
   TestNueFromNonswap->SetBranchStatus("*", true);
   TestNueFromNueswap->SetBranchStatus("*", true);
   TestNumuFromNonswap->SetBranchStatus("*", true);
   TestNumuFromTauswap->SetBranchStatus("*", true);
   TestNCTree->SetBranchStatus("*", true);
   TestNCTree->SetBranchAddress("isCC", &isCC);

   NuTauFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NuTauFromNueswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NueFromNonswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NueFromNueswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NumuFromNonswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NumuFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NCTree->Branch("POTScaledOscweight", &POTScaledOscweight);

   TestNuTauFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   TestNuTauFromNueswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   TestNueFromNonswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   TestNueFromNueswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   TestNumuFromNonswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   TestNumuFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   TestNCTree->Branch("POTScaledOscweight", &POTScaledOscweight);

   auto trainentries = NCInputTree->GetEntries();
   POTScaledOscweight = POTperYear / nonswapPOT;
   for (auto i = 0; i < trainentries; i++) {
      NCInputTree->GetEntry(i);
      if (!isCC)
         NCTree->Fill();
   }

   auto testentries = TestNCInputTree->GetEntries();
   for (auto i = 0; i < testentries; i++) {
      TestNCInputTree->GetEntry(i);
      if (!isCC)
         TestNCTree->Fill();
   }


   auto nentries = ch1->GetEntries();
   for (auto i = 0; i < nentries; i++) {
      ch1->GetEntry(i);
      if (isCC && (nuPDG == 16 || nuPDG == -16) && ch1->GetCurrentFile()->GetName() == tauswap) {
         POTScaledOscweight = OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / tauswapPOT);
         NuTauFromTauswap->Fill();
      }
      else if (isCC && (nuPDG == 16 || nuPDG == -16) && ch1->GetCurrentFile()->GetName() == nueswap) {
         POTScaledOscweight = OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nueswapPOT);
         NuTauFromNueswap->Fill();
      }
      else if (isCC && (nuPDG == 14 || nuPDG == -14) && ch1->GetCurrentFile()->GetName() == tauswap) {
         POTScaledOscweight = OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / tauswapPOT);
         NumuFromTauswap->Fill();
      }
      else if (isCC && (nuPDG == 14 || nuPDG == -14) && ch1->GetCurrentFile()->GetName() == nonswap) {
         POTScaledOscweight = OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nonswapPOT);
         NumuFromNonswap->Fill();
      }
      else if (isCC && (nuPDG == 12 || nuPDG == -12) && ch1->GetCurrentFile()->GetName() == nonswap) {
         POTScaledOscweight = OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nonswapPOT);
         NueFromNonswap->Fill();
      }
      else if (isCC && (nuPDG == 12 || nuPDG == -12) && ch1->GetCurrentFile()->GetName() == nueswap) {
         POTScaledOscweight = OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nueswapPOT);
         NueFromNueswap->Fill();
      }
      if (i % 10000 == 0) std::cerr << (float)(i + 1) / (float)nentries << std::endl;
   }

   auto nentries2 = ch2->GetEntries();
   for (auto i = 0; i < nentries2; i++) {
      ch2->GetEntry(i);
      if (isCC && (nuPDG == 16 || nuPDG == -16) && ch2->GetCurrentFile()->GetName() == tauswap) {
         POTScaledOscweight = OscWeight(ch2->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / tauswapPOT);
         TestNuTauFromTauswap->Fill();
      }
      else if (isCC && (nuPDG == 16 || nuPDG == -16) && ch2->GetCurrentFile()->GetName() == nueswap) {
         POTScaledOscweight = OscWeight(ch2->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nueswapPOT);
         TestNuTauFromNueswap->Fill();
      }
      else if (isCC && (nuPDG == 14 || nuPDG == -14) && ch2->GetCurrentFile()->GetName() == tauswap) {
         POTScaledOscweight = OscWeight(ch2->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / tauswapPOT);
         TestNumuFromTauswap->Fill();
      }
      else if (isCC && (nuPDG == 14 || nuPDG == -14) && ch2->GetCurrentFile()->GetName() == nonswap) {
         POTScaledOscweight = OscWeight(ch2->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nonswapPOT);
         TestNumuFromNonswap->Fill();
      }
      else if (isCC && (nuPDG == 12 || nuPDG == -12) && ch2->GetCurrentFile()->GetName() == nonswap) {
         POTScaledOscweight = OscWeight(ch2->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nonswapPOT);
         TestNueFromNonswap->Fill();
      }
      else if (isCC && (nuPDG == 12 || nuPDG == -12) && ch2->GetCurrentFile()->GetName() == nueswap) {
         POTScaledOscweight = OscWeight(ch2->GetCurrentFile()->GetName(), Ev, nuPDG) * (POTperYear / nueswapPOT);
         TestNueFromNueswap->Fill();
      }
      if (i % 10000 == 0) std::cerr << (float)(i + 1) / (float)nentries2 << std::endl;
   }

   outputFile->WriteObject(NuTauFromTauswap, "NuTauFromTauswap");
   outputFile->WriteObject(NuTauFromNueswap, "NuTauFromNueswap");
   outputFile->WriteObject(NueFromNonswap, "NueFromNonswap");
   outputFile->WriteObject(NueFromNueswap, "NueFromNueswap");
   outputFile->WriteObject(NumuFromNonswap, "NumuFromNonswap");
   outputFile->WriteObject(NumuFromTauswap, "NumuFromTauswap");
   outputFile->WriteObject(NCTree, "NCTree");
   outputFile->WriteObject(TestNuTauFromTauswap, "TestNuTauFromTauswap");
   outputFile->WriteObject(TestNuTauFromNueswap, "TestNuTauFromNueswap");
   outputFile->WriteObject(TestNueFromNonswap, "TestNueFromNonswap");
   outputFile->WriteObject(TestNueFromNueswap, "TestNueFromNueswap");
   outputFile->WriteObject(TestNumuFromNonswap, "TestNumuFromNonswap");
   outputFile->WriteObject(TestNumuFromTauswap, "TestNumuFromTauswap");
   outputFile->WriteObject(TestNCTree, "TestNCTree");
   outputFile->Close();

   NCFile->Close();

}
