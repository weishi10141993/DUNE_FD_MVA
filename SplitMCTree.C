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

TString FHCnonswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nonswap.root";
TString FHCnueswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nueswap.root";
TString FHCtauswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_tauswap.root";

TString RHCnonswap = "/storage/shared/wshi/CAFs/FD/FD_RHC_nonswap.root";
TString RHCnueswap = "/storage/shared/wshi/CAFs/FD/FD_RHC_nueswap.root";
TString RHCtauswap = "/storage/shared/wshi/CAFs/FD/FD_RHC_tauswap.root";

double POTperYear = 1.1e21;
double FHCnonswapPOT = 1.62824e24;
double FHCnueswapPOT = 1.64546e24;
double FHCtauswapPOT = 5.18551e24;

double RHCnonswapPOT = 3.27608e+24;
double RHCnueswapPOT = 3.24713e+24;
double RHCtauswapPOT = 8.58955e+24;

double scaleCorrection = 40. / 1.13;

double OscWeight(const TString filename, const double Ev, const int nuPDG, const double delmsq, const double theta_23, const double theta_13) {
  //First check flavor, then check which file it came from. Apply appropriate oscillation probability
  if (nuPDG == 12 || nuPDG == -12) {
      if (filename == FHCnonswap || filename == RHCnonswap) {
         return 1. - Power(Sin(2 * theta_13) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply nue disappearance
      }
      else if (filename == FHCnueswap || filename == RHCnueswap) {
         return Power(Sin(theta_23) * Sin(2 * theta_13) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply numu->nue
      }
      else if (filename == FHCtauswap || filename == RHCtauswap) {
         std::cerr << "This should not happen: (nu_e in tauswap) \n";
         return 0.;
      }
      else {
         std::cerr << "This should not happen: (File name not recognized) \n";
         return 0.;
      }
   }
   else if (nuPDG == 14 || nuPDG == -14) {
      if (filename == FHCnonswap || filename == RHCnonswap) {
         return 1. - Power(Sin(2 * theta_23) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply numu disappearance
      }
      else if (filename == FHCnueswap || filename == RHCnueswap) {
         std::cerr << "This should not happen: (nu_mu in nueswap) \n";
         return 0.;
      }
      else if (filename == FHCtauswap || filename == RHCtauswap) {
         return Power(Sin(2 * theta_13) * Sin(theta_23) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply nue->numu
      }
      else {
         std::cerr << "This should not happen: (File name not recognized) \n";
         return 0.;
      }
   }
   else if (nuPDG == 16 || nuPDG == -16) {
      if (filename == FHCnonswap || filename == RHCnonswap) {
         std:cerr << "This should not happen (nu_tau in nonswap) \n";
         return 0.;
      }
      else if (filename == FHCnueswap || filename == RHCnueswap) {
         return Power(Sin(2 * theta_13) * Cos(theta_23) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply nue->nutau
      }
      else if (filename == FHCtauswap || filename == RHCtauswap) {
         return Power(Power(Cos(theta_13), 2.) * Sin(2 * theta_23) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply numu->nutau
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

   //Change to true to use RHC files
   bool useRHC = true;
   TString nonswap = "";
   TString nueswap = "";
   TString tauswap = "";
   double nonswapPOT = 0.;
   double nueswapPOT = 0.;
   double tauswapPOT = 0.;

   if (!useRHC) {
      nonswap = FHCnonswap;
      nueswap = FHCnueswap;
      tauswap = FHCtauswap;

      nonswapPOT = FHCnonswapPOT;
      nueswapPOT = FHCnueswapPOT;
      tauswapPOT = FHCtauswapPOT;
   }
   else {
      nonswap = RHCnonswap;
      nueswap = RHCnueswap;
      tauswap = RHCtauswap;

      nonswapPOT = RHCnonswapPOT;
      nueswapPOT = RHCnueswapPOT;
      tauswapPOT = RHCtauswapPOT;
   }

   TString filename = ""; //This tracks if the event in the TChain is from non/nue/tauswap

   //This will be the tree we train with
   TChain* ch1 = new TChain("cafTree");
   ch1->Add(nonswap);
   ch1->Add(nueswap);
   ch1->Add(tauswap);


   //Use this file strictly to pull neutral current events
   TFile* NCFile = new TFile(nonswap, "OPEN");
   TTree* NCInputTree = (TTree*)NCFile->Get("cafTree");


   TFile* outputFile = NULL;
   if (!useRHC) {
      outputFile = new TFile("/storage/shared/ncchambe/FDMonteCarlo/FHCSplitMC.root", "RECREATE");
   }
   else {
      outputFile = new TFile("/storage/shared/ncchambe/FDMonteCarlo/RHCSplitMC.root", "RECREATE");
   }

   //Charged current events will be read from the chain
   //We clone the tchain but copy none of the entries to preserve the set of branches. We populate these trees below
   TTree* NuTauSignalFromTauswap = ch1->CloneTree(0);
   TTree* NuTauBkgdFromTauswap = ch1->CloneTree(0);
   TTree* NuTauFromNueswap = ch1->CloneTree(0);
   TTree* NueFromNonswap = ch1->CloneTree(0);
   TTree* NueFromNueswap = ch1->CloneTree(0);
   TTree* NumuFromNonswap = ch1->CloneTree(0);
   TTree* NumuFromTauswap = ch1->CloneTree(0);
   TTree* NCTree = NCInputTree->CloneTree(0);


   int nuPDG = 0.;
   int isCC = 0.;
   double Ev = 0.;
   double POTScaledOscweight = 0.;
   ch1->SetBranchStatus("*", true);
   ch1->SetBranchAddress("nuPDG", &nuPDG);
   ch1->SetBranchAddress("isCC", &isCC);
   ch1->SetBranchAddress("Ev", &Ev);

   NCInputTree->SetBranchStatus("*", true);
   NCInputTree->SetBranchAddress("nuPDG", &nuPDG);
   NCInputTree->SetBranchAddress("isCC", &isCC);
   NCInputTree->SetBranchAddress("Ev", &Ev);


   NuTauSignalFromTauswap->SetBranchStatus("*", true);
   NuTauBkgdFromTauswap->SetBranchStatus("*", true);
   NuTauFromNueswap->SetBranchStatus("*", true);
   NueFromNonswap->SetBranchStatus("*", true);
   NueFromNueswap->SetBranchStatus("*", true);
   NumuFromNonswap->SetBranchStatus("*", true);
   NumuFromTauswap->SetBranchStatus("*", true);
   NCTree->SetBranchStatus("*", true);


   NuTauSignalFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NuTauBkgdFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NuTauFromNueswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NueFromNonswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NueFromNueswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NumuFromNonswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NumuFromTauswap->Branch("POTScaledOscweight", &POTScaledOscweight);
   NCTree->Branch("POTScaledOscweight", &POTScaledOscweight);


   auto NCEntries = NCInputTree->GetEntries();
   POTScaledOscweight = scaleCorrection * POTperYear / nonswapPOT;
   for (auto i = 0; i < NCEntries; i++) {
      NCInputTree->GetEntry(i);
      if (!isCC) {
         NCTree->Fill();
      }
   }


   auto nentries = ch1->GetEntries();
   for (auto i = 0; i < nentries; i++) {
      ch1->GetEntry(i);
      filename = ch1->GetCurrentFile()->GetName();
      if (isCC) {
         if (Abs(nuPDG) == 16 && filename == tauswap) {
            POTScaledOscweight = OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scaleCorrection * POTperYear / tauswapPOT);
            if (!useRHC) {
               if (nuPDG == 16) NuTauSignalFromTauswap->Fill();
               else if (nuPDG == -16) NuTauBkgdFromTauswap->Fill();
            }
            else {
               if (nuPDG == -16) NuTauSignalFromTauswap->Fill();
               else if (nuPDG == 16) NuTauBkgdFromTauswap->Fill();
            }
         }
         else if (Abs(nuPDG) == 16 && filename == nueswap) {
            POTScaledOscweight = OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scaleCorrection * POTperYear / nueswapPOT);
            NuTauFromNueswap->Fill();
         }
         else if (Abs(nuPDG) == 14 && filename == tauswap) {
            POTScaledOscweight = OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scaleCorrection * POTperYear / tauswapPOT);
            NumuFromTauswap->Fill();
         }
         else if (Abs(nuPDG) == 14 && filename == nonswap) {
            POTScaledOscweight = OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scaleCorrection * POTperYear / nonswapPOT);
            NumuFromNonswap->Fill();
         }
         else if (Abs(nuPDG) == 12 && filename == nonswap) {
            POTScaledOscweight = OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scaleCorrection * POTperYear / nonswapPOT);
            NueFromNonswap->Fill();
         }
         else if (Abs(nuPDG) == 12 && filename == nueswap) {
            POTScaledOscweight = OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scaleCorrection * POTperYear / nueswapPOT);
            NueFromNueswap->Fill();
         }
      }
      if (i % 10000 == 0) {
          std::cerr << (float) i / (float)nentries << std::endl;
          std::cerr << "Current file: " << filename << std::endl;
      }
   }

   outputFile->WriteObject(NuTauSignalFromTauswap, "NuTauSignalFromTauswap");
   outputFile->WriteObject(NuTauBkgdFromTauswap, "NuTauBkgdFromTauswap");
   outputFile->WriteObject(NuTauFromNueswap, "NuTauFromNueswap");
   outputFile->WriteObject(NueFromNonswap, "NueFromNonswap");
   outputFile->WriteObject(NueFromNueswap, "NueFromNueswap");
   outputFile->WriteObject(NumuFromNonswap, "NumuFromNonswap");
   outputFile->WriteObject(NumuFromTauswap, "NumuFromTauswap");
   outputFile->WriteObject(NCTree, "NCTree");
   outputFile->Close();
   NCFile->Close();

}
