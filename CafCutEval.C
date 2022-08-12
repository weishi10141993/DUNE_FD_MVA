/*
   This code tests the sensitivity, efficiency, and purity of set of neutrino measurements
      given a particular cut.
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

//This function is what I used to calculate the POT numbers for each MC swapfile
void POTNumber(void) {


  TChain* ch1 = new TChain("meta");
  ch1->Add(FHCnonswap);
  ch1->Add(FHCnueswap);
  ch1->Add(FHCtauswap);

  int nentries = ch1->GetEntries();

  double nonswapPOT = 0;
  double nueswapPOT = 0;
  double tauswapPOT = 0;


  double POT = 0.;
  ch1->SetBranchStatus("*", kFALSE);
  ch1->SetBranchStatus("pot", kTRUE);

  ch1->SetBranchAddress("pot", &POT);

  for (auto i = 0; i < nentries; i++) {
    ch1->GetEntry(i);
    if (ch1->GetCurrentFile()->GetName() == RHCnonswap)
      nonswapPOT += POT;
    else if (ch1->GetCurrentFile()->GetName() == RHCnueswap)
      nueswapPOT += POT;
    else if (ch1->GetCurrentFile()->GetName() == RHCtauswap)
      tauswapPOT += POT;
    else
      std::cerr << "This should not happen (file name not recognized) \n";
  }

  std::cout << "nonswapPOT = " << nonswapPOT << '\n';
  std::cout << "nueswapPOT = " << nueswapPOT << '\n';
  std::cout << "tauswapPOT = " << tauswapPOT << '\n';
}


void CafCutEval(void) {

   bool useRHC = false;
   bool useEHadBinning = false;

   double nonswapPOT = 0.;
   double nueswapPOT = 0.;
   double tauswapPOT = 0.;

   TString nonswap = "";
   TString nueswap = "";
   TString tauswap = "";

   if (!useRHC) {
      nonswapPOT = FHCnonswapPOT;
      nueswapPOT = FHCnueswapPOT;
      tauswapPOT = FHCtauswapPOT;

      nonswap = FHCnonswap;
      nueswap = FHCnueswap;
      tauswap = FHCtauswap;

   }
   else {
      nonswapPOT = RHCnonswapPOT;
      nueswapPOT = RHCnueswapPOT;
      tauswapPOT = RHCtauswapPOT;

      nonswap = RHCnonswap;
      nueswap = RHCnueswap;
      tauswap = RHCtauswap;

   }

   double scalefactor = POTperYear * scaleCorrection * 3.5; //3.5 years of data taking


   gStyle->SetOptStat(0);
   TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
   c1->SetRightMargin(0.1);
   c1->SetLeftMargin(0.15);


   TChain* ch1 = new TChain("cafTree");


   ch1->Add(nonswap);
   ch1->Add(nueswap);
   ch1->Add(tauswap);

   /*
      These histograms are from when I made plots of efficiency/significance/purity as a function of energy
   TH1D hPurity = TH1D("hPurity", "hPurity", 25, 0., 10.);
   TH1D hEfficiency = TH1D("hEfficiency", "hEfficiency", 25, 0., 10.);
   TH1D hSignificance = TH1D("hSignificance", "hSignificance", 25, 0., 5.);

   TH1D* hTrueSignal = new TH1D("hTrueSignal", "hTrueSignal", 25, 0., 10.);
   TH1D* hTrueSignalCut = new TH1D("hTrueSignalCut", "hTrueSignalCut", 25, 0., 10.);
   TH1D* hCutCounter = new TH1D("hCutCounter", "hCutCounter", 25, 0., 10.);
   */


   double Ev = 0.;
   int nuPDG = 0.;
   int isCC = -1.;
   double cvnnue = 0.;
   double cvnnumu = 0.;
   double cvnnutau = 0.;
   double cvnnc = 0.;

   //Reconstructed Ev
   double eDepP = 0.;
   double eDepN = 0.;
   double eDepPip = 0.;
   double eDepPim = 0.;
   double eDepPi0 = 0.;
   double eDepOther = 0.;
   double eDepTotal = 0.;
   double DepPimRatio = 0.;
   double RecoLepEnNumu = 0.;
   double RecoLepEnNue = 0.;
   double eDepHad = 0.;

   double vtx_x = 0.;
   double vtx_y = 0.;
   double vtx_z = 0.;


   ch1->SetBranchStatus("*", false);
   ch1->SetBranchStatus("Ev", true);
   ch1->SetBranchAddress("Ev", &Ev);

   ch1->SetBranchStatus("nuPDG", true);
   ch1->SetBranchAddress("nuPDG", &nuPDG);

   ch1->SetBranchStatus("isCC", true);
   ch1->SetBranchAddress("isCC", &isCC);

   ch1->SetBranchStatus("cvnnue", true);
   ch1->SetBranchAddress("cvnnue", &cvnnue);

   ch1->SetBranchStatus("cvnnumu", true);
   ch1->SetBranchAddress("cvnnumu", &cvnnumu);

   ch1->SetBranchStatus("cvnnutau", true);
   ch1->SetBranchAddress("cvnnutau", &cvnnutau);

   ch1->SetBranchStatus("cvnnc", true);
   ch1->SetBranchAddress("cvnnc", &cvnnc);

   ch1->SetBranchStatus("eDepP", true);
   ch1->SetBranchAddress("eDepP", &eDepP);

   ch1->SetBranchStatus("eDepN", true);
   ch1->SetBranchAddress("eDepN", &eDepN);

   ch1->SetBranchStatus("eDepPip", true);
   ch1->SetBranchAddress("eDepPip", &eDepPip);

   ch1->SetBranchStatus("eDepPim", true);
   ch1->SetBranchAddress("eDepPim", &eDepPim);

   ch1->SetBranchStatus("eDepPi0", true);
   ch1->SetBranchAddress("eDepPi0", &eDepPi0);

   ch1->SetBranchStatus("eDepOther", true);
   ch1->SetBranchAddress("eDepOther", &eDepOther);

   ch1->SetBranchStatus("RecoLepEnNumu", true);
   ch1->SetBranchAddress("RecoLepEnNumu", &RecoLepEnNumu);

   ch1->SetBranchStatus("RecoLepEnNue", true);
   ch1->SetBranchAddress("RecoLepEnNue", &RecoLepEnNue);

   ch1->SetBranchStatus("vtx_x", true);
   ch1->SetBranchAddress("vtx_x", &vtx_x);
   ch1->SetBranchStatus("vtx_y", true);
   ch1->SetBranchAddress("vtx_y", &vtx_y);
   ch1->SetBranchStatus("vtx_z", true);
   ch1->SetBranchAddress("vtx_z", &vtx_z);


   int nentries = ch1->GetEntries();



   /*
    Purity: True signal that passes the cut / Events that pass the cut
    Efficiency: True signal that passes the cut/True signal
    Sensitivity: True signal that passes the cut /sqrt(Events the pass the cut)
   */
   double efficiency = 0.;
   double purity = 0.;
   double sensitivity = 0.;

   double totalNutauSignal = 0.; //CC Nutaus that come from tauswap
   double cutCounter = 0.; //Anything that passes the cut
   double selectedNutauSignal = 0.; //This is for true signal that passes the cut
   double bkgdCounter = 0.;

   double totalPOTfromFile = 0.;
   double *binVariable = NULL;
   if (!useEHadBinning) {
      binVariable = &Ev;
   }
   else {
      binVariable = &eDepHad;
   }

   TString filename = "";

   /*
   Unfortunately, I haven't configured a way to quickly change which cuts to test.
   Right now, this macro is configured to analyze nutau selection with CVN cuts
   To change this, change:
      the PDG and filename in line 344
      the filename in the if statements starting on line 352
      the PDG cut in line 364

   After these changes, you must check the relevant labeling in the legend of the stacked histogram
   */

   for (auto i = 0; i < nentries; i++) {
   ch1->GetEntry(i);
   filename = ch1->GetCurrentFile()->GetName();
   eDepHad = eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther;
   DepPimRatio = eDepPim / eDepHad;

   if (filename == nonswap) {
      totalPOTfromFile = nonswapPOT;
   }
   else if (filename == nueswap) {
      totalPOTfromFile = nueswapPOT;
   }
   else {
      totalPOTfromFile = tauswapPOT;
   }

     if (nuPDG == 16 && isCC == 1 && filename == tauswap && cvnnue > 0 && cvnnumu > 0 && cvnnutau > 0 && cvnnc > 0 &&
         (vtx_x > -310) && (vtx_x < 310) && (vtx_y > -550) && (vtx_y < 550) && (vtx_z > 50) && (vtx_z < 1244)) {
         totalNutauSignal += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
         //hTrueSignal->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
     }
     if (cvnnue > 0 && cvnnumu > 0 && cvnnutau > 0 && cvnnc > 0 && cvnnue < 0.85 && cvnnumu < 0.5 &&
         (vtx_x > -310) && (vtx_x < 310) && (vtx_y > -550) && (vtx_y < 550) && (vtx_z > 50) && (vtx_z < 1244)) {
        if (isCC) {
            if (filename == nueswap) {
               cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
               //hCutCounter->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
            }
            else if (filename == nonswap){
               cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
               //hCutCounter->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
            }
            else {
               cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
               //hCutCounter->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
               if (nuPDG == 16) {
                  selectedNutauSignal += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                  //hTrueSignalCut->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
               }
            }
        }
        else {
           if (filename == nonswap) {
             cutCounter += (scalefactor / totalPOTfromFile);
             //hCutCounter->Fill(*binVariable, (scalefactor / totalPOTfromFile));
          }
        }
     }
     else {
        if (isCC) {
           bkgdCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
        }
        else {
           if (filename == nonswap)
            bkgdCounter += (scalefactor / totalPOTfromFile);
        }
     }
   if (i % 10000 == 0)
    std::cerr << 100 * (float)(i + 1) / (float) nentries << "%" << '\n';
   }

   //Cutcounter is S + B
   efficiency = selectedNutauSignal / totalNutauSignal;
   purity = selectedNutauSignal / cutCounter;
   sensitivity = selectedNutauSignal / Sqrt(cutCounter - selectedNutauSignal); //S / Sqrt(B)

   /*hEfficiency = (* hTrueSignalCut) / (* hTrueSignal);
   hPurity = (* hTrueSignalCut) / (* hCutCounter);*/


   std::cout << "selectedNutauSignal = " << selectedNutauSignal << '\n';
   std::cout << "totalNutauSignal = " << totalNutauSignal << '\n';
   std::cout << "cutCounter = " << cutCounter << '\n';
   std::cout << "bkgdCounter = " << bkgdCounter << std::endl;


   std::cout << "Purity = " << purity << '\n';
   std::cout << "Efficiency = " << efficiency << '\n';
   std::cout << "S / Sqrt(B) = " << sensitivity << '\n';
   std::cout << "S / Sqrt(S + B) = " << selectedNutauSignal / Sqrt(cutCounter) << std::endl;

   c1->cd();
   /*
   hRecoEVis->SetLineWidth(2);
   hRecoEVis->SetTitle("Reconstruted Energy FHC");
   hRecoEVis->GetXaxis()->SetTitle("Reco E_vis (GeV)");
   hRecoEVis->Draw("HIST");
   c1->SaveAs("./Figures/RecoEnergyhisto.png");
   c1->Clear();
   */

   /*hEfficiency.SetLineWidth(2);
   hEfficiency.SetTitle("FHC #nu_{#tau} signal");
   hEfficiency.GetXaxis()->SetTitle("True E (GeV)");
   hEfficiency.GetYaxis()->SetTitle("Signal Efficiency");
   hEfficiency.Draw("HIST");
   c1->SaveAs("./Figures/NutauEfficiency.png");
   c1->Clear();
   c1->Update();

   hPurity.SetLineWidth(2);
   hPurity.SetTitle("FHC #nu_{#tau} signal");
   hPurity.GetXaxis()->SetTitle("True E (GeV)");
   hPurity.GetYaxis()->SetTitle("Signal Purity");
   hPurity.Draw("HIST");
   c1->SaveAs("./Figures/NutauPurity.png");
   c1->Clear();
   c1->Update();*/


}
