/*
This code takes the FD MC swap files and plots a variety of statistics

               FIX SAVING THE NORMALIZED SIGNAL/BKGD HISTOGRAMS ON LINE 404
*/


#include <iostream>
#include <cstdlib>
#include <string>

#include "TMVA/Types.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

TString FHCnonswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nonswap.root";
TString FHCnueswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nueswap.root";
TString FHCtauswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_tauswap.root";

TString RHCnonswap = "/storage/shared/wshi/CAFs/FD/FD_RHC_nonswap.root";
TString RHCnueswap = "/storage/shared/wshi/CAFs/FD/FD_RHC_nueswap.root";
TString RHCtauswap = "/storage/shared/wshi/CAFs/FD/FD_RHC_tauswap.root";

using namespace TMath;

const double DELMSQ_31 = 2.515e-3; //In eV^2
const double LOSC = 1300.; //In km

const double THETA_23 = 0.859;
const double THETA_13 = 0.150;

//Scaling
double scaleCorrection = 40. / 1.13;
double POTperYear = 1.1e21;
double FHCnonswapPOT = 1.62824e24;
double FHCnueswapPOT = 1.64546e24;
double FHCtauswapPOT = 5.18551e24;

double RHCnonswapPOT = 3.27608e+24;
double RHCnueswapPOT = 3.24713e+24;
double RHCtauswapPOT = 8.58955e+24;



/*
  Units: deltam^2 in eV^2, E in GeV, L in km
    This gives sin^2(deltam^2 L / 4E) --> sin^2(pi L(km) deltam^2 (eV^2) /2.48 E(GeV)) [Equation 15 in https://arxiv.org/pdf/1710.00715.pdf]
  numu disappearance: P = 1 - sin^2(2 theta_23)* sin^2(pi L deltam^2 / 2.48 E)
  numu->nue: P = sin^2(theta_23)sin^2(2 theta_13)*sin^2(pi L deltam^2 / 2.48 E)
  numu->nutau: P = cos^4(theta_13)sin^2(2 theta_23)*sin^2(pi L deltam^2 / 2.48 E)

  nue disappearance: P = 1 - sin^2(2 theta_13)*sin^2(pi L deltam^2 / 2.48 E)
  nue->numu: P = sin^2(2 theta_13)*sin^2(theta_23)*sin^2(...)
  nue->nutau P = sin^2(2 theta_13)*Cos^2(theta_23)*sin^2(...)
*/
double OscWeight(const TString filename, const double Ev, const int nuPDG, const double delmsq, const double theta_23, const double theta_13) {
  //First check flavor, then check which file it came from. Apply appropriate oscillation probability
  if (nuPDG == 12 || nuPDG == -12) {
    if (filename == FHCnonswap || filename == RHCnonswap)
      return 1. - Power(Sin(2 * theta_13) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply nue disappearance
    else if (filename == FHCnueswap || filename == RHCnueswap)
      return Power(Sin(theta_23) * Sin(2 * theta_13) * Sin((Pi() * delmsq * LOSC) / (2.48 * Ev)), 2.); //Apply numu->nue
    else if (filename == FHCtauswap || filename == RHCtauswap) {
      std::cerr << "This should not happen: (nu_e in tauswap) \n";
      return 0.;
    } else {
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

void CafCutStatPlot(void) {

   bool useRHC = false;

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



   TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
   c1->SetRightMargin(0.16);

   TChain* ch1 = new TChain("cafTree");
   ch1->Add(nonswap);
   ch1->Add(nueswap);
   ch1->Add(tauswap);


   TH1D* cvnnue_fromnue = new TH1D("cvnnue_fromnue", "cvnnue_fromnue", 50, 0., 1.);
   TH1D* cvnnue_fromnumu = new TH1D("cvnnue_fromnumu", "cvnnue_fromnumu", 50, 0., 1.);
   TH1D* cvnnue_fromnutau = new TH1D("cvnnue_fromnutau", "cvnnue_fromnutau", 50, 0., 1.);
   TH1D* cvnnue_fromnc = new TH1D("cvnnue_fromnc", "cvnnue_fromnc", 50, 0., 1.);

   TH1D* cvnnumu_fromnue = new TH1D("cvnnumu_fromnue", "cvnnumu_fromnue", 50, 0., 1.);
   TH1D* cvnnumu_fromnumu = new TH1D("cvnnumu_fromnumu", "cvnnumu_fromnumu", 50, 0., 1.);
   TH1D* cvnnumu_fromnutau = new TH1D("cvnnumu_fromnutau", "cvnnumu_fromnutau", 50, 0., 1.);
   TH1D* cvnnumu_fromnc = new TH1D("cvnnumu_fromnc", "cvnnumu_fromnc", 50, 0., 1.);

   TH1D* cvnnutau_fromnue = new TH1D("cvnnutau_fromnue", "cvnnutau_fromnue", 50, 0., 1.);
   TH1D* cvnnutau_fromnumu = new TH1D("cvnnutau_fromnumu", "cvnnutau_fromnumu", 50, 0., 1.);
   TH1D* cvnnutau_fromnutau = new TH1D("cvnnutau_fromnutau", "cvnnutau_fromnutau", 50, 0., 1.);
   TH1D* cvnnutau_fromnc = new TH1D("cvnnutau_fromnc", "cvnnutau_fromnc", 50, 0., 1.);

   TH1D* cvnnc_fromnue = new TH1D("cvnnc_fromnue", "cvnnc_fromnue", 50, 0., 1.);
   TH1D* cvnnc_fromnumu = new TH1D("cvnnc_fromnumu", "cvnnc_fromnumu", 50, 0., 1.);
   TH1D* cvnnc_fromnutau = new TH1D("cvnnc_fromnutau", "cvnnc_fromnutau", 50, 0., 1.);
   TH1D* cvnnc_fromnc = new TH1D("cvnnc_fromnc", "cvnnc_fromnc", 50, 0., 1.);



   TH1D* hnipimSignal = new TH1D("hnipimSignal", "hnipimSignal", 5, 0., 5);
   TH1D* hnipimBkgd = new TH1D("hnipimBkgd", "hnipimBkgd", 5, 0., 5);

   TH1D* hcvnnutauSignal = new TH1D("hcvnnutauSignal", "hcvnnutauSignal", 50, 0., 1.);
   TH1D* hcvnnutauBkgd = new TH1D("hcvnnutauBkgd", "hcvnnutauBkgd", 50, 0., 1.);

   TH1D* heDepHadSignal = new TH1D("heDepHadSignal", "heDepHadSignal", 25, 0., 10.);
   TH1D* heDepHadBkgd = new TH1D("heDepHadBkgd", "heDepHadBkgd", 25, 0., 10.);


   /*
   Energy resolution stuff that I made a while ago. These could easily be reimplemented into this analysis

   TH2D* hResolutionE = new TH2D("hResolutionE", "hResolutionE", 10, 0., 1., 25, 0., 10.);
   TH2D* hResolutionM = new TH2D("hResolutionM", "hResolutionM", 10, 0., 1., 25, 0., 10.);
   TH2D* hResolutionT = new TH2D("hResolutionT", "hResolutionT", 10, 0., 1., 25, 0., 10.);

   TH2D* hTruevRecoEE = new TH2D("hTruevRecoEE", "hTruevRecoEE", 25, 0., 10., 25, 0., 10.);
   TH2D* hTruevRecoEM = new TH2D("hTruevRecoEM", "hTruevRecoEM", 25, 0., 10., 25, 0., 10.);
   TH2D* hTruevRecoET = new TH2D("hTruevRecoET", "hTruevRecoET", 25, 0., 10., 25, 0., 10.);

   TH2D* hResolutionE2 = new TH2D("hResolutionE2", "hResolutionE2", 10, 0., 1., 25, 0., 10.);
   TH2D* hResolutionM2 = new TH2D("hResolutionM2", "hResolutionM2", 10, 0., 1., 25, 0., 10.);
   TH2D* hResolutionT2 = new TH2D("hResolutionT2", "hResolutionT2", 10, 0., 1., 25, 0., 10.);

   TH1D* hResolutionE1D = new TH1D("hResolutionE1D", "hResolutionE1D", 10, 0., 1.);
   TH1D* hResolutionM1D = new TH1D("hResolutionM1D", "hResolutionM1D", 10, 0., 1.);
   TH1D* hResolutionT1D = new TH1D("hResolutionT1D", "hResolutionT1D", 10, 0., 1.);
   */

   gStyle->SetOptStat(0);


   double Ev = 0.;
   int nuPDG = 0.;
   int isCC = -1.;
   double cvnnue = 0.;
   double cvnnumu = 0.;
   double cvnnutau = 0.;
   double cvnnc = 0.;

   double eDepP = 0.;
   double eDepN = 0.;
   double eDepPip = 0.;
   double eDepPim = 0.;
   double eDepPi0 = 0.;
   double eDepOther = 0.;
   double eDepTotalE = 0.;
   double eDepTotalM = 0.;
   double DepPimRatio = 0.;
   double RecoLepEnNumu = 0.;
   double RecoLepEnNue = 0.;
   double eDepHad = 0.;
   double eDepLep = 0.;

   double eRecoP = 0.;
   double eRecoN = 0.;
   double eRecoPip = 0.;
   double eRecoPim = 0.;
   double eRecoPi0 = 0.;
   double eRecoOther = 0.;
   double eRecoTotal = 0.;

   double vtx_x = 0.;
   double vtx_y = 0.;
   double vtx_z = 0.;

   int nipim = 0.;
   int nipip = 0.;

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

   ch1->SetBranchStatus("eRecoP", true);
   ch1->SetBranchAddress("eRecoP", &eRecoP);

   ch1->SetBranchStatus("eRecoN", true);
   ch1->SetBranchAddress("eRecoN", &eRecoN);

   ch1->SetBranchStatus("eRecoPip", true);
   ch1->SetBranchAddress("eRecoPip", &eRecoPip);

   ch1->SetBranchStatus("eRecoPim", true);
   ch1->SetBranchAddress("eRecoPim", &eRecoPim);

   ch1->SetBranchStatus("eRecoPi0", true);
   ch1->SetBranchAddress("eRecoPi0", &eRecoPi0);

   ch1->SetBranchStatus("eRecoOther", true);
   ch1->SetBranchAddress("eRecoOther", &eRecoOther);

   ch1->SetBranchStatus("RecoLepEnNumu", true);
   ch1->SetBranchAddress("RecoLepEnNumu", &RecoLepEnNumu);

   ch1->SetBranchStatus("RecoLepEnNue", true);
   ch1->SetBranchAddress("RecoLepEnNue", &RecoLepEnNue);

   ch1->SetBranchStatus("cvnnc", true);
   ch1->SetBranchAddress("cvnnc", &cvnnc);

   ch1->SetBranchStatus("cvnnue", true);
   ch1->SetBranchAddress("cvnnue", &cvnnue);

   ch1->SetBranchStatus("cvnnumu", true);
   ch1->SetBranchAddress("cvnnumu", &cvnnumu);

   ch1->SetBranchStatus("cvnnutau", true);
   ch1->SetBranchAddress("cvnnutau", &cvnnutau);

   ch1->SetBranchStatus("vtx_x", true);
   ch1->SetBranchAddress("vtx_x", &vtx_x);

   ch1->SetBranchStatus("vtx_y", true);
   ch1->SetBranchAddress("vtx_y", &vtx_y);

   ch1->SetBranchStatus("vtx_z", true);
   ch1->SetBranchAddress("vtx_z", &vtx_z);

   ch1->SetBranchStatus("nipim", true);
   ch1->SetBranchAddress("nipim", &nipim);

   ch1->SetBranchStatus("nipip", true);
   ch1->SetBranchAddress("nipip", &nipip);



   TString filename = "";

   int nentries = ch1->GetEntries();

   double totalPOTfromFile = 0.;

   for (auto i = 0; i < nentries; i++) {
      ch1->GetEntry(i);
      filename = ch1->GetCurrentFile()->GetName();
      eDepHad = eDepPip + eDepPi0 + eDepPim + eDepP + eDepN + eDepOther;

      if (filename == nonswap) {
         totalPOTfromFile = nonswapPOT;
      }
      else if (filename == nueswap) {
         totalPOTfromFile = nueswapPOT;
      }
      else {
         totalPOTfromFile = tauswapPOT;
      }

      if (cvnnue > 0 && cvnnumu > 0 && cvnnutau > 0 && cvnnc > 0 && (vtx_x > -310) &&
         (vtx_x < 310) && (vtx_y > -550) && (vtx_y < 550) && (vtx_z > 50) && (vtx_z < 1244)) {
         if (isCC) {
            if (nuPDG == 14 || nuPDG == -14) {
               cvnnue_fromnumu->Fill(cvnnue, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnumu_fromnumu->Fill(cvnnumu, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnutau_fromnumu->Fill(cvnnutau, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnc_fromnumu->Fill(cvnnc, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
            }
            else if (nuPDG == 12 || nuPDG == -12) {
               cvnnue_fromnue->Fill(cvnnue, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnumu_fromnue->Fill(cvnnumu, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnutau_fromnue->Fill(cvnnutau, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnc_fromnue->Fill(cvnnc, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
            }
            else {
               cvnnue_fromnutau->Fill(cvnnue, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnumu_fromnutau->Fill(cvnnumu, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnutau_fromnutau->Fill(cvnnutau, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               cvnnc_fromnutau->Fill(cvnnc, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
            }

            if (nuPDG == 16 && filename == tauswap) {
               heDepHadSignal->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               hcvnnutauSignal->Fill(cvnnutau, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               hnipimSignal->Fill(nipim, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
            }
            else {
               heDepHadBkgd->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               hcvnnutauBkgd->Fill(cvnnutau, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
               hnipimBkgd->Fill(nipim, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * scalefactor / totalPOTfromFile);
            }
         }
         else {
            if (filename == nonswap) {
               heDepHadBkgd->Fill(eDepHad, scalefactor / totalPOTfromFile);
               hcvnnutauBkgd->Fill(cvnnutau, scalefactor / totalPOTfromFile);
               hnipimBkgd->Fill(nipim, scalefactor / totalPOTfromFile);

               cvnnue_fromnc->Fill(cvnnue, scalefactor / totalPOTfromFile);
               cvnnumu_fromnc->Fill(cvnnumu, scalefactor / totalPOTfromFile);
               cvnnutau_fromnc->Fill(cvnnutau, scalefactor / totalPOTfromFile);
               cvnnc_fromnc->Fill(cvnnc, scalefactor / totalPOTfromFile);
            }
         }
      }

      if (i % 10000 == 0) std::cout << (float) (i + 1) * 100 / (float) nentries << '%' << '\n';
   }

   c1->cd();
   TH1D* heDepHadSignalNorm = (TH1D*) (heDepHadSignal->Clone("heDepHadSignal"));
   heDepHadSignalNorm->Scale(1. / heDepHadSignalNorm->Integral(), "width");

   TH1D* heDepHadBkgdNorm = (TH1D*) (heDepHadBkgd->Clone("heDepHadBkgd"));
   heDepHadBkgdNorm->Scale(1. / heDepHadBkgdNorm->Integral(), "width");

   TLegend* heDepHadleg = new TLegend(0.64, 0.8, 0.84, 0.9);
   heDepHadleg->AddEntry(heDepHadSignalNorm, "Signal");
   heDepHadleg->AddEntry(heDepHadBkgdNorm, "Background");

   heDepHadSignalNorm->SetFillColorAlpha(kBlue, 0.5);
   heDepHadSignalNorm->SetLineColor(kBlue + 1);
   heDepHadSignalNorm->SetLineWidth(2);
   heDepHadSignalNorm->SetFillStyle(3001);
   heDepHadSignalNorm->GetXaxis()->SetTitle("Reco E_{Had} (GeV)");
   heDepHadSignalNorm->GetYaxis()->SetTitle("1/N dN / dx");
   heDepHadSignalNorm->SetTitle("");


   heDepHadBkgdNorm->SetFillColorAlpha(kRed, 1.);
   heDepHadBkgdNorm->SetLineColor(kRed);
   heDepHadBkgdNorm->SetLineWidth(2);
   gStyle->SetHatchesLineWidth(2);
   heDepHadBkgdNorm->SetFillStyle(3554);

   heDepHadSignalNorm->Draw("HIST");
   heDepHadBkgdNorm->Draw("SAMEHIST");
   heDepHadleg->Draw();

   c1->SaveAs("./Figures/eDepHadDistrib.png");
   c1->Clear();
   c1->Update();

   TH1D* hcvnnutauSignalNorm = (TH1D*) (hcvnnutauSignal->Clone("hcvnnutauSignal"));
   hcvnnutauSignalNorm->Scale(1. / hcvnnutauSignalNorm->Integral(), "width");

   TH1D* hcvnnutauBkgdNorm = (TH1D*) (hcvnnutauBkgd->Clone("hcvnnutauBkgd"));
   hcvnnutauBkgdNorm->Scale(1. / hcvnnutauBkgdNorm->Integral(), "width");

   TLegend* hcvnnutauleg = new TLegend(0.64, 0.8, 0.84, 0.9);
   hcvnnutauleg->AddEntry(hcvnnutauSignalNorm, "Signal");
   hcvnnutauleg->AddEntry(hcvnnutauBkgdNorm, "Background");

   hcvnnutauSignalNorm->SetFillColorAlpha(kBlue, 0.5);
   hcvnnutauSignalNorm->SetLineColor(kBlue + 1);
   hcvnnutauSignalNorm->SetLineWidth(2);
   hcvnnutauSignalNorm->SetFillStyle(3001);
   hcvnnutauSignalNorm->GetXaxis()->SetTitle("CVN #nu_{#tau} score");
   hcvnnutauSignalNorm->GetYaxis()->SetTitle("1/N dN / dx");
   hcvnnutauSignalNorm->SetTitle("");


   hcvnnutauBkgdNorm->SetFillColorAlpha(kRed, 1.);
   hcvnnutauBkgdNorm->SetLineColor(kRed);
   hcvnnutauBkgdNorm->SetLineWidth(2);
   gStyle->SetHatchesLineWidth(2);
   hcvnnutauBkgdNorm->SetFillStyle(3554);

   hcvnnutauSignalNorm->Draw("HIST");
   hcvnnutauBkgdNorm->Draw("SAMEHIST");
   hcvnnutauleg->Draw();

   c1->SaveAs("./Figures/CVNNutauDistrib.png");
   c1->Clear();
   c1->Update();

   TH1D* hnipimSignalNorm = (TH1D*) (hnipimSignal->Clone("hnipimSignal"));
   hnipimSignalNorm->Scale(1. / hnipimSignalNorm->Integral(), "width");

   TH1D* hnipimBkgdNorm = (TH1D*) (hnipimBkgd->Clone("hnipimBkgd"));
   hnipimBkgdNorm->Scale(1. / hnipimBkgdNorm->Integral(), "width");

   TLegend* hnipimleg = new TLegend(0.64, 0.8, 0.84, 0.9);
   hnipimleg->AddEntry(hnipimSignalNorm, "Signal");
   hnipimleg->AddEntry(hnipimBkgdNorm, "Background");

   hnipimSignalNorm->SetFillColorAlpha(kBlue, 0.5);
   hnipimSignalNorm->SetLineColor(kBlue + 1);
   hnipimSignalNorm->SetLineWidth(2);
   hnipimSignalNorm->SetFillStyle(3001);
   hnipimSignalNorm->GetXaxis()->SetTitle("Reco E_{Had} (GeV)");
   hnipimSignalNorm->GetYaxis()->SetTitle("1/N dN / dx");
   hnipimSignalNorm->SetTitle("");


   hnipimBkgdNorm->SetFillColorAlpha(kRed, 1.);
   hnipimBkgdNorm->SetLineColor(kRed);
   hnipimBkgdNorm->SetLineWidth(2);
   gStyle->SetHatchesLineWidth(2);
   hnipimBkgdNorm->SetFillStyle(3554);

   hnipimSignalNorm->Draw("HIST");
   hnipimBkgdNorm->Draw("SAMEHIST");
   hnipimleg->Draw();

   c1->SaveAs("./Figures/NPimDistrib.png");
   c1->Clear();
   c1->Update();

   /*outputFile->WriteObject(wcvnnue_fromnue, "wcvnnue_fromnue");
   outputFile->WriteObject(wcvnnue_fromnumu, "wcvnnue_fromnumu");
   outputFile->WriteObject(wcvnnue_fromnutau, "wcvnnue_fromnutau");
   outputFile->WriteObject(wcvnnue_fromnc, "wcvnnue_fromnc");

   outputFile->WriteObject(wcvnnumu_fromnue, "wcvnnumu_fromnue");
   outputFile->WriteObject(wcvnnumu_fromnumu, "wcvnnumu_fromnumu");
   outputFile->WriteObject(wcvnnumu_fromnutau, "wcvnnumu_fromnutau");
   outputFile->WriteObject(wcvnnumu_fromnc, "wcvnnumu_fromnc");

   outputFile->WriteObject(wcvnnutau_fromnue, "wcvnnutau_fromnue");
   outputFile->WriteObject(wcvnnutau_fromnumu, "wcvnnutau_fromnumu");
   outputFile->WriteObject(wcvnnutau_fromnutau, "wcvnnutau_fromnutau");
   outputFile->WriteObject(wcvnnutau_fromnc, "wcvnnutau_fromnc");

   outputFile->WriteObject(wcvnnc_fromnue, "wcvnnc_fromnue");
   outputFile->WriteObject(wcvnnc_fromnumu, "wcvnnc_fromnumu");
   outputFile->WriteObject(wcvnnc_fromnutau, "wcvnnc_fromnutau");
   outputFile->WriteObject(wcvnnc_fromnc, "wcvnnc_fromnc");

   outputFile->WriteObject(nuetonue, "nuetonuehisto");
   outputFile->WriteObject(nuetonumu, "nuetonumuhisto");
   outputFile->WriteObject(nuetonutau, "nuetonutauhisto");

   outputFile->WriteObject(numutonue, "numutonuehisto");
   outputFile->WriteObject(numutonumu, "numutonumuhisto");
   outputFile->WriteObject(numutonutau, "numutonutauhisto");*/



   cvnnue_fromnue->SetLineWidth(2);
   cvnnue_fromnue->SetTitle("FHC CVN nu_e scores from nu_e events");
   cvnnue_fromnue->GetXaxis()->SetTitle("CVN score");
   cvnnue_fromnue->Draw("HIST");
   c1->SaveAs("./Figures/cvnnue_fromnue.png");
   c1->Clear();
   c1->Update();

   cvnnue_fromnumu->SetLineWidth(2);
   cvnnue_fromnumu->SetTitle("FHC CVN nu_e scores from nu_mu events");
   cvnnue_fromnumu->GetXaxis()->SetTitle("CVN score");
   cvnnue_fromnumu->Draw("HIST");
   c1->SaveAs("./Figures/cvnnue_fromnumu.png");
   c1->Clear();
   c1->Update();

   cvnnue_fromnutau->SetLineWidth(2);
   cvnnue_fromnutau->SetTitle("FHC CVN nu_e scores from nu_tau events");
   cvnnue_fromnutau->GetXaxis()->SetTitle("CVN score");
   cvnnue_fromnutau->Draw("HIST");
   c1->SaveAs("./Figures/cvnnue_fromnutau.png");
   c1->Clear();
   c1->Update();

   cvnnue_fromnc->SetLineWidth(2);
   cvnnue_fromnc->SetTitle("FHC CVN nu_e scores from Neutral Current events");
   cvnnue_fromnc->GetXaxis()->SetTitle("CVN score");
   cvnnue_fromnc->Draw("HIST");
   c1->SaveAs("./Figures/cvnnue_fromnc.png");
   c1->Clear();
   c1->Update();

   cvnnumu_fromnue->SetLineWidth(2);
   cvnnumu_fromnue->SetTitle("FHC CVN nu_mu scores from nu_e events");
   cvnnumu_fromnue->GetXaxis()->SetTitle("CVN score");
   cvnnumu_fromnue->Draw("HIST");
   c1->SaveAs("./Figures/cvnnumu_fromnue.png");
   c1->Clear();
   c1->Update();

   cvnnumu_fromnumu->SetLineWidth(2);
   cvnnumu_fromnumu->SetTitle("FHC CVN nu_mu scores from nu_mu events");
   cvnnumu_fromnumu->GetXaxis()->SetTitle("CVN score");
   cvnnumu_fromnumu->Draw("HIST");
   c1->SaveAs("./Figures/cvnnumu_fromnumu.png");
   c1->Clear();
   c1->Update();

   cvnnumu_fromnutau->SetLineWidth(2);
   cvnnumu_fromnutau->SetTitle("FHC CVN nu_mu scores from nu_tau events");
   cvnnumu_fromnutau->GetXaxis()->SetTitle("CVN score");
   cvnnumu_fromnutau->Draw("HIST");
   c1->SaveAs("./Figures/cvnnumu_fromnutau.png");
   c1->Clear();
   c1->Update();

   cvnnumu_fromnc->SetLineWidth(2);
   cvnnumu_fromnc->SetTitle("FHC CVN nu_mu scores from Neutral Current events");
   cvnnumu_fromnc->GetXaxis()->SetTitle("CVN score");
   cvnnumu_fromnc->Draw("HIST");
   c1->SaveAs("./Figures/cvnnumu_fromnc.png");
   c1->Clear();
   c1->Update();

   cvnnutau_fromnue->SetLineWidth(2);
   cvnnutau_fromnue->SetTitle("FHC CVN nu_tau scores from nu_e events");
   cvnnutau_fromnue->GetXaxis()->SetTitle("CVN score");
   cvnnutau_fromnue->Draw("HIST");
   c1->SaveAs("./Figures/cvnnutau_fromnue.png");
   c1->Clear();
   c1->Update();

   cvnnutau_fromnumu->SetLineWidth(2);
   cvnnutau_fromnumu->SetTitle("FHC CVN nu_tau scores from nu_mu events");
   cvnnutau_fromnumu->GetXaxis()->SetTitle("CVN score");
   cvnnutau_fromnumu->Draw("HIST");
   c1->SaveAs("./Figures/cvnnutau_fromnumu.png");
   c1->Clear();
   c1->Update();

   cvnnutau_fromnutau->SetLineWidth(2);
   cvnnutau_fromnutau->SetTitle("FHC CVN nu_tau scores from nu_tau events");
   cvnnutau_fromnutau->GetXaxis()->SetTitle("CVN score");
   cvnnutau_fromnutau->Draw("HIST");
   c1->SaveAs("./Figures/cvnnutau_fromnutau.png");
   c1->Clear();
   c1->Update();

   cvnnutau_fromnc->SetLineWidth(2);
   cvnnutau_fromnc->SetTitle("FHC CVN nu_tau scores from Neutral Current events");
   cvnnutau_fromnc->GetXaxis()->SetTitle("CVN score");
   cvnnutau_fromnc->Draw("HIST");
   c1->SaveAs("./Figures/cvnnutau_fromnc.png");
   c1->Clear();
   c1->Update();

   cvnnc_fromnue->SetLineWidth(2);
   cvnnc_fromnue->SetTitle("FHC CVN Neutral Current scores from nu_e events");
   cvnnc_fromnue->GetXaxis()->SetTitle("CVN score");
   cvnnc_fromnue->Draw("HIST");
   c1->SaveAs("./Figures/cvnnc_fromnue.png");
   c1->Clear();
   c1->Update();

   cvnnc_fromnumu->SetLineWidth(2);
   cvnnc_fromnumu->SetTitle("FHC CVN Neutral Current scores from nu_mu events");
   cvnnc_fromnumu->GetXaxis()->SetTitle("CVN score");
   cvnnc_fromnumu->Draw("HIST");
   c1->SaveAs("./Figures/cvnnc_fromnumu.png");
   c1->Clear();
   c1->Update();

   cvnnc_fromnutau->SetLineWidth(2);
   cvnnc_fromnutau->SetTitle("FHC CVN Neutral Current scores from nu_tau events");
   cvnnc_fromnutau->GetXaxis()->SetTitle("CVN score");
   cvnnc_fromnutau->Draw("HIST");
   c1->SaveAs("./Figures/cvnnc_fromnutau.png");
   c1->Clear();
   c1->Update();

   cvnnc_fromnc->SetLineWidth(2);
   cvnnc_fromnc->SetTitle("FHC CVN Neutral Current scores from Neutral Current events");
   cvnnc_fromnc->GetXaxis()->SetTitle("CVN score");
   cvnnc_fromnc->Draw("HIST");
   c1->SaveAs("./Figures/cvnnc_fromnc.png");
   c1->Clear();
   c1->Update();

   /*
   hResolutionE->Draw("COLZ");
   hResolutionE->SetTitle("FHC Fractional Energy Resolution for #nu_{e}");
   hResolutionE->GetXaxis()->SetTitle("(True E - Reco. E) / True E");
   hResolutionE->GetYaxis()->SetTitle("True E (GeV)");
   c1->SaveAs("./Figures/ResolutionE.png");
   c1->Clear();
   c1->Update();

   hResolutionM->Draw("COLZ");
   hResolutionM->SetTitle("FHC Fractional Energy Resolution for #nu_{#mu}");
   hResolutionM->GetXaxis()->SetTitle("(True E - Reco. E) / True E");
   hResolutionM->GetYaxis()->SetTitle("True E (GeV)");
   c1->SaveAs("./Figures/ResolutionM.png");
   c1->Clear();
   c1->Update();

   /*hResolutionT->Draw("COLZ");
   hResolutionT->SetTitle("FHC Fractional Energy Resolution for #nu_{#tau}");
   hResolutionT->GetXaxis()->SetTitle("(True E - Reco. E) / True E");
   hResolutionT->GetYaxis()->SetTitle("True E (GeV)");
   c1->SaveAs("./Figures/ResolutionT.png");
   c1->Clear();
   c1->Update();


   hResolutionE2->Draw("COLZ");
   hResolutionE2->SetTitle("FHC Energy Resolution for #nu_{e}");
   hResolutionE2->GetXaxis()->SetTitle("Reco. E / True E");
   hResolutionE2->GetYaxis()->SetTitle("True E (GeV)");
   c1->SaveAs("./Figures/ResolutionE2.png");
   c1->Clear();
   c1->Update();

   hResolutionM2->Draw("COLZ");
   hResolutionM2->SetTitle("FHC Energy Resolution for #nu_{#mu}");
   hResolutionM2->GetXaxis()->SetTitle("Reco. E / True E");
   hResolutionM2->GetYaxis()->SetTitle("True E (GeV)");
   c1->SaveAs("./Figures/ResolutionM2.png");
   c1->Clear();
   c1->Update();

   /*hResolutionT2->Draw("COLZ");
   hResolutionT2->SetTitle("FHC Energy Resolution for #nu_{#tau}");
   hResolutionT2->GetXaxis()->SetTitle("Reco. E / True E");
   hResolutionT2->GetYaxis()->SetTitle("True E (GeV)");
   c1->SaveAs("./Figures/ResolutionT2.png");
   c1->Clear();
   c1->Update();

   hResolutionE1D->SetLineWidth(2);
   hResolutionE1D->SetTitle("FHC Energy Resolution for #nu_{e} (2 GeV < Ev < 3 GeV)");
   hResolutionE1D->GetXaxis()->SetTitle("Reco. E / True E");
   hResolutionE1D->Draw("HIST");
   c1->SaveAs("./Figures/ResolutionE1D.png");
   c1->Clear();
   c1->Update();

   hResolutionM1D->SetLineWidth(2);
   hResolutionM1D->SetTitle("FHC Energy Resolution for #nu_{#mu} (2 GeV < Ev < 3 GeV)");
   hResolutionM1D->GetXaxis()->SetTitle("Reco. E / True E");
   hResolutionM1D->Draw("HIST");
   c1->SaveAs("./Figures/ResolutionM1D.png");
   c1->Clear();
   c1->Update();

   /*hResolutionT1D->SetLineWidth(2);
   hResolutionT1D->SetTitle("FHC Energy Resolution for #nu_{#tau} (3.5 GeV < Ev < 4.5 GeV)");
   hResolutionT1D->GetXaxis()->SetTitle("Reco. E / True E");
   hResolutionT1D->Draw("HIST");
   c1->SaveAs("./Figures/ResolutionT1D.png");
   c1->Clear();
   c1->Update();*/



   //outputFile->Close();

}
