/*
Used to find optimal cuts output by the BDT for implementation into CutTest.C
*/

#include <iostream>
#include <cstdlib>
#include <string>

#include "TMVA/Types.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

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


void MVASignalEval(void) {

   //Uncomment to save graphs in a ROOT file
   //TFile* outputFile = new TFile("MVAEvaluatedMCData.root", "RECREATE")

   bool useRHC = false;

   double nonswapPOT = 0.;
   double nueswapPOT = 0.;
   double tauswapPOT = 0.;

   TString nonswap = "";
   TString nueswap = "";
   TString tauswap = "";

   //Change this whenever any change is made to the training sample/training procedure
   double optimalBDTCut = 0.;

   if (!useRHC) {
      nonswapPOT = FHCnonswapPOT;
      nueswapPOT = FHCnueswapPOT;
      tauswapPOT = FHCtauswapPOT;

      nonswap = FHCnonswap;
      nueswap = FHCnueswap;
      tauswap = FHCtauswap;

      optimalBDTCut = 0.1206;
   }
   else {
      nonswapPOT = RHCnonswapPOT;
      nueswapPOT = RHCnueswapPOT;
      tauswapPOT = RHCtauswapPOT;

      nonswap = RHCnonswap;
      nueswap = RHCnueswap;
      tauswap = RHCtauswap;

      optimalBDTCut = 0.1132;
   }

   double scalefactor = POTperYear * scaleCorrection * 3.5; //3.5 years of data taking

   gStyle->SetOptStat(0);
   TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
   c1->SetRightMargin(0.1);
   c1->SetLeftMargin(0.15);

   TChain* ch1 = new TChain("cafTree");
   ch1->Add(nonswap);
   ch1->Add(nueswap);
   ch1->Add(tauswap);

   THStack* hStack = new THStack("hStack", "");
   TLegend *leg = new TLegend(0.65,0.65, 0.9, 0.9);
   leg->SetBorderSize(0);

   THStack* hStackReco = new THStack("hStackReco", "");
   TLegend* legReco = new TLegend(0.65, 0.65, 0.9, 0.9);
   legReco->SetBorderSize(0);

   THStack* hStackRecoVis = new THStack("hStackRecoVis", "");
   TLegend* legRecoVis = new TLegend(0.65, 0.65, 0.9, 0.9);
   legRecoVis->SetBorderSize(0);

   TH1D* hSignal = new TH1D("NuTau Signal", "NuTau Signal", 25, 0., 10.);
   TH1D* hNCBkgd = new TH1D("NCBkgd", "NCBkgd", 25, 0., 10.);
   TH1D* hWLBkgdE = new TH1D("WLBkgdE", "WLBkgdE", 25, 0., 10.);
   TH1D* hWLBkgdM = new TH1D("WLBkgdM", "WLBkgdM", 25, 0., 10.);
   TH1D* hWSBkgdT = new TH1D("WSBkgdT", "WSBkgdT", 25, 0., 10.);
   TH1D* hIntBkgdT = new TH1D("WLBkgdT", "WLBkgdT", 25, 0., 10.);


   TH1D* hSignalReco = new TH1D("NuTau Signal Reco", "NuTau Signal Reco", 25, 0., 10.);
   TH1D* hNCBkgdReco = new TH1D("NCBkgd Reco", "NCBkgd Reco", 25, 0., 10.);
   TH1D* hWLBkgdEReco = new TH1D("WLBkgdE Reco", "WLBkgdE Reco", 25, 0., 10.);
   TH1D* hWLBkgdMReco = new TH1D("WLBkgdM Reco", "WLBkgdM Reco", 25, 0., 10.);
   TH1D* hWSBkgdTReco = new TH1D("WSBkgdT Reco", "WsBkgdT Reco", 25, 0., 10.);
   TH1D* hIntBkgdTReco = new TH1D("WLBkgdT Reco", "WLBkgdT Reco", 25, 0., 10.);

   TH1D* hSignalRecoVis = new TH1D("NuTau Signal Reco Vis", "NuTau Signal Reco Vis", 25, 0., 10.);
   TH1D* hNCBkgdRecoVis = new TH1D("NCBkgd Reco Vis", "NCBkgd Reco Vis", 25, 0., 10.);
   TH1D* hWLBkgdERecoVis = new TH1D("WLBkgdE Reco Vis", "WLBkgdE Reco Vis", 25, 0., 10.);
   TH1D* hWLBkgdMRecoVis = new TH1D("WLBkgdM Reco Vis", "WLBkgdM Reco Vis", 25, 0., 10.);
   TH1D* hWSBkgdTRecoVis = new TH1D("WSBkgdT Reco Vis", "WsBkgdT Reco Vis", 25, 0., 10.);
   TH1D* hIntBkgdTRecoVis = new TH1D("WLBkgdT Reco Vis", "WLBkgdT Reco Vis", 25, 0., 10.);

   TMVA::Reader* reader = new TMVA::Reader();

   float Ev = 0.;
   int nuPDG = 0.;
   int isCC = -1.;
   float cvnnue = 0.;
   float cvnnumu = 0.;
   float cvnnutau = 0.;
   float cvnnc = 0.;
   float EPimRatio = 0.;

   double eDepP = 0.;
   double eDepN = 0.;
   double eDepPip = 0.;
   float eDepPim = 0.;
   double eDepPi0 = 0.;
   double eDepOther = 0.;
   double eDepTotal = 0.;
   double DepPimRatio = 0.;
   double RecoLepEnNumu = 0.;
   double RecoLepEnNue = 0.;
   float eDepHad = 0.;
   double eDepLep = 0.;
   float nipim = 0.;
   float POTScaledOscweight = 0.;
   float nipip;


   double dcvnnue = 0.;
   double dcvnnumu = 0.;
   double dcvnnutau = 0.;
   double dcvnnc = 0.;
   double deDepPim = 0.;
   int inipim = 0;
   double dEv = 0.;
   int inipip = 0;
   double eDepHadVis = 0.;
   double RecoEVisE = 0.;
   double RecoEVisM = 0.;
   double RecoEVisT = 0.;

   double vtx_x = 0.;
   double vtx_y = 0.;
   double vtx_z = 0.;

   TString filename = "";

   ch1->SetBranchStatus("*", false);
   ch1->SetBranchStatus("isCC", true);
   ch1->SetBranchAddress("isCC", &isCC);
   ch1->SetBranchStatus("cvnnue", true);
   ch1->SetBranchAddress("cvnnue", &dcvnnue);
   ch1->SetBranchStatus("cvnnutau", true);
   ch1->SetBranchAddress("cvnnutau", &dcvnnutau);
   ch1->SetBranchStatus("cvnnumu", true);
   ch1->SetBranchAddress("cvnnumu", &dcvnnumu);
   ch1->SetBranchStatus("cvnnc", true);
   ch1->SetBranchAddress("cvnnc", &dcvnnc);
   ch1->SetBranchStatus("eDepPim", true);
   ch1->SetBranchAddress("eDepPim", &deDepPim);
   ch1->SetBranchStatus("eDepP", true);
   ch1->SetBranchAddress("eDepP", &eDepP);
   ch1->SetBranchStatus("eDepN", true);
   ch1->SetBranchAddress("eDepN", &eDepN);
   ch1->SetBranchStatus("eDepPip", true);
   ch1->SetBranchAddress("eDepPip", &eDepPip);
   ch1->SetBranchStatus("eDepPi0", true);
   ch1->SetBranchAddress("eDepPi0", &eDepPi0);
   ch1->SetBranchStatus("eDepOther", true);
   ch1->SetBranchAddress("eDepOther", &eDepOther);
   ch1->SetBranchStatus("nipim", true);
   ch1->SetBranchAddress("nipim", &inipim);
   ch1->SetBranchStatus("nuPDG", true);
   ch1->SetBranchAddress("nuPDG", &nuPDG);
   ch1->SetBranchStatus("Ev", true);
   ch1->SetBranchAddress("Ev", &dEv);
   ch1->SetBranchStatus("nipip", true);
   ch1->SetBranchAddress("nipip", &inipip);
   ch1->SetBranchStatus("vtx_x", true);
   ch1->SetBranchAddress("vtx_x", &vtx_x);
   ch1->SetBranchStatus("vtx_y", true);
   ch1->SetBranchAddress("vtx_y", &vtx_y);
   ch1->SetBranchStatus("vtx_z", true);
   ch1->SetBranchAddress("vtx_z", &vtx_z);
   ch1->SetBranchStatus("RecoLepEnNue", true);
   ch1->SetBranchAddress("RecoLepEnNue", &RecoLepEnNue);
   ch1->SetBranchStatus("RecoLepEnNumu", true);
   ch1->SetBranchAddress("RecoLepEnNumu", &RecoLepEnNumu);



   reader->AddVariable("cvnnue", &cvnnue);
   reader->AddVariable("cvnnumu", &cvnnumu);
   reader->AddVariable("cvnnutau", &cvnnutau);
   reader->AddVariable("EHad := eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther", &eDepHad);
   reader->AddVariable("eDepPim", &eDepPim);
   //My current version of the xml file doesn't have the eDepPip variable in the training. Uncomment this if you run your own training/use RHC.
   //reader->AddVariable("eDepPip", &eDepPip);
   reader->AddVariable("nipim", &nipim);
   reader->AddVariable("cvnnc", &cvnnc);
   reader->AddVariable("nipip", &nipip);

   reader->AddSpectator("Ev", &Ev);
   reader->AddSpectator("POTScaledOscweight", &POTScaledOscweight);

   TString FHCweightfile = "./dataset/weights/FHCNuTauMVA_BDTA.weights.xml";
   TString RHCweightfile = "./dataset/weights/RHCNuTauMVA_BDTA.weights.xml";

   if (!useRHC) {
      reader->BookMVA("BDT", FHCweightfile);
   }
   else {
      reader->BookMVA("BDT", RHCweightfile);
   }

   float BDTOutput = 0.;
   double cutCounter = 0.;
   double totalNutauSignal = 0.;
   double selectedNutauSignal = 0.;
   double bkgdCounter = 0.;

   int nentries = ch1->GetEntries();

   double totalPOTfromFile = 0.;

   for (auto j = 0; j < nentries; j++) {
      ch1->GetEntry(j);
      filename = ch1->GetCurrentFile()->GetName();

      cvnnue = (float) dcvnnue;
      cvnnumu = (float) dcvnnumu;
      cvnnutau = (float) dcvnnutau;
      cvnnc = (float) dcvnnc;
      Ev = (float) dEv;
      eDepPim = (float) deDepPim;
      nipim = (float) inipim;
      nipip = (float) inipip;
      eDepHad = (eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther);

      /*
      Preliminary code to implement fully reconstructed energy. Right now we only look at eDepHad
      RecoEVisE = eDepHad + RecoLepEnNue;
      RecoEVisM = eDepHad + RecoLepEnNumu;
      RecoEVisT = eDepHad;
      */


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
         if (isCC && nuPDG == 16 && filename == tauswap)
            totalNutauSignal += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);

         //optimalBDTCut comes from the optimal BDT cut given by the TMVA output ROOT file (F/R)HC.MVA
         if (reader->EvaluateMVA("BDT") >= optimalBDTCut) {
            if (isCC) {
               if (filename == nonswap) {
                  cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                  if (nuPDG == 14 || nuPDG == -14) {
                     hWLBkgdM->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hWLBkgdMReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hWLBkgdMRecoVis->Fill(RecoEVisM, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else {
                     hWLBkgdE->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hWLBkgdEReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hWLBkgdERecoVis->Fill(RecoEVisE, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
               }
               else if (filename == nueswap) {
                  cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                  if (nuPDG == 12 || nuPDG == -12) {
                     hWLBkgdE->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hWLBkgdEReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hWLBkgdERecoVis->Fill(RecoEVisE, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else {
                     hIntBkgdT->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hIntBkgdTReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hIntBkgdTRecoVis->Fill(RecoEVisT, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
               }
               else {
                  if (nuPDG == 16) {
                     selectedNutauSignal += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                     cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                     hSignal->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hSignalReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hSignalRecoVis->Fill(RecoEVisT, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else if (nuPDG == -16){
                     cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                     hWSBkgdT->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hWSBkgdTReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hWSBkgdTRecoVis->Fill(RecoEVisT, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else {
                     cutCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                     hWLBkgdM->Fill(Ev, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     hWLBkgdMReco->Fill(eDepHad, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     //hWLBkgdMRecoVis->Fill(RecoEVisM, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
               }
            }
            else {
               if (filename == nonswap) {
                  cutCounter += (scalefactor / totalPOTfromFile);
                  if (nuPDG == 16 || nuPDG == -16) {
                     hNCBkgd->Fill(Ev, (scalefactor / totalPOTfromFile));
                     hNCBkgdReco->Fill(eDepHad, (scalefactor / totalPOTfromFile));
                     //hNCBkgdRecoVis->Fill(RecoEVisT, (scalefactor / totalPOTfromFile));
                  }
                  else if (nuPDG == 14 || nuPDG == -14) {
                     hNCBkgd->Fill(Ev, (scalefactor / totalPOTfromFile));
                     hNCBkgdReco->Fill(eDepHad, (scalefactor / totalPOTfromFile));
                     //hNCBkgdRecoVis->Fill(RecoEVisM, (scalefactor / totalPOTfromFile));
                  }
                  else {
                     hNCBkgd->Fill(Ev, (scalefactor / totalPOTfromFile));
                     hNCBkgdReco->Fill(eDepHad, (scalefactor / totalPOTfromFile));
                     //hNCBkgdRecoVis->Fill(RecoEVisE, (scalefactor / totalPOTfromFile));
                  }
               }
            }
         }
         else {
            if (isCC) {
                if (filename == nonswap)
                   bkgdCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                else if (filename == nueswap)
                   bkgdCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
                else if (filename == tauswap)
                   bkgdCounter += OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile);
            }
            else {
               if (filename == nonswap)
                bkgdCounter += (scalefactor / totalPOTfromFile);
            }
         }
      }
      if (j % 10000 == 0) std::cout << j << "/" << nentries << std::endl;
   }

   double purity = selectedNutauSignal / cutCounter;
   double efficiency = selectedNutauSignal / totalNutauSignal;
   double significance = selectedNutauSignal / Sqrt(cutCounter - selectedNutauSignal);

   std::cout << "totalNutauSignal = " << totalNutauSignal << std::endl;
   std::cout << "selectedNutauSignal = " << selectedNutauSignal << std::endl;
   std::cout << "cutCounter = " << cutCounter << std::endl;
   std::cout << "bkgdCounter = " << bkgdCounter << std::endl;

   std::cout << "Purity = " << purity << std::endl;
   std::cout << "Efficiency = " << efficiency << std::endl;
   std::cout << "S / Sqrt(B) = " << significance << std::endl;
   std::cout << "S / Sqrt(S+B) = " << selectedNutauSignal / Sqrt(cutCounter) << std::endl;

   std::cout << "Signal = " << hSignal->Integral() << '\t' << "Background = " << hNCBkgd->Integral() + hWLBkgdE->Integral() + hWLBkgdM->Integral() +
      hIntBkgdT->Integral() + hWSBkgdT->Integral() << std::endl;



   hSignal->SetLineWidth(2);
   hSignal->SetFillStyle(3001);
   hSignal->SetFillColor(kRed);

   hNCBkgd->SetLineWidth(2);
   hNCBkgd->SetFillStyle(3001);
   hNCBkgd->SetFillColor(kYellow);

   hWLBkgdE->SetLineWidth(2);
   hWLBkgdE->SetFillStyle(3001);
   hWLBkgdE->SetFillColor(kGreen);

   hWLBkgdM->SetLineWidth(2);
   hWLBkgdM->SetFillStyle(3001);
   hWLBkgdM->SetFillColor(kBlue);

   hIntBkgdT->SetLineWidth(2);
   hIntBkgdT->SetFillStyle(3001);
   hIntBkgdT->SetFillColor(kOrange);

   hWSBkgdT->SetLineWidth(2);
   hWSBkgdT->SetFillStyle(3001);
   hWSBkgdT->SetFillColor(kOrange + 2);

   hStack->Add(hNCBkgd, "HIST");
   hStack->Add(hWLBkgdE, "HIST");
   hStack->Add(hWLBkgdM, "HIST");
   hStack->Add(hWSBkgdT, "HIST");
   hStack->Add(hIntBkgdT, "HIST");
   hStack->Add(hSignal, "HIST");


   if (!useRHC) {
      leg->AddEntry(hSignal, "#nu_{#tau} Signal");
      leg->AddEntry(hWSBkgdT, "WS #bar{#nu_{#tau}} Background");
      hStack->SetTitle("FD Events that pass cut of BDT >= 0.1206 (FHC -- 3.5 POT Years)");
   }
   else {
      leg->AddEntry(hSignal, "#bar{#nu_{#tau}} Signal");
      leg->AddEntry(hWSBkgdT, "WS #nu_{#tau} Background");
      hStack->SetTitle("FD Events that pass cut of BDT >= 0.1132 (RHC -- 3.5 POT Years)");
   }

   leg->AddEntry(hNCBkgd, "NC Background");
   leg->AddEntry(hWLBkgdE, "WL #nu_{e} Background");
   leg->AddEntry(hWLBkgdM, "WL #nu_{#mu} Background");
   leg->AddEntry(hIntBkgdT, "Intrinsic Background");

   hStack->Draw("HIST");
   hStack->GetXaxis()->SetTitle("True E_{#nu} (GeV)");
   hStack->GetYaxis()->SetTitle("Weighted Event Count");
   leg->Draw();
   if (!useRHC) {
      c1->SaveAs("./Figures/FHCMVAStackedHistTrueE.png");
   }
   else {
      c1->SaveAs("./Figures/RHCMVAStackedHistTrueE.png");
   }

   c1->Clear();
   c1->Update();

   hSignalReco->SetLineWidth(2);
   hSignalReco->SetFillStyle(3001);
   hSignalReco->SetFillColor(kRed);

   hNCBkgdReco->SetLineWidth(2);
   hNCBkgdReco->SetFillStyle(3001);
   hNCBkgdReco->SetFillColor(kYellow);

   hWLBkgdEReco->SetLineWidth(2);
   hWLBkgdEReco->SetFillStyle(3001);
   hWLBkgdEReco->SetFillColor(kGreen);

   hWLBkgdMReco->SetLineWidth(2);
   hWLBkgdMReco->SetFillStyle(3001);
   hWLBkgdMReco->SetFillColor(kBlue);

   hIntBkgdTReco->SetLineWidth(2);
   hIntBkgdTReco->SetFillStyle(3001);
   hIntBkgdTReco->SetFillColor(kOrange);

   hWSBkgdTReco->SetLineWidth(2);
   hWSBkgdTReco->SetFillStyle(3001);
   hWSBkgdTReco->SetFillColor(kOrange + 2);

   hStackReco->Add(hNCBkgdReco, "HIST");
   hStackReco->Add(hWLBkgdEReco, "HIST");
   hStackReco->Add(hWLBkgdMReco, "HIST");
   hStackReco->Add(hWSBkgdTReco, "HIST");
   hStackReco->Add(hIntBkgdTReco, "HIST");
   hStackReco->Add(hSignalReco, "HIST");


   if (!useRHC) {
      legReco->AddEntry(hSignal, "#nu_{#tau} Signal");
      legReco->AddEntry(hWSBkgdT, "WS #bar{#nu_{#tau}} Background");
      hStackReco->SetTitle("FD Events that pass cut of BDT >= 0.1206 (FHC -- 3.5 POT Years)");
   }
   else {
      legReco->AddEntry(hSignal, "#bar{#nu_{#tau}} Signal");
      legReco->AddEntry(hWSBkgdT, "WS #nu_{#tau} Background");
      hStackReco->SetTitle("FD Events that pass cut of BDT >= 0.1132 (RHC -- 3.5 POT Years)");
   }

   legReco->AddEntry(hNCBkgdReco, "NC Background");
   legReco->AddEntry(hWLBkgdEReco, "WL #nu_{e} Background");
   legReco->AddEntry(hWLBkgdMReco, "WL #nu_{#mu} Background");
   legReco->AddEntry(hIntBkgdTReco, "Intrinsic Background");

   hStackReco->Draw("HIST");
   hStackReco->GetXaxis()->SetTitle("Reco E_{had} (GeV)");
   hStackReco->GetYaxis()->SetTitle("Weighted Event Count");
   leg->Draw();
   if (!useRHC) {
      c1->SaveAs("./Figures/FHCMVAStackedHistRecoEHad.png");
   }
   else {
      c1->SaveAs("./Figures/RHCMVAStackedHistRecoEHad.png");
   }
   c1->Clear();
   c1->Update();



   /*hSignalRecoVis->SetLineWidth(2);
   hSignalRecoVis->SetFillStyle(3001);
   hSignalRecoVis->SetFillColor(kRed);

   hNCBkgdRecoVis->SetLineWidth(2);
   hNCBkgdRecoVis->SetFillStyle(3001);
   hNCBkgdRecoVis->SetFillColor(kYellow);

   hWLBkgdERecoVis->SetLineWidth(2);
   hWLBkgdERecoVis->SetFillStyle(3001);
   hWLBkgdERecoVis->SetFillColor(kGreen);

   hWLBkgdMRecoVis->SetLineWidth(2);
   hWLBkgdMRecoVis->SetFillStyle(3001);
   hWLBkgdMRecoVis->SetFillColor(kBlue);

   hIntBkgdTRecoVis->SetLineWidth(2);
   hIntBkgdTRecoVis->SetFillStyle(3001);
   hIntBkgdTRecoVis->SetFillColor(kOrange);

   hWSBkgdTRecoVis->SetLineWidth(2);
   hWSBkgdTRecoVis->SetFillStyle(3001);
   hWSBkgdTRecoVis->SetFillColor(kOrange + 2);

   hStackRecoVis->Add(hNCBkgdRecoVis, "HIST");
   hStackRecoVis->Add(hWLBkgdERecoVis, "HIST");
   hStackRecoVis->Add(hWLBkgdMRecoVis, "HIST");
   hStackRecoVis->Add(hWSBkgdTRecoVis, "HIST");
   hStackRecoVis->Add(hIntBkgdTRecoVis, "HIST");
   hStackRecoVis->Add(hSignalRecoVis, "HIST");

   if (!useRHC) {
      legRecoVis->AddEntry(hSignal, "#nu_{#tau} Signal");
      legRecoVis->AddEntry(hWSBkgdT, "WS #bar{#nu_{#tau}} Background");
      hStackRecoVis->SetTitle("FD Events that pass cut of BDT >= 0.1206 (FHC -- 3.5 POT Years)");
   }
   else {
      legRecoVis->AddEntry(hSignal, "#bar{#nu_{#tau}} Signal");
      legRecoVis->AddEntry(hWSBkgdT, "WS #nu_{#tau} Background");
      hStackRecoVis->SetTitle("FD Events that pass cut of BDT >= 0.1132 (RHC -- 3.5 POT Years)");
   }

   legRecoVis->AddEntry(hNCBkgdRecoVis, "NC Background");
   legRecoVis->AddEntry(hWLBkgdERecoVis, "WL #nu_{e} Background");
   legRecoVis->AddEntry(hWLBkgdMRecoVis, "WL #nu_{#mu} Background");
   legRecoVis->AddEntry(hIntBkgdTRecoVis, "Intrinsic Background");

   hStackRecoVis->Draw("HIST");
   hStackRecoVis->GetXaxis()->SetTitle("Reco E_{vis} (GeV)");
   hStackRecoVis->GetYaxis()->SetTitle("Weighted Event Count");
   leg->Draw();

   if (!useRHC) {
      c1->SaveAs("./Figures/FHCMVAStackedHistRecoEVis.png");
   }
   else {
      c1->SaveAs("./Figures/RHCMVAStackedHistRecoEVis.png");
   }

   c1->Clear();
   c1->Update();*/


}
