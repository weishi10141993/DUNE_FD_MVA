
#include <iostream>
#include <cstdlib>
#include <string>

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

void MVAChi2Test(void) {

   bool useRHC = false;
   bool useEHadBinning = true;

   double nonswapPOT = 0.;
   double nueswapPOT = 0.;
   double tauswapPOT = 0.;

   TString nonswap = "";
   TString nueswap = "";
   TString tauswap = "";

   double optimalBDTCut = 0.;

   if (!useRHC) {
      nonswapPOT = FHCnonswapPOT;
      nueswapPOT = FHCnueswapPOT;
      tauswapPOT = FHCtauswapPOT;

      nonswap = FHCnonswap;
      nueswap = FHCnueswap;
      tauswap = FHCtauswap;

      optimalBDTCut = 0.1132;
   }
   else {
      nonswapPOT = RHCnonswapPOT;
      nueswapPOT = RHCnueswapPOT;
      tauswapPOT = RHCtauswapPOT;

      nonswap = RHCnonswap;
      nueswap = RHCnueswap;
      tauswap = RHCtauswap;

      optimalBDTCut = 0.1206;
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

   float Ev = 0.;
   int nuPDG = 0.;
   int isCC = -1.;

   float cvnnue = 0.;
   float cvnnumu = 0.;
   float cvnnutau = 0.;
   float cvnnc = 0.;
   float eDepPim = 0.;
   float eDepHad = 0.;
   float nipim = 0.;
   float POTScaledOscweight = 0.;
   float nipip;


   double dEv = 0.;
   double dcvnnue = 0.;
   double dcvnnumu = 0.;
   double dcvnnutau = 0.;
   double dcvnnc = 0.;
   double deDepPim = 0.;
   double deDepHad = 0.;
   int inipim = 0;
   int inipip = 0;

   double eDepP = 0.;
   double eDepN = 0.;
   double eDepPip = 0.;
   double eDepPi0 = 0.;
   double eDepOther = 0.;
   double eDepTotal = 0.;
   double DepPimRatio = 0.;
   double RecoLepEnNumu = 0.;
   double RecoLepEnNue = 0.;

   double eDepLep = 0.;
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


   std::vector<double> binContents = {};
   int nentries = ch1->GetEntries();

   TMVA::Reader* reader = new TMVA::Reader();

   reader->AddVariable("cvnnue", &cvnnue);
   reader->AddVariable("cvnnumu", &cvnnumu);
   reader->AddVariable("cvnnutau", &cvnnutau);
   reader->AddVariable("EHad := eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther", &eDepHad);
   reader->AddVariable("eDepPim", &eDepPim);
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

   float BDToutput = 0.;

   std::vector<double> theta_23Values = {};
   std::vector<double> delmsqValues = {};
   double LHRatio = 0.;

   for (int i = 0; i <= 24; i++) {
      theta_23Values.push_back(0.689 + (0.908 - 0.689) * ((double)(i) / 24.));
      delmsqValues.push_back(2.431e-3 + (2.599e-3 - 2.431e-3) * ((double)(i) / 24.));
   }

   std::vector<THStack> stacked_hists = {};
   std::vector<TLegend> legends = {};
   std::vector<TH1D> signal_hists = {};
   std::vector<TH1D> WLBkgdTwo_hists = {};
   std::vector<TH1D> WLBkgdOne_hists = {};
   std::vector<TH1D> WSBkgd_hists = {};
   std::vector<TH1D> IntBkgd_hists = {};
   std::vector<TH1D> NCBkgd_hists = {};
   std::vector<TH1D> Total_hists = {};
   std::vector<TH1D> Bkgd_hists = {};

   THStack* hDefaultStack = new THStack("TotalStack", "");
   TLegend* leg = new TLegend(0.65,0.65, 0.9, 0.9);
   TH1D* hDefaultSignal = new TH1D("TotalSignal", "TotalSignal", 75, 0., 30.);
   TH1D* hDefaultWLBkgdTwo = new TH1D("TotalWLBkgdTwo", "TotalWLBkgdTwo", 75, 0., 30.);
   TH1D* hDefaultWLBkgdOne = new TH1D("TotalWLBkgdOne", "TotalWLBkgdOne", 75, 0., 30.);
   TH1D* hDefaultWSBkgd = new TH1D("TotalWSBkgd", "TotalWSBkgd", 75, 0., 30.);
   TH1D* hDefaultIntBkgd = new TH1D("TotalIntBkgd", "TotalIntBkgd", 75, 0., 30.);
   TH1D* hDefaultNCBkgd = new TH1D("TotalNCBkgd", "TotalNCBkgd", 75, 0., 30.);
   TH1D hDefaultHist = TH1D("TotalHist", "TotalHist", 75, 0., 30.);



   //Entries 0 through 24 vary delmsq and 25 through 49 vary theta_23
   for (int i = 0; i < 50; i++) {
      THStack* dummyStack = new THStack(TString::Format("Stack%d", i), "");
      stacked_hists.push_back(*dummyStack);

      TH1D* dummySignal = new TH1D(TString::Format("Signal%d", i), TString::Format("Signal%d", i), 75, 0., 30.);
      signal_hists.push_back(*dummySignal);

      TH1D* dummyWLBkgdTwo = new TH1D(TString::Format("WLBkgdTwo%d", i), TString::Format("WLBkgdTwo%d", i), 75, 0., 30.);
      WLBkgdTwo_hists.push_back(*dummyWLBkgdTwo);

      TH1D* dummyWLBkgdOne = new TH1D(TString::Format("WLBkgdOne%d", i), TString::Format("WLBkgdOne%d", i), 75, 0., 30.);
      WLBkgdOne_hists.push_back(*dummyWLBkgdOne);

      TH1D* dummyWSBkgd = new TH1D(TString::Format("WSBkgd%d", i), TString::Format("WSBkgd%d", i), 75, 0., 30.);
      WSBkgd_hists.push_back(*dummyWSBkgd);

      TH1D* dummyIntBkgd = new TH1D(TString::Format("IntBkgd%d", i), TString::Format("IntBkgd%d", i), 75, 0., 30.);
      IntBkgd_hists.push_back(*dummyIntBkgd);

      TH1D* dummyNCBkgd = new TH1D(TString::Format("NCBkgd%d", i), TString::Format("NCBkgd%d", i), 75, 0., 30.);
      NCBkgd_hists.push_back(*dummyNCBkgd);

      TH1D* dummyPlaceholder = new TH1D(TString::Format("Total%d", i), TString::Format("Total%d", i), 75, 0., 30.);
      Total_hists.push_back(*dummyPlaceholder);

      TH1D* dummyBkgdHist = new TH1D(TString::Format("Bkgd%d", i), TString::Format("Bkgd%d", i), 75, 0., 30.);
      Bkgd_hists.push_back(*dummyBkgdHist);

      TLegend* dummyleg = new TLegend(0.65,0.65, 0.9, 0.9);
      legends.push_back(*dummyleg);
   }

   float* binVariable = NULL;
   if (!useEHadBinning) {
      binVariable = &Ev;
   }
   else {
      binVariable = &eDepHad;
   }

   double totalPOTfromFile = 0.;

   for (auto j = 0; j < nentries; j++) {
      ch1->GetEntry(j);

      cvnnue = (float) dcvnnue;
      cvnnumu = (float) dcvnnumu;
      cvnnutau = (float) dcvnnutau;
      cvnnc = (float) dcvnnc;
      Ev = (float) dEv;
      eDepPim = (float) deDepPim;
      nipim = (float) inipim;
      nipip = (float) inipip;
      eDepHad = (eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther);

      filename = ch1->GetCurrentFile()->GetName();

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

         if (reader->EvaluateMVA("BDT") >= optimalBDTCut) {
            if (isCC) {
               if (filename == nonswap) {
                  if (nuPDG == 14 || nuPDG == -14) {
                     for (int i = 0; i < 25; i++) {
                        WLBkgdOne_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        WLBkgdOne_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultWLBkgdOne->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else {
                     for (int i = 0; i < 25; i++) {
                        WLBkgdTwo_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        WLBkgdTwo_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultWLBkgdTwo->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
               }
               else if (filename == nueswap) {
                  if (nuPDG == 12 || nuPDG == -12) {
                     for (int i = 0; i < 25; i++) {
                        WLBkgdTwo_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        WLBkgdTwo_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultWLBkgdTwo->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else {
                     for (int i = 0; i < 25; i++) {
                        IntBkgd_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        IntBkgd_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultIntBkgd->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
               }
               else {
                  if (nuPDG == 16) {
                     for (int i = 0; i < 25; i++) {
                        signal_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        signal_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultSignal->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else if (nuPDG == -16){
                     for (int i = 0; i < 25; i++) {
                        WSBkgd_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        WSBkgd_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultWSBkgd->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
                  else {
                     for (int i = 0; i < 25; i++) {
                        WLBkgdOne_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, delmsqValues[i], THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     for (int i = 25; i < 50; i++) {
                        WLBkgdOne_hists[i].Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, theta_23Values[i - 25], THETA_13) * (scalefactor / totalPOTfromFile));
                     }
                     hDefaultWLBkgdOne->Fill(*binVariable, OscWeight(filename, Ev, nuPDG, DELMSQ_31, THETA_23, THETA_13) * (scalefactor / totalPOTfromFile));
                  }
               }
            }
            else {
               //This was split up by PDG because I once had different definitions for Reco E for each PDG. Now Im too lazy to change it.
               if (filename == nonswap) {
                  hDefaultNCBkgd->Fill(*binVariable, (scalefactor / totalPOTfromFile));
                  for (int i = 0; i < 25; i++) {
                     NCBkgd_hists[i].Fill(*binVariable, scalefactor / totalPOTfromFile);
                  }
                  for (int i = 25; i < 50; i++) {
                     NCBkgd_hists[i].Fill(*binVariable, scalefactor / totalPOTfromFile);
                  }
               }
            }
         }
      }
      if (j % 10000 == 0) std::cout << j << "/" << nentries << std::endl;
   }

   hDefaultHist = (*hDefaultSignal) + (*hDefaultWLBkgdTwo) + (*hDefaultWLBkgdOne) + (*hDefaultWSBkgd) + (*hDefaultIntBkgd) + (*hDefaultNCBkgd);

   for (int bin = 0; bin <= hDefaultHist.GetNcells(); bin++) {
      /*
      std::cout << hDefaultSignal->GetBinContent(bin) << std::endl;
      std::cout << hDefaultWLBkgdTwo->GetBinContent(bin) << std::endl;
      std::cout << hDefaultWLBkgdOne->GetBinContent(bin) << std::endl;
      std::cout << hDefaultWSBkgd->GetBinContent(bin) << std::endl;
      std::cout << hDefaultIntBkgd->GetBinContent(bin) << std::endl;
      std::cout << hDefaultNCBkgd->GetBinContent(bin) << std::endl;
      */

     binContents.push_back(hDefaultHist.GetBinContent(bin));
   }

   hDefaultSignal->SetLineWidth(2);
   hDefaultSignal->SetFillStyle(3001);
   hDefaultSignal->SetFillColor(kRed);

   hDefaultNCBkgd->SetLineWidth(2);
   hDefaultNCBkgd->SetFillStyle(3001);
   hDefaultNCBkgd->SetFillColor(kYellow);

   hDefaultWLBkgdTwo->SetLineWidth(2);
   hDefaultWLBkgdTwo->SetFillStyle(3001);
   hDefaultWLBkgdTwo->SetFillColor(kGreen);

   hDefaultWLBkgdOne->SetLineWidth(2);
   hDefaultWLBkgdOne->SetFillStyle(3001);
   hDefaultWLBkgdOne->SetFillColor(kBlue);

   hDefaultIntBkgd->SetLineWidth(2);
   hDefaultIntBkgd->SetFillStyle(3001);
   hDefaultIntBkgd->SetFillColor(kOrange);

   hDefaultWSBkgd->SetLineWidth(2);
   hDefaultWSBkgd->SetFillStyle(3001);
   hDefaultWSBkgd->SetFillColor(kOrange + 2);


   hDefaultStack->Add(hDefaultNCBkgd, "HIST");
   hDefaultStack->Add(hDefaultWLBkgdTwo, "HIST");
   hDefaultStack->Add(hDefaultWLBkgdOne, "HIST");
   //hDefaultStack->Add(hDefaultIntBkgd, "HIST");
   hDefaultStack->Add(hDefaultWSBkgd, "HIST");
   hDefaultStack->Add(hDefaultSignal, "HIST");

   if (!useRHC) {
      leg->AddEntry(hDefaultSignal, "#nu_{#tau} Signal");
      leg->AddEntry(hDefaultWSBkgd, "WS #bar{#nu_{#tau}} Background");
      hDefaultStack->SetTitle("FD Events that pass cut of BDT >= 0.1206 (FHC -- 3.5 POT Years)");
   }
   else {
      leg->AddEntry(hDefaultSignal, "#bar{#nu_{#tau}} Signal");
      leg->AddEntry(hDefaultWSBkgd, "WS #nu_{#tau} Background");
      hDefaultStack->SetTitle("FD Events that pass cut of BDT >= 0.1132 (RHC -- 3.5 POT Years)");
   }

   leg->AddEntry(hDefaultNCBkgd, "NC Background");
   leg->AddEntry(hDefaultWLBkgdTwo, "WL #nu_{e} Background");
   leg->AddEntry(hDefaultWLBkgdOne, "WL #nu_{#mu} Background");
   //leg->AddEntry(hDefaultIntBkgd, "Intrinsic Background");

   hDefaultStack->Draw("HIST");
   if (!useRHC) {
      hDefaultStack->SetTitle("#nu_{#tau} Selection with BDT (FHC -- 3.5 POT Years)");
   }
   else {
      hDefaultStack->SetTitle("#bar{#nu_{#tau}} Selection with BDT (RHC -- 3.5 POT Years)");
   }

   if (!useEHadBinning) {
      hDefaultStack->GetXaxis()->SetTitle("True E_{#nu} (GeV)");
   }
   else {
      hDefaultStack->GetXaxis()->SetTitle("TReco. E_{had} (GeV)");
   }

   hDefaultStack->GetYaxis()->SetTitle("Counts");

   //These values are done in a different macro, so I had to input these in manually
   //It wouldn't be too hard to extend this macro to compute these with each run

   TLatex *tb1, *tb2, *tb3 = NULL;

   if (!useRHC) {
      tb1 = new TLatex(5.8, 32, "#scale[0.7]{Efficiency: 0.41387}");
      tb2 = new TLatex(5.8, 28, "#scale[0.7]{Purity: 0.213314}");
      tb3 = new TLatex(5.8, 24, "#scale[0.7]{Total S / #sqrt{B}: 5.4697}");
   }
   else {
      tb1 = new TLatex(0.5, 10, "#scale[0.5]{Efficiency: 0.325302}");
      tb2 = new TLatex(0.5, 9.5, "#scale[0.5]{Purity: 0.150121}");
      tb3 = new TLatex(0.5, 9, "#scale[0.5]{Total S / #sqrt{B}: 2.2944}");
   }


   tb1->Draw();
   tb2->Draw();
   tb3->Draw();

   leg->Draw();
   if (!useRHC) {
      if (!useEHadBinning) {
         c1->SaveAs("./Figures/FHCMVAStackedHistTrueE.png");
      }
      else {
         c1->SaveAs("./Figures/FHCMVAStackedHistRecoEHad.png");
      }
   }
   else {
      if (!useEHadBinning) {
         c1->SaveAs("./Figures/RHCMVAStackedHistTrueE.png");
      }
      else {
         c1->SaveAs("./Figures/FHCMVAStackedHistRecoEHad.png");
      }
   }
   c1->Clear();
   c1->Update();

   for (int i = 0; i < 50; i++) {
      Total_hists[i] = signal_hists[i] + WLBkgdTwo_hists[i] + WLBkgdOne_hists[i] + WSBkgd_hists[i] + IntBkgd_hists[i] + NCBkgd_hists[i];
      Bkgd_hists[i] = WLBkgdTwo_hists[i] + WLBkgdOne_hists[i] + WSBkgd_hists[i] + IntBkgd_hists[i] + NCBkgd_hists[i];

      signal_hists[i].SetLineWidth(2);
      signal_hists[i].SetFillStyle(3001);
      signal_hists[i].SetFillColor(kRed);

      NCBkgd_hists[i].SetLineWidth(2);
      NCBkgd_hists[i].SetFillStyle(3001);
      NCBkgd_hists[i].SetFillColor(kYellow);

      WLBkgdTwo_hists[i].SetLineWidth(2);
      WLBkgdTwo_hists[i].SetFillStyle(3001);
      WLBkgdTwo_hists[i].SetFillColor(kGreen);

      WLBkgdOne_hists[i].SetLineWidth(2);
      WLBkgdOne_hists[i].SetFillStyle(3001);
      WLBkgdOne_hists[i].SetFillColor(kBlue);

      IntBkgd_hists[i].SetLineWidth(2);
      IntBkgd_hists[i].SetFillStyle(3001);
      IntBkgd_hists[i].SetFillColor(kOrange);

      WSBkgd_hists[i].SetLineWidth(2);
      WSBkgd_hists[i].SetFillStyle(3001);
      WSBkgd_hists[i].SetFillColor(kOrange + 2);

      stacked_hists[i].Add((TH1D*)&NCBkgd_hists[i], "HIST");
      stacked_hists[i].Add((TH1D*)&IntBkgd_hists[i], "HIST");
      stacked_hists[i].Add((TH1D*)&WSBkgd_hists[i], "HIST");
      stacked_hists[i].Add((TH1D*)&WLBkgdOne_hists[i], "HIST");
      stacked_hists[i].Add((TH1D*)&WLBkgdTwo_hists[i], "HIST");
      stacked_hists[i].Add((TH1D*)&signal_hists[i], "HIST");

      legends[i].AddEntry((TH1D*)&signal_hists[i], "#nu_{#tau} Signal");
      legends[i].AddEntry((TH1D*)&NCBkgd_hists[i], "NC Background");
      legends[i].AddEntry((TH1D*)&WLBkgdTwo_hists[i], "WL #nu_{e} Background");
      legends[i].AddEntry((TH1D*)&WLBkgdOne_hists[i], "WL #nu_{#mu} Background");
      legends[i].AddEntry((TH1D*)&WSBkgd_hists[i], "WS #bar{#nu_{#tau}} Background");
      legends[i].AddEntry((TH1D*)&IntBkgd_hists[i], "Intrinsic Background");
   }

   TGraph* delmsqChi2Graph = new TGraph();
   TGraph* delmsqZscore = new TGraph();
   double mu2 = 0.;
   double logval = 0.;
   double n = 0.;
   double pvalue = 0.;
   double Zscore = 0.;
   for (int i = 0; i < 25; i++) {
     std::cout << "delmsq = " << delmsqValues[i] << "\t" << "theta_23 = DEFAULT (0.859)" << std::endl;
     for (int j = 1; j <= 75; j++) {
         logval = 0.;
         mu2 = Total_hists[i].GetBinContent(j);
         n = Bkgd_hists[i].GetBinContent(j);
         if (mu2 == 0.)
            logval = 0.;
         else
            logval = Log(n / mu2);
         //std::cout << "mu2 = " << mu2 << std::endl;
         //std::cout << "binContents[" << j << "] = " << n << std::endl;
         if (n != 0.)
            LHRatio += 2. * (mu2 - n + (n * logval));
         else
            LHRatio += 2. * mu2;


         //std::cout << "Chi2 = " << LHRatio << '\t' << "Pvalue = " << pvalue << std::endl;
     }
     pvalue = Prob(LHRatio, 75);
     Zscore = Sqrt2() * ErfInverse(1 - 2 * pvalue);
     delmsqChi2Graph->SetPoint(i, delmsqValues[i], LHRatio);
     delmsqZscore->SetPoint(i, delmsqValues[i], Zscore);
     LHRatio = 0.;
     pvalue = 0.;
     Zscore = 0.;
   }

   TGraph* theta_23Chi2Graph = new TGraph();
   TGraph* theta_23Zscore = new TGraph();
   double mu1 = 0.;
   for (int i = 25; i < 50; i++) {
     std::cout << "delmsq = DEFUALT (2.515e-3)"  << "\t" << "theta_23 = " << theta_23Values[i - 25] << std::endl;
     for (int j = 1; j <= 75; j++) {
         logval = 0.;
         mu1 = Total_hists[i].GetBinContent(j);
         n = Bkgd_hists[i].GetBinContent(j);
         if (mu1 == 0.)
            logval = 0.;
         else
            logval = Log(n / mu1);
         //std::cout << "mu1 = " << mu1 << std::endl;
         //std::cout << "binContents[" << j << "] = " << n << std::endl;
         if (n != 0.)
            LHRatio += 2. * (mu1 - n + (n * logval));
         else
            LHRatio += 2. * mu1;

         std::cout << "Chi2 = " << LHRatio << '\t' << "Pvalue = " << pvalue << std::endl;
     }
     pvalue = Prob(LHRatio, 75);
     Zscore = Sqrt2() * ErfInverse(1 - 2. * pvalue);
     theta_23Chi2Graph->SetPoint(i - 25, Power(Sin(theta_23Values[i - 25]), 2.), LHRatio);
     theta_23Zscore->SetPoint(i - 25, Power(Sin(theta_23Values[i - 25]), 2.), Zscore);
     LHRatio = 0.;
     pvalue = 0.;
     Zscore = 0.;
   }


   if (!useRHC) {
      if (!useEHadBinning) {
         theta_23Chi2Graph->SetTitle("FHC #Delta #chi^{2} (BDT cut -- True E. Binning)");
      }
      else {
         theta_23Chi2Graph->SetTitle("FHC #Delta #chi^{2} (BDT cut -- E_{Had} Binning)");
      }
   }
   else {
      if (!useEHadBinning) {
         theta_23Chi2Graph->SetTitle("RHC #Delta #chi^{2} (BDT cut -- True E. Binning)");
      }
      else {
         theta_23Chi2Graph->SetTitle("RHC #Delta #chi^{2} (BDT cut -- E_{Had} Binning)");
      }
   }

   theta_23Chi2Graph->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
   theta_23Chi2Graph->GetYaxis()->SetTitle("#chi^{2}");
   theta_23Chi2Graph->SetMarkerSize(1);
   theta_23Chi2Graph->SetMarkerStyle(21);
   theta_23Chi2Graph->Draw("AP");

   if (!useRHC) {
      if (!useEHadBinning) {
         c1->SaveAs("./Figures/Chi2/MVAFHCTheta23Chi2True.png");
      }
      else {
         c1->SaveAs("./Figures/Chi2/MVAFHCTheta23Chi2Reco.png");
      }
   }
   else {
      if (!useEHadBinning) {
         c1->SaveAs("./Figures/Chi2/MVARHCTheta23Chi2Reco.png");
      }
      else {
         c1->SaveAs("./Figures/Chi2/MVARHCTheta23Chi2Reco.png");
      }
   }

   c1->Clear();
   c1->Update();


   if (!useRHC) {
      if (!useEHadBinning) {
         delmsqChi2Graph->SetTitle("FHC #Delta #chi^{2} (BDT cut -- True E. Binning)");
      }
      else {
         delmsqChi2Graph->SetTitle("FHC #Delta #chi^{2} (BDT cut -- E_{Had} Binning)");
      }
   }
   else {
      if (!useEHadBinning) {
         delmsqChi2Graph->SetTitle("RHC #Delta #chi^{2} (BDT cut -- True E. Binning)");
      }
      else {
         delmsqChi2Graph->SetTitle("RHC #Delta #chi^{2} (BDT cut -- E_{Had} Binning)");
      }
   }

   delmsqChi2Graph->GetXaxis()->SetTitle("#Delta m^{2} (eV^{2})");
   delmsqChi2Graph->GetXaxis()->SetNdivisions(506);
   delmsqChi2Graph->GetYaxis()->SetTitle("#chi^{2}");
   delmsqChi2Graph->SetMarkerSize(1);
   delmsqChi2Graph->SetMarkerStyle(21);
   delmsqChi2Graph->Draw("AP");

   if (!useRHC) {
      if (!useEHadBinning) {
         c1->SaveAs("./Figures/Chi2/MVAFHCDeltaM2Chi2True.png");
      }
      else {
         c1->SaveAs("./Figures/Chi2/MVAFHCDeltaM2Chi2Reco.png");
      }
   }
   else {
      if (!useEHadBinning) {
         c1->SaveAs("./Figures/Chi2/MVARHCDeltaM2Chi2Reco.png");
      }
      else {
         c1->SaveAs("./Figures/Chi2/MVARHCDeltaM2Chi2Reco.png");
      }
   }


   c1->Clear();
   c1->Update();

   //Need to add if statements to change titles and filenames for Zscore plots

   /*theta_23Zscore->SetTitle("FHC Z-Score (BDT Cut -- 3.5 POT Years)");
   theta_23Zscore->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
   theta_23Zscore->GetYaxis()->SetTitle("Z-Score");
   theta_23Zscore->SetMarkerStyle(20);
   theta_23Zscore->Draw("AP");
   c1->SaveAs("./Figures/Chi2/MVAFHCTheta23ZScore.png");
   c1->Clear();
   c1->Update();

   delmsqZscore->SetTitle("FHC Z-Score (BDT Cut -- 3.5 POT Years)");
   delmsqZscore->GetXaxis()->SetTitle("#Delta m^{2} (eV^{2})");
   delmsqZscore->GetXaxis()->SetNdivisions(506);
   delmsqZscore->GetYaxis()->SetTitle("Z-Score");
   delmsqZscore->SetMarkerStyle(20);
   delmsqZscore->Draw("AP");
   c1->SaveAs("./Figures/Chi2/MVAFHCDeltaM2ZScore.png");
   c1->Clear();
   c1->Update();*/

   /*for (int i = 0; i < 25; i++) {
      stacked_hists[i].SetTitle(TString::Format("Placeholder%d -- theta_23 = 0.859, #Delta m^{2} = %f", i, delmsqValues[i]));
      //stacked_hists[i].GetXaxis()->SetTitle("True E (GeV)");
      stacked_hists[i].Draw("HIST");
      legends[i].Draw();
      c1->SaveAs(TString::Format("./Figures/Chi2/dummyhists/MVAStackedHist%d.png", i));
      c1->Clear();
      c1->Update();

   }

   for (int i = 25; i < 50; i++) {
      stacked_hists[i].SetTitle(TString::Format("Placeholder%d -- theta_23 = %f, #Delta m^{2} = 2.515e-3", i, theta_23Values[i - 25]));
      //stacked_hists[i].GetXaxis()->SetTitle("True E (GeV)");
      stacked_hists[i].Draw("HIST");
      legends[i].Draw();
      c1->SaveAs(TString::Format("./Figures/Chi2/dummyhists/MVAStackedHist%d.png", i));
      c1->Clear();
      c1->Update();

   }*/

}
