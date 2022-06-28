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

TString nonswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nonswap.root";
TString nueswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nueswap.root";
TString tauswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_tauswap.root";



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

//This function is what I used to calculate the POT numbers for each MC swapfile
void POTNumber(void) {


  TChain* ch1 = new TChain("meta");
  ch1->Add(nonswap);
  ch1->Add(nueswap);
  ch1->Add(tauswap);

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
    if (ch1->GetCurrentFile()->GetName() == nonswap)
      nonswapPOT += POT;
    else if (ch1->GetCurrentFile()->GetName() == nueswap)
      nueswapPOT += POT;
    else if (ch1->GetCurrentFile()->GetName() == tauswap)
      tauswapPOT += POT;
    else
      std::cerr << "This should not happen (file name not recognized) \n";
  }

  std::cout << "nonswapPOT = " << nonswapPOT << '\n';
  std::cout << "nueswapPOT = " << nueswapPOT << '\n';
  std::cout << "tauswapPOT = " << tauswapPOT << '\n';
}


void CutTest(void) {
   gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetRightMargin(0.1);
  c1->SetLeftMargin(0.15);


  TChain* ch1 = new TChain("cafTree");


  ch1->Add(nonswap);
  ch1->Add(nueswap);
  ch1->Add(tauswap);

  THStack* hStack = new THStack("hStack", "");
  TLegend *leg = new TLegend(0.5,0.45, 0.8, 0.8);
  leg->SetBorderSize(0);

  TH1D* hSignal = new TH1D("NuTau Signal", "NuTau Signal", 25, 0., 10.);
  TH1D* hNCBkgd = new TH1D("NCBkgd", "NCBkgd", 25, 0., 10.);
  TH1D* hWLBkgdE = new TH1D("WLBkgdE", "WLBkgdE", 25, 0., 10.);
  TH1D* hWLBkgdM = new TH1D("WLBkgdM", "WLBkgdM", 25, 0., 10.);
  TH1D* hWLBkgdT = new TH1D("WLBkgdT", "WLBkgdT", 25, 0., 10.);

  TH1D hPurity = TH1D("hPurity", "hPurity", 25, 0., 10.);
  TH1D hEfficiency = TH1D("hEfficiency", "hEfficiency", 25, 0., 10.);
  TH1D hSignificance = TH1D("hSignificance", "hSignificance", 25, 0., 5.);

  TH1D* hTrueSignal = new TH1D("hTrueSignal", "hTrueSignal", 25, 0., 10.);
  TH1D* hTrueSignalCut = new TH1D("hTrueSignalCut", "hTrueSignalCut", 25, 0., 10.);
  TH1D* hCutCounter = new TH1D("hCutCounter", "hCutCounter", 25, 0., 10.);


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

  //Scaling
  double years = 3.5; //Exposure time
  double POTperYear = 1.1e21;
  double nonswapPOT = 1.62824e24;
  double nueswapPOT = 1.64546e24;
  double tauswapPOT = 5.18551e24;

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


  int nentries = ch1->GetEntries();



  /*
    Purity: True signal that passes the cut / Events that pass the cut
    Efficiency: True signal that passes the cut/True signal
    Sensitivity: True signal that passes the cut /sqrt(Events the pass the cut)
  */
  double efficiency = 0.;
  double purity = 0.;
  double sensitivity = 0.;

  double trueSignal = 0.; //CC Nutaus that come from tauswap
  double cutCounter = 0.; //Anything that passes the cut
  double trueSignalCut = 0.; //This is for true signal that passes the cut

  for (auto i = 0; i < nentries; i++) {
    ch1->GetEntry(i);
    eDepTotal = eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther + RecoLepEnNumu + RecoLepEnNue;
    eDepHad = eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther;
    DepPimRatio = eDepPim / eDepHad;

    if (cvnnue < 0.85 && cvnnumu < 0.5) {
      if (ch1->GetCurrentFile()->GetName() == nonswap) {
        if (isCC) {
          if (nuPDG == 12 || nuPDG == -12)
            hWLBkgdE->Fill(eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nonswapPOT));
          else
            hWLBkgdM->Fill(eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nonswapPOT));
        }
        else
          hNCBkgd->Fill(eDepTotal, (1./3.) * (years * POTperYear / nonswapPOT));
      }
      else if (ch1->GetCurrentFile()->GetName() == nueswap) {
        if (isCC) {
          if (nuPDG == 12 || nuPDG == -12)
            hWLBkgdE->Fill(eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nueswapPOT));
          else
            hWLBkgdT->Fill(eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nueswapPOT));
        }
        else
          hNCBkgd->Fill(eDepTotal, (1./3.) * (years * POTperYear / nueswapPOT));
      }
      else if (ch1->GetCurrentFile()->GetName() == tauswap) {
        if (isCC) {
          if (nuPDG == 16 || nuPDG == -16)
            hSignal->Fill(eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT));
          else
            hWLBkgdM->Fill(eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT));
        }
        else
          hNCBkgd->Fill(eDepTotal, (1./3.) * (years * POTperYear / tauswapPOT));
      }
      else
        std::cerr << "This should not happen (file name not recognized) \n";
    }

    //Currently configured for nu_es

     if (nuPDG == 16 && isCC == 1 && ch1->GetCurrentFile()->GetName() == tauswap) {
         trueSignal += OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT);
         hTrueSignal->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT));
     }
     if (cvnnue < 0.85 && cvnnumu < 0.5) {
       if (ch1->GetCurrentFile()->GetName() == nonswap) {
         if (isCC){
           cutCounter += OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nonswapPOT);
           hCutCounter->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nonswapPOT));
        }
         else {
           cutCounter += (1./3.) * (years * POTperYear / tauswapPOT);
           hCutCounter->Fill(Ev, (1./3.) * (years * POTperYear / tauswapPOT));
         }
       }
       else if (ch1->GetCurrentFile()->GetName() == nueswap) {
         if (isCC) {
           cutCounter += OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nueswapPOT);
           hCutCounter->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / nueswapPOT));
         }
         else {
           cutCounter += (1./3.) * (years * POTperYear / nueswapPOT);
           hCutCounter->Fill(Ev, (1./3.) * (years * POTperYear / nueswapPOT));
        }
       }
       else if (ch1->GetCurrentFile()->GetName() == tauswap) {
         if (isCC) {
           if (nuPDG == 16){
             trueSignalCut += OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT);
             cutCounter += OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT);

             hTrueSignalCut->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT));
             hCutCounter->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT));
           }
           else
             cutCounter += OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT);
             hCutCounter->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG) * (years * POTperYear / tauswapPOT));
         }
         else
           cutCounter += (1. / 3.) * (years * POTperYear / tauswapPOT);
           hCutCounter->Fill(Ev, (1. / 3.) * (years * POTperYear / tauswapPOT));
       }
       else {
         std::cerr << "This should not happen (file name not recognized)";
       }
     }
  if (i % 10000 == 0)
    std::cerr << 100 * (float)(i + 1) / (float) nentries << "%" << '\n';
  }

  cutCounter *= 40. / 1.13;
  trueSignal *= 40. / 1.13;
  trueSignalCut *= 40. / 1.13;

//Cutcounter is S + B
  efficiency = trueSignalCut / trueSignal;
  purity = trueSignalCut / cutCounter;
  sensitivity = trueSignalCut / Sqrt(cutCounter - trueSignalCut); //S / Sqrt(B)

  hEfficiency = (* hTrueSignalCut) / (* hTrueSignal);
  hPurity = (* hTrueSignalCut) / (* hCutCounter);


  std::cout << "trueSignalCut = " << trueSignalCut << '\n';
  std::cout << "trueSignal = " << trueSignal << '\n';
  std::cout << "cutCounter = " << cutCounter << '\n';


  std::cout << "Purity = " << purity << '\n';
  std::cout << "Efficiency = " << efficiency << '\n';
  std::cout << "Sensitivity = " << sensitivity << '\n';

  c1->cd();
  /*
  hRecoEVis->SetLineWidth(2);
  hRecoEVis->SetTitle("Reconstruted Energy FHC");
  hRecoEVis->GetXaxis()->SetTitle("Reco E_vis (GeV)");
  hRecoEVis->Draw("HIST");
  c1->SaveAs("./Figures/RecoEnergyhisto.png");
  c1->Clear();
  */

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

  hWLBkgdT->SetLineWidth(2);
  hWLBkgdT->SetFillStyle(3001);
  hWLBkgdT->SetFillColor(kOrange);

  hStack->Add(hNCBkgd, "HIST");
  hStack->Add(hWLBkgdE, "HIST");
  hStack->Add(hWLBkgdM, "HIST");
  hStack->Add(hWLBkgdT, "HIST");
  hStack->Add(hSignal, "HIST");



  leg->AddEntry(hSignal, "Nu_tau Signal");
  leg->AddEntry(hNCBkgd, "NC Background");
  leg->AddEntry(hWLBkgdE, "WL #e^{-} Background");
  leg->AddEntry(hWLBkgdM, "WL #mu^{-} Background");
  leg->AddEntry(hWLBkgdT, "WL #tau^{-} Background");

  hStack->Draw("HIST");
  hStack->SetTitle("FHC MC Events that pass cut of (cvnnue < 0.85 && cvnnumu < 0.5)");
  hStack->GetXaxis()->SetTitle("Reco E (GeV)");
  leg->Draw();
  c1->SaveAs("./Figures/StackedHistTest.png");
  c1->Clear();
  c1->Update();

  hEfficiency.SetLineWidth(2);
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
  c1->Update();

}
