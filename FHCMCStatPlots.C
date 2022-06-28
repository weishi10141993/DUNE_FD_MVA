/*
This code takes the FD MC swap files and plots a variety of statistics


*/







#include <iostream>
#include <string>

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

void FHCMCStatPlots(void) {

  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->SetRightMargin(0.16);

  TChain* ch1 = new TChain("cafTree");
  ch1->Add(nonswap);
  ch1->Add(nueswap);
  ch1->Add(tauswap);



  /*
  TFile* nonswapFile = new TFile(nonswap, "READ");
  TFile* nueswapFile = new TFile(nueswap, "READ");
  TFile* tauswapFile = new TFile(tauswap, "READ");
  */

  TFile* outputFile = new TFile("WeightedCVNFHC.root", "RECREATE");

  TH1D* wcvnnue_fromnue = new TH1D("wcvnnue_fromnue", "wcvnnue_fromnue", 50, 0., 1.);
  TH1D* wcvnnue_fromnumu = new TH1D("wcvnnue_fromnumu", "wcvnnue_fromnumu", 50, 0., 1.);
  TH1D* wcvnnue_fromnutau = new TH1D("wcvnnue_fromnutau", "wcvnnue_fromnutau", 50, 0., 1.);
  TH1D* wcvnnue_fromnc = new TH1D("wcvnnue_fromnc", "wcvnnue_fromnc", 50, 0., 1.);

  TH1D* wcvnnumu_fromnue = new TH1D("wcvnnumu_fromnue", "wcvnnumu_fromnue", 50, 0., 1.);
  TH1D* wcvnnumu_fromnumu = new TH1D("wcvnnumu_fromnumu", "wcvnnumu_fromnumu", 50, 0., 1.);
  TH1D* wcvnnumu_fromnutau = new TH1D("wcvnnumu_fromnutau", "wcvnnumu_fromnutau", 50, 0., 1.);
  TH1D* wcvnnumu_fromnc = new TH1D("wcvnnumu_fromnc", "wcvnnumu_fromnc", 50, 0., 1.);

  TH1D* wcvnnutau_fromnue = new TH1D("wcvnnutau_fromnue", "wcvnnutau_fromnue", 50, 0., 1.);
  TH1D* wcvnnutau_fromnumu = new TH1D("wcvnnutau_fromnumu", "wcvnnutau_fromnumu", 50, 0., 1.);
  TH1D* wcvnnutau_fromnutau = new TH1D("wcvnnutau_fromnutau", "wcvnnutau_fromnutau", 50, 0., 1.);
  TH1D* wcvnnutau_fromnc = new TH1D("wcvnnutau_fromnc", "wcvnnutau_fromnc", 50, 0., 1.);

  TH1D* wcvnnc_fromnue = new TH1D("wcvnnc_fromnue", "wcvnnc_fromnue", 50, 0., 1.);
  TH1D* wcvnnc_fromnumu = new TH1D("wcvnnc_fromnumu", "wcvnnc_fromnumu", 50, 0., 1.);
  TH1D* wcvnnc_fromnutau = new TH1D("wcvnnc_fromnutau", "wcvnnc_fromnutau", 50, 0., 1.);
  TH1D* wcvnnc_fromnc = new TH1D("wcvnnc_fromnc", "wcvnnc_fromnc", 50, 0., 1.);

  TH1D* weDepPimNue = new TH1D("weDepPimNue", "weDepPimNue", 25, 0., 10.);
  TH1D* weDepPimNumu = new TH1D("weDepPimNumu", "weDepPimNumu", 25, 0., 10.);
  TH1D* weDepPimNutau = new TH1D("weDepPimNutau", "weDepPimNutau", 25, 0., 10.);
  TH1D* weDepPimNC= new TH1D("weDepPimNC", "weDepPimNC", 25, 0., 10.);

  TH1D* weDepElseNue = new TH1D("weDepElseNue", "weDepElseNue", 25, 0., 10.);
  TH1D* weDepElseNumu = new TH1D("weDepElseNumu", "weDepElseNumu", 25, 0., 10.);
  TH1D* weDepElseNutau = new TH1D("weDepElseNutau", "weDepElseNutau", 25, 0., 10.);
  TH1D* weDepElseNC = new TH1D("weDepElseNC", "weDepElseNC", 25, 0., 10.);

  TH2D* hDepPimRatioE = new TH2D("hDepPimRatioE", "hDepPimRatioE", 25, 0., 10, 20, 0., 1.);
  TH2D* hDepPimRatioM = new TH2D("hDepPimRatioM", "hDepPimRatioM", 25, 0., 10, 20, 0., 1.);
  TH2D* hDepPimRatioT = new TH2D("hDepPimRatioT", "hDepPimRatioT", 25, 0., 10, 20, 0., 1.);
  TH2D* hDepPimRatioNC = new TH2D("hDepPimRatioNC", "hDepPimRatioNC", 25, 0., 10, 20, 0., 1.);

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

  gStyle->SetOptStat(0);

  //TTree* outputTree = new TTree("WeightedCVNFHC", "WeightedCVNFHC");

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
  double eDepTotal = 0.;
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


  int nentries = ch1->GetEntries();

  /*
  outputTree->Branch("wcvnnue_fromnue", &wcvnnue_fromnue);
  outputTree->Branch("wcvnnue_fromnumu", &wcvnnue_fromnumu);
  outputTree->Branch("wcvnnue_fromnutau", &wcvnnue_fromnutau);
  outputTree->Branch("wcvnnue_fromnc", &wcvnnue_fromnc);

  outputTree->Branch("wcvnnumu_fromnue", &wcvnnue_fromnue);
  outputTree->Branch("wcvnnumu_fromnumu", &wcvnnue_fromnumu);
  outputTree->Branch("wcvnnumu_fromnutau", &wcvnnue_fromnutau);
  outputTree->Branch("wcvnnumu_fromnc", &wcvnnue_fromnc);

  outputTree->Branch("wcvnnutau_fromnue", &wcvnnue_fromnue);
  outputTree->Branch("wcvnnutau_fromnumu", &wcvnnue_fromnumu);
  outputTree->Branch("wcvnnutau_fromnutau", &wcvnnue_fromnutau);
  outputTree->Branch("wcvnnutau_fromnc", &wcvnnue_fromnc);
  */

  //Now create 6 2-D histograms of osc prob vs true neutrino energy
  TH2D* nuetonue = new TH2D("nuetonue", "nuetonue", 25, 0., 10., 20, 0., 1.);
  TH2D* nuetonumu = new TH2D("nuetonumu", "nuetonumu", 25, 0., 10., 20, 0., 1.);
  TH2D* nuetonutau = new TH2D("nuetonutau", "nuetonutau", 25, 0., 10., 20, 0., 1.);

  TH2D* numutonue = new TH2D("numutonue", "numutonue", 25, 0., 10., 20, 0., 1.);
  TH2D* numutonumu = new TH2D("numutonumu", "numutonumu", 25, 0., 10., 20, 0., 1.);
  TH2D* numutonutau = new TH2D("numutonutau", "numutonutau", 25, 0., 10., 20, 0., 1.);




  for (auto i = 0; i < nentries; i++) {
    ch1->GetEntry(i);
    eDepHad = eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther; //Attempt at reconstructed energy
    eDepLep = RecoLepEnNumu + RecoLepEnNue;
    eDepTotal = eDepLep + eDepHad;
    //eRecoTotal = eRecoP + eRecoN + eRecoPip + eRecoPim + eRecoPi0 + eRecoOther;
    DepPimRatio = eDepPim / eDepHad;
//    if (eDepTotal <= 10) {
      if (isCC) {
        //Check Flavor
        if (nuPDG == 12 || nuPDG == -12) {
          wcvnnue_fromnue->Fill(cvnnue, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnumu_fromnue->Fill(cvnnumu, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnutau_fromnue->Fill(cvnnutau, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnc_fromnue->Fill(cvnnc, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          weDepPimNue->Fill(eDepPim, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          weDepElseNue->Fill(eDepHad - eDepPim, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hDepPimRatioE->Fill(Ev, DepPimRatio, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hResolutionE->Fill((Ev - eDepTotal)/Ev, Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hTruevRecoEE->Fill(Ev, eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hResolutionE2->Fill(eDepTotal/Ev, Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

         if (2. <= Ev && Ev <= 3.)
            hResolutionE1D->Fill(eDepTotal/Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));


          if (ch1->GetCurrentFile()->GetName() ==  nonswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside nue nonswap\n";
            nuetonue->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG), OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          }
          else if (ch1->GetCurrentFile()->GetName() == nueswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside nue nueswap\n";
            numutonue->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG), OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          }
          else if (ch1->GetCurrentFile()->GetName() == tauswap) {
            std::cerr << "This should not happen (nue in tauswap) \n";
          }
          else {
            std::cerr << "This should not happen (file name not recognized) \n";
          }

        } else if (nuPDG == 14 || nuPDG == -14) {
          wcvnnue_fromnumu->Fill(cvnnue, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnumu_fromnumu->Fill(cvnnumu, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnutau_fromnumu->Fill(cvnnutau, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnc_fromnumu->Fill(cvnnc, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          weDepPimNumu->Fill(eDepPim, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          weDepElseNumu->Fill(eDepHad - eDepPim, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hDepPimRatioM->Fill(Ev, DepPimRatio, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hResolutionM->Fill((Ev - eDepTotal)/Ev, Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hTruevRecoEM->Fill(Ev, eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hResolutionM2->Fill(eDepTotal/Ev, Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          if (2. <= Ev && Ev <= 3.)
            hResolutionM1D->Fill(eDepTotal/Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          if (ch1->GetCurrentFile()->GetName() ==  nonswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside numu nonswap\n";
            numutonumu->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG), OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          }
          else if (ch1->GetCurrentFile()->GetName() == nueswap) {
            //std::cerr << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG;
            std::cerr << "This should not happen (numu in nueswap) \n";
          }
          else if (ch1->GetCurrentFile()->GetName() == tauswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside numu tauswap\n";
            nuetonumu->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG), OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          }
          else {
            std::cerr << "This should not happen (file name not recognized) \n";
          }

        } else if (nuPDG == 16 || nuPDG == -16) {
          wcvnnue_fromnutau->Fill(cvnnue, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnumu_fromnutau->Fill(cvnnumu, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnutau_fromnutau->Fill(cvnnutau, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          wcvnnc_fromnutau->Fill(cvnnc, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          weDepPimNutau->Fill(eDepPim, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          weDepElseNutau->Fill(eDepHad - eDepPim, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hDepPimRatioT->Fill(Ev, DepPimRatio, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hResolutionT->Fill((Ev - eDepTotal)/Ev, Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hTruevRecoET->Fill(Ev, eDepTotal, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          hResolutionT2->Fill(eDepTotal/Ev, Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          if (3.5 <= Ev && Ev <= 4.5)
            hResolutionT1D->Fill(eDepTotal/Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));

          if (ch1->GetCurrentFile()->GetName() == nonswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside nutau nonswap\n";
            std::cerr << "This should not happen (nutau in nonswap) \n";
          }
          else if (ch1->GetCurrentFile()->GetName() == nueswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside nutau nueswap\n";
            nuetonutau->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG), OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          }
          else if (ch1->GetCurrentFile()->GetName() == tauswap) {
            //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << nuPDG << "Inside nutau tauswap\n";
            numutonutau->Fill(Ev, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG), OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
          }
          else {
            std::cerr << "This should not happen (file name not recognized) \n";
          }
        } else {
          std::cerr << "This should not happen: (nuPDG not recognized) \n";
        }

      } else {
        //Neutral current events are the same in each swap file, divide by 3
        //std::cout << "Testing..." << ch1->GetCurrentFile()->GetName() << ' ' << isCC << "Inside neutral current\n";
        wcvnnue_fromnc->Fill(cvnnue, (1. / 3.));
        wcvnnumu_fromnc->Fill(cvnnumu, (1. / 3.));
        wcvnnutau_fromnc->Fill(cvnnutau, (1. / 3.));
        wcvnnc_fromnc->Fill(cvnnc, (1. / 3.));

        weDepPimNC->Fill(eDepPim, (1. / 3.));
        weDepElseNC->Fill(eDepPim, (1. / 3.));

        hDepPimRatioNC->Fill(Ev, DepPimRatio, OscWeight(ch1->GetCurrentFile()->GetName(), Ev, nuPDG));
      }
//    }
    if (i % 10000 == 0) std::cout << (float) (i + 1) * 100 / (float) nentries << '%' << '\n';
  }

  outputFile->WriteObject(wcvnnue_fromnue, "wcvnnue_fromnue");
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
  outputFile->WriteObject(numutonutau, "numutonutauhisto");

  c1->cd();

  wcvnnue_fromnue->SetLineWidth(2);
  wcvnnue_fromnue->SetTitle("FHC CVN nu_e scores from nu_e events");
  wcvnnue_fromnue->GetXaxis()->SetTitle("CVN score");
  wcvnnue_fromnue->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnue_fromnue.png");
  c1->Clear();
  c1->Update();

  wcvnnue_fromnumu->SetLineWidth(2);
  wcvnnue_fromnumu->SetTitle("FHC CVN nu_e scores from nu_mu events");
  wcvnnue_fromnumu->GetXaxis()->SetTitle("CVN score");
  wcvnnue_fromnumu->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnue_fromnumu.png");
  c1->Clear();
  c1->Update();

  wcvnnue_fromnutau->SetLineWidth(2);
  wcvnnue_fromnutau->SetTitle("FHC CVN nu_e scores from nu_tau events");
  wcvnnue_fromnutau->GetXaxis()->SetTitle("CVN score");
  wcvnnue_fromnutau->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnue_fromnutau.png");
  c1->Clear();
  c1->Update();

  wcvnnue_fromnc->SetLineWidth(2);
  wcvnnue_fromnc->SetTitle("FHC CVN nu_e scores from Neutral Current events");
  wcvnnue_fromnc->GetXaxis()->SetTitle("CVN score");
  wcvnnue_fromnc->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnue_fromnc.png");
  c1->Clear();
  c1->Update();

  wcvnnumu_fromnue->SetLineWidth(2);
  wcvnnumu_fromnue->SetTitle("FHC CVN nu_mu scores from nu_e events");
  wcvnnumu_fromnue->GetXaxis()->SetTitle("CVN score");
  wcvnnumu_fromnue->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnumu_fromnue.png");
  c1->Clear();
  c1->Update();

  wcvnnumu_fromnumu->SetLineWidth(2);
  wcvnnumu_fromnumu->SetTitle("FHC CVN nu_mu scores from nu_mu events");
  wcvnnumu_fromnumu->GetXaxis()->SetTitle("CVN score");
  wcvnnumu_fromnumu->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnumu_fromnumu.png");
  c1->Clear();
  c1->Update();

  wcvnnumu_fromnutau->SetLineWidth(2);
  wcvnnumu_fromnutau->SetTitle("FHC CVN nu_mu scores from nu_tau events");
  wcvnnumu_fromnutau->GetXaxis()->SetTitle("CVN score");
  wcvnnumu_fromnutau->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnumu_fromnutau.png");
  c1->Clear();
  c1->Update();

  wcvnnumu_fromnc->SetLineWidth(2);
  wcvnnumu_fromnc->SetTitle("FHC CVN nu_mu scores from Neutral Current events");
  wcvnnumu_fromnc->GetXaxis()->SetTitle("CVN score");
  wcvnnumu_fromnc->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnumu_fromnc.png");
  c1->Clear();
  c1->Update();

  wcvnnutau_fromnue->SetLineWidth(2);
  wcvnnutau_fromnue->SetTitle("FHC CVN nu_tau scores from nu_e events");
  wcvnnutau_fromnue->GetXaxis()->SetTitle("CVN score");
  wcvnnutau_fromnue->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnutau_fromnue.png");
  c1->Clear();
  c1->Update();

  wcvnnutau_fromnumu->SetLineWidth(2);
  wcvnnutau_fromnumu->SetTitle("FHC CVN nu_tau scores from nu_mu events");
  wcvnnutau_fromnumu->GetXaxis()->SetTitle("CVN score");
  wcvnnutau_fromnumu->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnutau_fromnumu.png");
  c1->Clear();
  c1->Update();

  wcvnnutau_fromnutau->SetLineWidth(2);
  wcvnnutau_fromnutau->SetTitle("FHC CVN nu_tau scores from nu_tau events");
  wcvnnutau_fromnutau->GetXaxis()->SetTitle("CVN score");
  wcvnnutau_fromnutau->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnutau_fromnutau.png");
  c1->Clear();
  c1->Update();

  wcvnnutau_fromnc->SetLineWidth(2);
  wcvnnutau_fromnc->SetTitle("FHC CVN nu_tau scores from Neutral Current events");
  wcvnnutau_fromnc->GetXaxis()->SetTitle("CVN score");
  wcvnnutau_fromnc->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnutau_fromnc.png");
  c1->Clear();
  c1->Update();

  wcvnnc_fromnue->SetLineWidth(2);
  wcvnnc_fromnue->SetTitle("FHC CVN Neutral Current scores from nu_e events");
  wcvnnc_fromnue->GetXaxis()->SetTitle("CVN score");
  wcvnnc_fromnue->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnc_fromnue.png");
  c1->Clear();
  c1->Update();

  wcvnnc_fromnumu->SetLineWidth(2);
  wcvnnc_fromnumu->SetTitle("FHC CVN Neutral Current scores from nu_mu events");
  wcvnnc_fromnumu->GetXaxis()->SetTitle("CVN score");
  wcvnnc_fromnumu->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnc_fromnumu.png");
  c1->Clear();
  c1->Update();

  wcvnnc_fromnutau->SetLineWidth(2);
  wcvnnc_fromnutau->SetTitle("FHC CVN Neutral Current scores from nu_tau events");
  wcvnnc_fromnutau->GetXaxis()->SetTitle("CVN score");
  wcvnnc_fromnutau->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnc_fromnutau.png");
  c1->Clear();
  c1->Update();

  wcvnnc_fromnc->SetLineWidth(2);
  wcvnnc_fromnc->SetTitle("FHC CVN Neutral Current scores from Neutral Current events");
  wcvnnc_fromnc->GetXaxis()->SetTitle("CVN score");
  wcvnnc_fromnc->Draw("HIST");
  c1->SaveAs("./Figures/wcvnnc_fromnc.png");
  c1->Clear();
  c1->Update();

  weDepPimNue->SetLineWidth(2);
  weDepPimNue->SetTitle("FHC #pi#^- Deposited Energy Dist. from #nu_e");
  weDepPimNue->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepPimNue->Draw("HIST");
  c1->SaveAs("./Figures/weDepPimNue.png");
  c1->Clear();
  c1->Update();

  weDepPimNumu->SetLineWidth(2);
  weDepPimNumu->SetTitle("FHC #pi#^- Deposited Energy Dist. from #nu_#mu");
  weDepPimNumu->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepPimNumu->Draw("HIST");
  c1->SaveAs("./Figures/weDepPimNumu.png");
  c1->Clear();
  c1->Update();

  weDepPimNutau->SetLineWidth(2);
  weDepPimNutau->SetTitle("FHC #pi#^- Deposited Energy Dist. from #nu_#tau");
  weDepPimNutau->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepPimNutau->Draw("HIST");
  c1->SaveAs("./Figures/weDepPimNutau.png");
  c1->Clear();
  c1->Update();

  weDepPimNC->SetLineWidth(2);
  weDepPimNC->SetTitle("FHC #pi#^- Deposited Energy Dist. from NC events");
  weDepPimNC->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepPimNC->Draw("HIST");
  c1->SaveAs("./Figures/weDepPimNC.png");
  c1->Clear();
  c1->Update();

  weDepElseNue->SetLineWidth(2);
  weDepElseNue->SetTitle("FHC Remaining Deposited Energy Dist. from #nu_e");
  weDepElseNue->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepElseNue->Draw("HIST");
  c1->SaveAs("./Figures/weDepElseNue.png");
  c1->Clear();
  c1->Update();

  weDepElseNumu->SetLineWidth(2);
  weDepElseNumu->SetTitle("FHC Remaining Deposited Energy Dist. from #nu_#mu");
  weDepElseNumu->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepElseNumu->Draw("HIST");
  c1->SaveAs("./Figures/weDepElseNumu.png");
  c1->Clear();
  c1->Update();

  weDepElseNutau->SetLineWidth(2);
  weDepElseNutau->SetTitle("FHC Remaining Deposited Energy Dist. from #nu_#tau");
  weDepElseNutau->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepElseNutau->Draw("HIST");
  c1->SaveAs("./Figures/weDepElseNutau.png");
  c1->Clear();
  c1->Update();

  weDepElseNC->SetLineWidth(2);
  weDepElseNC->SetTitle("FHC Remaining Deposited Energy Dist. from NC events");
  weDepElseNC->GetXaxis()->SetTitle("Dep. E (GeV)");
  weDepElseNC->Draw("HIST");
  c1->SaveAs("./Figures/weDepElseNC.png");
  c1->Clear();
  c1->Update();

  nuetonue->SetLineWidth(2);
  nuetonue->Draw("COLZ");
  nuetonue->SetTitle("P(nu_e -> nu_e) events");
  nuetonue->GetXaxis()->SetTitle("True Energy (GeV)");
  nuetonue->GetYaxis()->SetTitle("Oscillation Probability");
  c1->SaveAs("./Figures/nuetonuehisto.png");
  c1->Clear();
  c1->Update();

  nuetonumu->SetLineWidth(2);
  nuetonumu->Draw("COLZ");
  nuetonumu->SetTitle("P(nu_e -> nu_mu) events");
  nuetonumu->GetXaxis()->SetTitle("True Energy (GeV)");
  nuetonumu->GetYaxis()->SetTitle("Oscillation Probability");
  c1->SaveAs("./Figures/nuetonumuhisto.png");
  c1->Clear();
  c1->Update();

  nuetonutau->SetLineWidth(2);
  nuetonutau->Draw("COLZ");
  nuetonutau->SetTitle("P(nu_e -> nu_tau) events");
  nuetonutau->GetXaxis()->SetTitle("True Energy (GeV)");
  nuetonutau->GetYaxis()->SetTitle("Oscillation Probability");
  c1->SaveAs("./Figures/nuetonutauhisto.png");
  c1->Clear();
  c1->Update();

  numutonue->SetLineWidth(2);
  numutonue->Draw("COLZ");
  numutonue->SetTitle("P(nu_mu -> nu_e) events");
  numutonue->GetXaxis()->SetTitle("True Energy (GeV)");
  numutonue->GetYaxis()->SetTitle("Oscillation Probability");
  c1->SaveAs("./Figures/numutonuehisto.png");
  c1->Clear();
  c1->Update();

  numutonumu->SetLineWidth(2);
  numutonumu->Draw("COLZ");
  numutonumu->SetTitle("P(nu_mu -> nu_mu) events");
  numutonumu->GetXaxis()->SetTitle("True Energy (GeV)");
  numutonumu->GetYaxis()->SetTitle("Oscillation Probability");
  c1->SaveAs("./Figures/numutonumuhisto.png");
  c1->Clear();
  c1->Update();

  numutonutau->SetLineWidth(2);
  numutonutau->Draw("COLZ");
  numutonutau->SetTitle("P(nu_mu -> nu_tau) events");
  numutonutau->GetXaxis()->SetTitle("True Energy (GeV)");
  numutonutau->GetYaxis()->SetTitle("Oscillation Probability");
  c1->SaveAs("./Figures/numutonutauhisto.png");
  c1->Clear();
  c1->Update();

  hDepPimRatioE->Draw("COLZ");
  hDepPimRatioE->SetTitle("FHC E_{#pi^{-}} / Dep. E for #nu_{e} events");
  hDepPimRatioE->GetXaxis()->SetTitle("True Energy (GeV)");
  hDepPimRatioE->GetYaxis()->SetTitle("E_{#pi^{-}} / Dep. E");
  c1->SaveAs("./Figures/DepPimRatioE.png");
  c1->Clear();
  c1->Update();

  hDepPimRatioM->Draw("COLZ");
  hDepPimRatioM->SetTitle("FHC E_{#pi^{-}} / Dep. E for #nu_{#mu} events");
  hDepPimRatioM->GetXaxis()->SetTitle("True Energy (GeV)");
  hDepPimRatioM->GetYaxis()->SetTitle("E_{#pi^{-}} / Dep. E");
  c1->SaveAs("./Figures/DepPimRatioM.png");
  c1->Clear();
  c1->Update();

  hDepPimRatioT->Draw("COLZ");
  hDepPimRatioT->SetTitle("FHC E_{#pi^{-}} / Dep. E for #nu_{#tau} events");
  hDepPimRatioT->GetXaxis()->SetTitle("True Energy (GeV)");
  hDepPimRatioT->GetYaxis()->SetTitle("E_{#pi^{-}} / Dep. E");
  c1->SaveAs("./Figures/DepPimRatioT.png");
  c1->Clear();
  c1->Update();

  hDepPimRatioNC->Draw("COLZ");
  hDepPimRatioNC->SetTitle("FHC E_{#pi^{-}} / Dep. E for NC events");
  hDepPimRatioNC->GetXaxis()->SetTitle("True Energy (GeV)");
  hDepPimRatioNC->GetYaxis()->SetTitle("E_{#pi^{-}} / Dep. E");
  c1->SaveAs("./Figures/DepPimRatioNC.png");
  c1->Clear();
  c1->Update();

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

  hResolutionT->Draw("COLZ");
  hResolutionT->SetTitle("FHC Fractional Energy Resolution for #nu_{#tau}");
  hResolutionT->GetXaxis()->SetTitle("(True E - Reco. E) / True E");
  hResolutionT->GetYaxis()->SetTitle("True E (GeV)");
  c1->SaveAs("./Figures/ResolutionT.png");
  c1->Clear();
  c1->Update();

  hTruevRecoEE->Draw("COLZ");
  hTruevRecoEE->SetTitle("FHC #nu_{e} energy");
  hTruevRecoEE->GetXaxis()->SetTitle("True E (GeV)");
  hTruevRecoEE->GetYaxis()->SetTitle("Reco. E (GeV)");
  c1->SaveAs("./Figures/TruevRecoEVE.png");
  c1->Clear();
  c1->Update();

  hTruevRecoEM->Draw("COLZ");
  hTruevRecoEM->SetTitle("FHC #nu_{#mu} energy");
  hTruevRecoEM->GetXaxis()->SetTitle("True E (GeV)");
  hTruevRecoEM->GetYaxis()->SetTitle("Reco. E (GeV)");
  c1->SaveAs("./Figures/TruevRecoEVM.png");
  c1->Clear();
  c1->Update();

  hTruevRecoET->Draw("COLZ");
  hTruevRecoET->SetTitle("FHC #nu_{#tau} energy");
  hTruevRecoET->GetXaxis()->SetTitle("True E (GeV)");
  hTruevRecoET->GetYaxis()->SetTitle("Reco. E (GeV)");
  c1->SaveAs("./Figures/TruevRecoEVT.png");
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

  hResolutionT2->Draw("COLZ");
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

  hResolutionT1D->SetLineWidth(2);
  hResolutionT1D->SetTitle("FHC Energy Resolution for #nu_{#tau} (3.5 GeV < Ev < 4.5 GeV)");
  hResolutionT1D->GetXaxis()->SetTitle("Reco. E / True E");
  hResolutionT1D->Draw("HIST");
  c1->SaveAs("./Figures/ResolutionT1D.png");
  c1->Clear();
  c1->Update();



  outputFile->Close();

}
