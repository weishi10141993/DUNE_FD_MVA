/*
   Implementation of machine learning algorithms to optimize cutting on nutau signal/bkgd
*/


#include <iostream>
#include <cstdlib>

#include "TMVA/Types.h"


void MVA(void) {

   bool useRHC = false;

   TString FHCSplitFile = "/storage/shared/ncchambe/FDMonteCarlo/FHCSplitMC.root";
   TString RHCSplitFile = "/storage/shared/ncchambe/FDMonteCarlo/RHCSplitMC.root";

   TFile* inputFile = NULL;
   if (!useRHC) {
      inputFile = new TFile(FHCSplitFile, "READ");
   }
   else {
      inputFile = new TFile(RHCSplitFile, "READ");
   }

   TTree* NuTauSignalFromTauswap = (TTree*)inputFile->Get("NuTauSignalFromTauswap");
   TTree* NuTauBkgdFromTauswap = (TTree*)inputFile->Get("NuTauBkgdFromTauswap");
   TTree* NuTauFromNueswap = (TTree*)inputFile->Get("NuTauFromNueswap");
   TTree* NueFromNonswap = (TTree*)inputFile->Get("NueFromNonswap");
   TTree* NueFromNueswap = (TTree*)inputFile->Get("NueFromNueswap");
   TTree* NumuFromNonswap = (TTree*)inputFile->Get("NumuFromNonswap");
   TTree* NumuFromTauswap = (TTree*)inputFile->Get("NumuFromTauswap");
   TTree* NCTree = (TTree*)inputFile->Get("NCTree");



   std::vector<TTree*> TreeArray = {NuTauSignalFromTauswap, NuTauBkgdFromTauswap, NuTauFromNueswap, NueFromNonswap, NueFromNueswap,
      NumuFromNonswap, NumuFromTauswap, NCTree};

   TString FHCOutFileName = "/storage/shared/ncchambe/FDMonteCarlo/MVA/FHCMVA.root";
   TString RHCOutFileName = "/storage/shared/ncchambe/FDMonteCarlo/MVA/RHCMVA.root";

   TFile* outputFile = NULL;
   if (!useRHC) {
      outputFile = new TFile(FHCOutFileName, "RECREATE");
   }
   else {
      outputFile = new TFile(RHCOutFileName, "RECREATE");
   }

   TMVA::Factory* factory = NULL;
   if (!useRHC) {
      factory = new TMVA::Factory("FHCNuTauMVA", outputFile,
         "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
   }
   else {
      factory = new TMVA::Factory("RHCNuTauMVA", outputFile,
         "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
   }

   TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

   for (auto i : TreeArray) {
      i->SetBranchStatus("*", false);
      i->SetBranchStatus("cvnnue", true);
      i->SetBranchStatus("cvnnumu", true);
      i->SetBranchStatus("cvnnutau", true);
      i->SetBranchStatus("POTScaledOscweight", true);
      i->SetBranchStatus("isCC", true);
      i->SetBranchStatus("eDepPim", true);
      i->SetBranchStatus("eDepP", true);
      i->SetBranchStatus("eDepN", true);
      i->SetBranchStatus("eDepPip", true);
      i->SetBranchStatus("eDepPi0", true);
      i->SetBranchStatus("eDepOther", true);
      i->SetBranchStatus("nipim", true);
      i->SetBranchStatus("eRecoPim", true);
      i->SetBranchStatus("eRecoP", true);
      i->SetBranchStatus("eRecoN", true);
      i->SetBranchStatus("eRecoPip",true);
      i->SetBranchStatus("eRecoPi0", true);
      i->SetBranchStatus("eRecoOther", true);
      i->SetBranchStatus("Ev", true);
      i->SetBranchStatus("nuPDG", true);
      i->SetBranchStatus("cvnnc", true);
      i->SetBranchStatus("nik0", true);
      i->SetBranchStatus("nikp", true);
      i->SetBranchStatus("nikm", true);
      i->SetBranchStatus("niother", true);
      i->SetBranchStatus("nipip", true);
      i->SetBranchStatus("niem", true);
      i->SetBranchStatus("RecoLepEnNue", true);
      i->SetBranchStatus("RecoLepEnNumu", true);
      i->SetBranchStatus("vtx_x", true);
      i->SetBranchStatus("vtx_y", true);
      i->SetBranchStatus("vtx_z", true);
   }

   dataloader->SetWeightExpression("POTScaledOscweight");

   double POTYears = 3.5;


   dataloader->AddSignalTree(NuTauSignalFromTauswap, POTYears);
   dataloader->AddBackgroundTree(NuTauBkgdFromTauswap, POTYears);
   dataloader->AddBackgroundTree(NuTauFromNueswap, POTYears);
   dataloader->AddBackgroundTree(NumuFromTauswap, POTYears);
   dataloader->AddBackgroundTree(NumuFromNonswap, POTYears);
   dataloader->AddBackgroundTree(NueFromNueswap, POTYears);
   dataloader->AddBackgroundTree(NueFromNonswap, POTYears);
   dataloader->AddBackgroundTree(NCTree, POTYears);

   //For some reason, the variables are converted to floats and must be accessed as floats in MVAStats.C
   dataloader->AddVariable("cvnnue", 'F');
   dataloader->AddVariable("cvnnumu", 'F');
   dataloader->AddVariable("cvnnutau", 'F');
   dataloader->AddVariable("EHad := eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther", 'F');
   dataloader->AddVariable("eDepPim", 'F');
   dataloader->AddVariable("eDepPip", 'F');
   dataloader->AddVariable("nipim", 'I');
   dataloader->AddVariable("cvnnc", 'F');
   dataloader->AddVariable("nipip", 'I');

   dataloader->AddSpectator("Ev", 'F');
   dataloader->AddSpectator("POTScaledOscweight", 'F');

   //Remove events with anomalous CVN scores and apply fiducial volume cut
   TCut cut = "(cvnnue > 0.) && (cvnnumu > 0.) && (cvnnutau > 0.) && (cvnnc > 0.) &&"
      "(vtx_x > -310) && (vtx_x < 310) && (vtx_y > -550) && (vtx_y < 550) && (vtx_z > 50) && (vtx_z < 1244)";


   dataloader->PrepareTrainingAndTestTree(cut, cut, "");

   /*
      Various MVA methods. Adaptive BDT does the best. I only tested BDT methods and 'Cuts' method. Other methods were too expensive to run on NNGroup
      machine
   */

   //factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts",
   //      "!H:!V:FitMethod=MC:EffSel:VarProp=FSmart");
   //factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERS",
   //      "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");
   //factory->BookMethod(dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.1:VarTransform=Norm");
   //factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
   //   "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");
   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTA",
      "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
   //factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTB",
   //   "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");
   //factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
   //   "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");
   //factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTF",
   //   "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");
   /*factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU",
         "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:"
         "Layout=TANH|128,TANH|128,TANH|128,LINEAR:"
         "TrainingStrategy=LearningRate=1e-2,Momentum=0.9,"
         "ConvergenceSteps=20,BatchSize=100,TestRepetitions=1,"
         "WeightDecay=1e-4,Regularization=None,"
         "DropConfig=0.0+0.5+0.5+0.5"
         ":Architecture=CPU");*/


   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   outputFile->Close();
}
