/*
   Implementation of machine learning algorithms to optimize cutting on nutau signal/bkgd
*/


#include <iostream>
#include <cstdlib>

#include "TMVA/Types.h"

TString nonswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nonswap.root";
TString nueswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_nueswap.root";
TString tauswap = "/storage/shared/wshi/CAFs/FD/FD_FHC_tauswap.root";



void NNTest(void) {
   TFile* inputFile = new TFile("/storage/shared/ncchambe/FDMonteCarlo/FHCSplitMC.root", "READ");
   TTree* NuTauFromTauswap = (TTree*)inputFile->Get("NuTauFromTauswap");
   TTree* NuTauFromNueswap = (TTree*)inputFile->Get("NuTauFromNueswap");
   TTree* NueFromNonswap = (TTree*)inputFile->Get("NueFromNonswap");
   TTree* NueFromNueswap = (TTree*)inputFile->Get("NueFromNueswap");
   TTree* NumuFromNonswap = (TTree*)inputFile->Get("NumuFromNonswap");
   TTree* NumuFromTauswap = (TTree*)inputFile->Get("NumuFromTauswap");
   TTree* NCTree = (TTree*)inputFile->Get("NCTree");

   /*TTree* TestNuTauFromTauswap = (TTree*)inputFile->Get("TestNuTauFromTauswap");
   TTree* TestNuTauFromNueswap = (TTree*)inputFile->Get("TestNuTauFromNueswap");
   TTree* TestNueFromNonswap = (TTree*)inputFile->Get("TestNueFromNonswap");
   TTree* TestNueFromNueswap = (TTree*)inputFile->Get("TestNueFromNueswap");
   TTree* TestNumuFromNonswap = (TTree*)inputFile->Get("TestNumuFromNonswap");
   TTree* TestNumuFromTauswap = (TTree*)inputFile->Get("TestNumuFromTauswap");
   TTree* TestNCTree = (TTree*)inputFile->Get("TestNCTree");*/

   std::vector<TTree*> TreeArray = {NuTauFromTauswap, NuTauFromNueswap, NueFromNonswap, NueFromNueswap,
      NumuFromNonswap, NumuFromTauswap, NCTree,
   //    TestNuTauFromTauswap, TestNuTauFromNueswap, TestNueFromNonswap, TestNueFromNueswap,
   //      TestNumuFromNonswap, TestNumuFromTauswap, TestNCTree
    };

   TFile* outputFile = new TFile("/storage/shared/ncchambe/FDMonteCarlo/NNTest2.root", "RECREATE");

   TMVA::Factory* factory = new TMVA::Factory("NuTauClassification", outputFile,
         "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

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
   }

   dataloader->SetWeightExpression("POTScaledOscweight");

   double signalweight = 3.5;
   double bkgdweight = 3.5; //These variables are the years of running the exp


   dataloader->AddSignalTree(NuTauFromTauswap, signalweight
   //    ,TMVA::Types::kTraining
   );
   dataloader->AddBackgroundTree(NuTauFromNueswap, bkgdweight
       //   ,TMVA::Types::kTraining
   );
   dataloader->AddBackgroundTree(NumuFromTauswap, bkgdweight
       //   ,TMVA::Types::kTraining
   );
   dataloader->AddBackgroundTree(NumuFromNonswap, bkgdweight
       //   ,TMVA::Types::kTraining
   );
   dataloader->AddBackgroundTree(NueFromNueswap, bkgdweight
       //   ,TMVA::Types::kTraining
   );
   dataloader->AddBackgroundTree(NueFromNonswap, bkgdweight
       //   ,TMVA::Types::kTraining
   );
   dataloader->AddBackgroundTree(NCTree, bkgdweight
       //   ,TMVA::Types::kTraining
   );

   //Optional tree for testing
   /*dataloader->AddSignalTree(TestNuTauFromTauswap, signalweight, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(TestNuTauFromNueswap, bkgdweight, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(TestNumuFromTauswap, bkgdweight, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(TestNumuFromNonswap, bkgdweight, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(TestNueFromNueswap, bkgdweight, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(TestNueFromNonswap, bkgdweight, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(TestNCTree, bkgdweight, TMVA::Types::kTesting);*/

   dataloader->AddVariable("cvnnue", 'F');
   dataloader->AddVariable("cvnnumu", 'F');
   dataloader->AddVariable("cvnnutau", 'F');
   dataloader->AddVariable("EPimRatio := eDepPim / (eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther)", 'F');
   dataloader->AddVariable("nipim", 'I');
   dataloader->AddVariable("cvnnc", 'F');

   dataloader->AddSpectator("RecoEVis := eRecoPim + eRecoPip + eRecoPi0 + eRecoP + eRecoN + eRecoOther", 'F');
   dataloader->AddSpectator("Ev", 'F');
   dataloader->AddSpectator("POTScaledOscweight", 'F');

   TCut mycut = "(cvnnue > 0.) && (cvnnumu > 0.) && (cvnnutau > 0.) && (cvnnc > 0.) && (eDepP + eDepN + eDepPip + eDepPim + eDepPi0 + eDepOther <= 10.)";

   dataloader->PrepareTrainingAndTestTree(mycut, mycut, "");

   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "");
   //factory->BookMethod(dataloader, TMVA::Types::kSVM, "SVM", "");
   /*factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU",
         "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:"
         "Layout=TANH|128,TANH|128,TANH|128,LINEAR:"
         "TrainingStrategy=LearningRate=1e-2,Momentum=0.9,"
         "ConvergenceSteps=20,BatchSize=100,TestRepetitions=1,"
         "WeightDecay=1e-4,Regularization=None,"
         "DropConfig=0.0+0.5+0.5+0.5"
         ":Architecture=CPU");*/
   //factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts",
   //      "!H:!V:FitMethod=MC:EffSel:VarProp=FSmart");
   //factory->BookMethod(dataloader, TMVA::Types::kPDERS, "PDERS",
   //      "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   outputFile->Close();
}
