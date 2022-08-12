# FDMonteCarlo
Stony Brook REU Project

			Running the BDT for Nutau selection

Macros:
SplitMCTree.C		--Organizes FD CAF files into signal/background types for training

MVA.C			--Trains the BDT with the files (FHC/RHC)SplitMC.root. Outputs XML file containing the weights used in the decision tree boosting. 
  Also outputs (FHC/RHC)MVA.root containing graphs depicting training results

MVASignalEval.C		--Implements BDT training for nutau selection and evaluates significance/purity/efficiency. Creates stacked histogram of events passing the cut. 
  Input is XML file mentioned previously.

MVAChi2Test.C		--Uses the same structure as MVASignalEval.C, but creates the Chi2 graphs as another metric of signal significance.

CafCutEval.C		--Same analysis as in MVASignalEval.C, but configured with the naive cut using CVN variables. Also outputs stacked histogram and 
significance/purity/efficiency

CafChi2Test.C		--Same analysis as MVAChi2Test.C, but configured with naive cut using CVN variables.

CafCutStatPlot.C  --This file is what I used to create the various plots showing the distribution of signal/background events according to training variables I used.
  Additionally, I used this file to create 2D histograms of Reco/True Energy to study reconstruction for nutau as well as plot the distribution of CVN scores for 
  each type of event.


I/O files and data:
FD_(FHC/RHC)_(non/nue/tau)swap.root -- Raw FD MC data

(FHC/RHC)Split.root	--Output from SplitMCTree.C. The raw FD CAF files must be organized into signal/background during the TMVA training. 
  These files have several trees comprised of the total events summed up from non/nue/tauswap.root. Additionally, oscillation weight is 
  implemented into these files, which is not present in FD CAF files

(FHC/RHC)NuTauMVA.weights.xml --Output from MVA.C containing BDT training information necessary for implementation



Instructions:
If starting from scratch with FD MC files (non/nue/tau)swap, begin the analysis by running SplitMCTree.C. In order to train the BDT using the TMVA class, all events
must be separated into signal/background trees. The SplitMCTree file takes all 3 FD Caf files and combines the events into a TChain before separating the events into 
a handful of trees according to their event type (CC numu/nue/nutau, NC, CC nutau_bar, etc.). The output of this macro (FHC/RHC)split.root is read into the BDT
training. It should be noted that I have computed the oscillation probability for each event as well as scaled each event according to the POT scaling for 1 year of
data taking. To change the number of years of data taking, one must adjust the variable 'POTYears' in MVA.C as well as the variable 'scalefactor' in the analysis 
files.

In each of the files where it is relevant, I have attempted to implement the option of running FHC/RHC in the same file. When I was performing the analysis, I just
duplicates of the same type of file for FHC/RHC. You'll find 'bool useRHC' at the beginning of each main function, so set this to 'true' and it should perform the 
analysis with RHC Caf files. I have not tested the 'useRHC' implementation in each of the files, so there may be issues there.

After running SplitMCTree.C, the next step is to train the BDT in MVA.C. Events are trained with a weight corresponding to their oscillation probability (stored in the
POTScaledOscWeight branch in (FHC/RHC)Split.root). At the bottom of the MVA.C file, you will find where I have commented out other MVA algorithms. The BDT algorithm
with adaptive boosting (BDTA) performed the best of all the BDT algorithms. Other algorithms such as PDE-RS and SVM took too much memory and could not be run on NNHome
in compute nodes in either batch mode or interactive mode. BDT training takes ~20 minutes. 

Save MVA.C output to a log file. In the output, you will find towards the beginning a stat labelled 'Sum of weights:' for both signal and background. These are, of 
course, the weighted sum of signal/background events. See the screenshot SumOfWeights.png in the git repository. These numbers are what you put into the BDT root file
to extract the optimal BDT cut obtained by the algorithm. There's also a command in the TMVA class that does this, I believe. After the training is complete, run

> root -l -e 'TMVA::TMVAGui("[FHC/RHC]MVA.root")'

to open the training stats output by the TMVA class for this BDT. Click the tab labeled 'Classifier cut efficiencies' and a graph displaying significance, efficiency,
and purity as a function of the BDT cut will appear. Another box will appear in the corner with 'Signal Events' and 'Background events'. Input the 'sum of weights' 
numbers into those fields and click draw. The graphs will change and the optimal BDT cut will appear in the bottom left corner. This is what I put into the subsequent
analysis files.

After this, use MVASignalEval.C to obtain the efficiency/purity/significance of the BDT cut. You'll notice in loading the BDT output information, we include the 
weights from the XML file output by the training. Additionally, TMVA requires that all of the placeholder variables used to evaluate the BDT response of each event 
must be floats. In the for loop over the events, we typecast all of the doubles to floats before using these variables to evaluate the BDT response. 

Finally, MVAChi2Test.C generates the Chi2 graphs and the stacked histogram of events passing the BDT cut. In this file, I included the option of binning according to 
Reco E_Had. I should note that you should double check the legend labels and saved file names. 

The files CafCutEval.C, CafChi2Test.C, CafCutStatPlot.C perform similar analysis on the naive CVN cut, if you're interested.



