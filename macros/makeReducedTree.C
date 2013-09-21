//--------------------------------------------------------------------------------------------------
// Reduce the ntuples and compute all the necessary weights for further studies and plots
//
// Authors: L. Di Matteo                                                                  (Aug 2013)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TFile.h"
#include "TSystem.h"
#include "TSystem.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitMonoPhoton/Utils/interface/TreeReducer.h"
#endif

using namespace std;
using namespace mithep;

TString getEnv(const char* name);

//==================================================================================================
void makeReducedTree()
{
  // read all environment variables
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_ANA_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((mitHgg + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < *samples->NDataSamples(); iSample++) listOfSamples.push_back(samples->GetDataSample(iSample));
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  
  
  // define outfile
  TString outfileName = prdCfg + "_reduced_new.root";
  TFile* outfile = new TFile(outfileName,"RECREATE");
  // define infolder
  TString sampleBaseDir = *samples->Dir();
  std::cout << "sampleBaseDir " << sampleBaseDir << std::endl;
  // define normalized target PU distribution
  TFile* pufile   = new TFile("PileUpHistograms/MyDataPileupHistogram.root","READ");
  TH1D*  putarget = (TH1D*) pufile -> Get("pileup");
  putarget -> Scale(1.0/putarget->GetSumOfWeights());
  TFile* pufileup   = new TFile("PileUpHistograms/MyDataPileupHistogramUp.root","READ");
  TH1D*  putargetup = (TH1D*) pufileup -> Get("pileup");
  putargetup -> Scale(1.0/putargetup->GetSumOfWeights());
  TFile* pufiledown   = new TFile("PileUpHistograms/MyDataPileupHistogramDown.root","READ");
  TH1D*  putargetdown = (TH1D*) pufiledown -> Get("pileup");
  putargetdown -> Scale(1.0/putargetdown->GetSumOfWeights());

  // generate reduced trees
  // loop through the samples and produce the reduced trees
  for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {
    TreeReducer  *thisReducer = new TreeReducer(listOfSamples.at(iSample));
    thisReducer -> SetVerbose(true);
    thisReducer -> SetPUTarget(putarget);
    thisReducer -> SetPUTargetUp(putargetup);
    thisReducer -> SetPUTargetDown(putargetdown);
    thisReducer -> SetInputBaseDir(sampleBaseDir);
    thisReducer -> SetOutput(outfile);
    thisReducer -> SetLumi(19789.);
    thisReducer -> MakeTree();
    delete thisReducer;
    
  }//end loop on samples
    
  outfile   -> Close();
  
  return;
}

TString getEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}
