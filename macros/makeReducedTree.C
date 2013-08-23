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
  
  // define outfile
  TFile* outfile = new TFile("out.root","RECREATE");
  // define infolder
  TString sampleBaseDir = *samples->Dir();
  std::cout << "sampleBaseDir " << sampleBaseDir << std::endl;
  // define normalized target PU distribution
  TFile* pufile   = new TFile("MyDataPileupHistogram.root","READ");
  TH1D*  putarget = (TH1D*) pufile -> Get("pileup");
  putarget -> Scale(1.0/putarget->GetSumOfWeights());

  // generate reduced trees
  // loop through the samples and produce the reduced trees
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) {
    const Sample *thisSample = samples->GetSample(iSample);    
    TreeReducer  *thisReducer = new TreeReducer(thisSample);
    thisReducer -> SetPUTarget(putarget);
    thisReducer -> SetInputBaseDir(sampleBaseDir);
    thisReducer -> SetOutput(outfile);
    thisReducer -> SetLumi(19500.);
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
