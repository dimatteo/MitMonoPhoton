//--------------------------------------------------------------------------------------------------
// Compute the systematic uncertainties from the reduced trees
//
// Authors: L. Di Matteo                                                                  (Sep 2013)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TFile.h"
#include "TSystem.h"
#include "TSystem.h"
#include <TTree.h>
#include <THStack.h>
#include <TLegend.h>
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitMonoPhoton/Core/MitGPTreeReduced.h"
#endif

using namespace std;
using namespace mithep;

TString getEnv(const char* name);

//==================================================================================================
void makeSyst()
{
  // read all environment variables
  TString here   = getEnv("PWD");
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_PLOT_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((mitHgg + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  
  vector<const Sample*> listOfDatasets;
  for (UInt_t iSample=0; iSample < *samples->NDataSamples(); iSample++) listOfDatasets.push_back(samples->GetDataSample(iSample));
  // define infolder
  TString inFileName = here + "/monoph-2013-July9_reduced.root";
  std::cout << "inFileName " << inFileName << std::endl;

  // prepare the systematic names : keep an eye to the ordering !
  vector<TString>  v_systNames;
  v_systNames.push_back("Nominal");
  v_systNames.push_back("PuUp");
  v_systNames.push_back("PuDown");
  const int nSyst = v_systNames.size();

  // prepare the systematic matrix to store all the yields
  vector<vector<float> > m_systYields; 

  // prepare the shortened sample list to organize the output in an easy way
  vector<TString>  v_groupSampleName;
  
  // loop through the samples and compute the yields for different conditions
  int   theGroupSampleCounter = 0;
  vector<float> v_thisYield;
  //for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {
  for (UInt_t iSample=0; iSample < 10; iSample++) {

    //Say which sample we are processing
    cout << "Processing sample " << *(listOfSamples.at(iSample)->Name()) << endl;

    //determine if the sample is first of the series
    bool isFirstOfSerie = (*listOfSamples.at(iSample)->Legend()).CompareTo(" ");
    bool isLastOfSerie = false;
    if (iSample == listOfSamples.size() - 1) isLastOfSerie = true;
    if (iSample < listOfSamples.size() - 1 && (*listOfSamples.at(iSample+1)->Legend()).CompareTo(" ") != 0) isLastOfSerie = true;
    
    //if sample first of the list insert the shortname in the appropriate vector and refresh the yield
    if (isFirstOfSerie) {
      v_groupSampleName.push_back(*(listOfSamples.at(iSample)->Name()));
      for (int iSyst=0; iSyst < nSyst; iSyst++) 
        v_thisYield.push_back(0.);
    }

    //get the sample scale 
    float thisScale = *(listOfSamples.at(iSample)->Scale());
 
    //get a nice GPTree for this sample 
    MitGPTreeReduced thistree;
    thistree.LoadTree(inFileName,*(listOfSamples.at(iSample)->Name()));
    thistree.InitTree();

    // Loop on the tree and make all the computations
    Int_t nEntries = thistree.tree_->GetEntries();
    for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
      thistree.tree_-> GetEntry(iEntry);
      //nominal
      v_thisYield[0] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_;
      //PuUp
      v_thisYield[1] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weightup_;
      //PuDown
      v_thisYield[2] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weightdo_;
    }
    
    cout << "this yield is " <<  v_thisYield[0] << " " << v_thisYield[1] << " " << v_thisYield[2] << endl;
    
    //add the yields to the syst matrix if this sample is last of the series:
    //either last sample or ~ sample followed by non ~ sample
    if (isLastOfSerie) {
       //Add me to the matrix!
       //clear the yield vector
       m_systYields.push_back(v_thisYield);
       v_thisYield.clear();
       theGroupSampleCounter++;
    }
    
  }//end loop on samples

  //test the matrix
  for (int iGroupSample = 0; iGroupSample < m_systYields.size(); iGroupSample++)
    cout << "the matrix yield is " << m_systYields.at(iGroupSample)[0] << " " << m_systYields.at(iGroupSample)[1] << " " << m_systYields.at(iGroupSample)[2] << endl;
        
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
