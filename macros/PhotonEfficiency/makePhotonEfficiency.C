//--------------------------------------------------------------------------------------------------
// Compute photon efficiencies in Zmumu + gamma events
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
#include <TRandom3.h>
#include <TVector2.h>
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitMonoPhoton/Core/MitGPTree.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
#include "TNtuple.h"
#endif

using namespace std;
using namespace mithep;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

TString getEnv(const char* name);

//==================================================================================================
void makePhotonEfficiency()
{
  // read all environment variables
  TString here   = getEnv("PWD");
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_EFF_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");
  
  // define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((mitHgg + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  
  vector<const Sample*> listOfDatasets;
  for (UInt_t iSample=0; iSample < *samples->NDataSamples(); iSample++) listOfDatasets.push_back(samples->GetDataSample(iSample));
  
  //define a very simple ntuple
  TNtuple* ZTuple = new TNtuple("ZTuple","Ntuple for Photon Efficiency Studies",
  "run:lumi:event:ZMass:lepid:ZEta:ZPhi:ZPt:ZEnergy:GEta:GPhi:GPt:GEnergy");

  //define an output file
  TString outFileName  = "ZTuple.root";
    
  for (UInt_t iDataset=0; iDataset < listOfDatasets.size(); iDataset++) {

    //Say which sample we are processing
    cout << "Processing sample " << *(listOfDatasets.at(iDataset)->Name()) << endl;
    //get the input file name
    TString inFileName = hstDir+"/"+*(listOfDatasets.at(iDataset)->File());
    //get a nice GPTree for this sample 
    MitGPTree thistree;
    thistree.LoadTree(inFileName,2);
    thistree.InitTree(0);

    // Loop on the tree and make all the computations
    float theMassCut = 0;
    int   theLepId = 13;

    Int_t nEntries = thistree.tree_->GetEntries();
    for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
      thistree.tree_-> GetEntry(iEntry);
      
      //Step1: select lep+lep- events
      if (thistree.nlep_ < 2) continue;
      if (fabs(thistree.lid1_) != fabs(thistree.lid2_)) continue;
      if (thistree.lid1_ * thistree.lid2_ > 0) continue;

      //Step2: select Zmumu events
      LorentzVector thisZBoson;
      thisZBoson = thistree.lep1_ + thistree.lep2_;

      //Step3: select the recoil photon
      if (thistree.nphotons_ < 1) continue;
      LorentzVector thisGamma;
      thisGamma = thistree.pho1_;

      //Step3: fill nutple
      ZTuple->Fill(thistree.event_,thistree.lumi_,thistree.run_,thisZBoson.M(),thistree.lid1_,
                   thisZBoson.Eta(),thisZBoson.Phi(),thisZBoson.Pt(), thisZBoson.E(),
                   thisGamma.Eta(),thisGamma.Phi(),thisGamma.Pt(), thisGamma.E());
    }
    
  }//end loop on samples
        
  //prepare output file
  TFile* outFile = new TFile(outFileName,"RECREATE");
  ZTuple->Write();
  outFile->Close();
  
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

