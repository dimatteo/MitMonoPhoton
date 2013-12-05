//--------------------------------------------------------------------------------------------------
// Print out tables of numbers and plots for synchronization purposes
//
// Authors: L. Di Matteo                                                                  (Sep 2013)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TSystem.h"
#include "TSystem.h"
#include <TTree.h>
#include <THStack.h>
#include <TLegend.h>
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitMonoPhoton/Core/MitGPTree.h"
#endif

using namespace std;
using namespace mithep;

TString getEnv(const char* name);

//==================================================================================================
void makeSyncTables()
{
  // read all environment variables
  TString here   = getEnv("PWD");
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_PLOT_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // Input file
  //TString sampleName = "s12-zjets-ptz100-v7a";
  TString sampleName = "s12-addmpho-md3_d2-v7a";
  TString inFileName = here + "/sync_" + sampleName + "_noskim_0000.root";
  std::cout << "inFileName " << inFileName << std::endl;
  inFileName = "monojet-2013-Oct10_s12-zjets-ptz100-v7a_noskim_0000.root";
  TFile* inFileRoot = new TFile(inFileName,"READ");
    
  //Output txt file
  ofstream txtoutfile;
  txtoutfile.open (sampleName+".txt");
  txtoutfile << "event#   photonEt   photonEta   MET   leadJetPt   leadJetEta\n";
  txtoutfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    
  //get the processed number of events histo
  TString allEvts = "hDAllEvents";  
  TH1D* hDAllEvents = dynamic_cast<TH1D*>(inFileRoot->FindObjectAny(allEvts));

  //get a nice GPTree for this sample 
  MitGPTree thistree;
  thistree.LoadTree(inFileName,0);
  thistree.InitTree(0);
  
  //prepare the counters
  float counts[7]; 
  counts[0] = hDAllEvents -> GetEntries();
  counts[1] = 0;
  counts[2] = 0;
  counts[3] = 0;
  counts[4] = 0;
  counts[5] = 0;
  counts[6] = 0;
  counts[7] = 0;
  counts[8] = 0;

  //prepare the histograms
  TH1F* photonEt = new TH1F("photonEt","photonEt",100,0,1000);
  TH1F* photonEta = new TH1F("photonEta","photonEta",50,-3.0,3.0);
  TH1F* MET = new TH1F("MET","MET",100,0,1000);
  
  // Loop on the tree and make all the computations
  Int_t nEntries = thistree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
    thistree.tree_-> GetEntry(iEntry);
    //Only for the first 10 events printout relevant quantites
    if (iEntry < 20) {
      txtoutfile << showpos << fixed << setprecision(3) 
		 << setw(9) << thistree.event_       << " "   
                 << setw(9) << thistree.pho1_.Et()   << " "
		 << setw(9) << thistree.pho1_.Eta()  << " "
 		 << setw(9) << thistree.metCor_         << " "
		 << setw(9) << thistree.jet1_.Pt()   << " "
		 << setw(9) << thistree.jet1_.Eta()  << "\n";
    }
    //good photon
    if (thistree.nphotons_ < 1) continue;
    counts[1] += 1.;
    //fill the distributions in case there is at least one photon
    photonEt -> Fill(thistree.pho1_.Et());
    photonEta -> Fill(thistree.pho1_.Eta());
    MET -> Fill(thistree.met_);
    //barrel photon
    if (fabs(thistree.pho1_.Eta()) > 1.479) continue;
    counts[2] += 1.;
    ////isolated photon: default category 1
    //float cutIso1 = 6.0; 
    //float cutIso2 = 10.0;
    //float cutIso3 = 3.8; 
    ////change category based on R9
    //if ( thistree.phoR9_a1_ < 0.94 ) {
      //cutIso1 = 4.7;
      //cutIso2 = 6.5;
      //cutIso3 = 2.5;
    //}
    //isolated photon: default category 1
    float cutIso1 = 1.5; 
    float cutIso2 = 1.+0.04*thistree.pho1_.Pt();
    float cutIso3 = 0.7+0.005*thistree.pho1_.Pt(); 
    if (thistree.phoCombIso1_a1_ > cutIso1) continue;
    if (thistree.phoCombIso2_a1_ > cutIso2) continue;
    if (thistree.phoCombIso3_a1_ > cutIso3) continue;
    counts[3] += 1.;
    //hard photon
    if (thistree.pho1_.Et() < 160) continue;
    counts[4] += 1.;
    //hard met
    if (thistree.metCor_ < 140) continue;
    counts[5] += 1.;
    //other photon cuts
    if (thistree.phoMipTotEnergy_a1_ > 6.3) continue;
    if (fabs(thistree.phoSeedTime_a1_) > 3.) continue;
    if (fabs(thistree.phoLeadTimeSpan_a1_) > 5.) continue;
    if (thistree.phoCoviEtaiEta_a1_ < 0.001) continue;
    if (thistree.phoCoviPhiiPhi_a1_ < 0.001) continue;
    if (thistree.phoR9_a1_ >= 1.) continue;
    if (thistree.phoPassEleVeto_a1_ == 0) continue;
    counts[6] += 1.;
    //jet veto
    if (thistree.jet1_.Pt() >= 100) continue;
    counts[7] += 1.;
    //cosmic veto
    if (thistree.ncosmics_ > 0) continue;
    counts[8] += 1.;
  }
  txtoutfile.close();
  
  cout << fixed << setprecision(3) << "Total               | " << (int)counts[0] << " | " << counts[0]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "nPhoton > 0         | " << (int)counts[1] << " | " << counts[1]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "photon[0] in Barrel | " << (int)counts[2] << " | " << counts[2]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "photon[0] iso       | " << (int)counts[3] << " | " << counts[3]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "photon[0] Et > 160  | " << (int)counts[4] << " | " << counts[4]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "met > 140           | " << (int)counts[5] << " | " << counts[5]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "photon[0] is good   | " << (int)counts[6] << " | " << counts[6]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "leadJet < 100       | " << (int)counts[7] << " | " << counts[7]/counts[0] * 100. << "%" << endl; 
  cout << fixed << setprecision(3) << "no cosmics          | " << (int)counts[8] << " | " << counts[8]/counts[0] * 100. << "%" << endl; 

  // Output root file
  TFile* outfile = new TFile(sampleName+".root","RECREATE");
  outfile -> cd();
  photonEt -> Write();
  photonEta -> Write();
  MET -> Write();
  outfile -> Close();

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
