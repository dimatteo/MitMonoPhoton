//--------------------------------------------------------------------------------------------------
// Make the final tables from the reduced trees
// Macro based on the plotting utility (1-bin plots)
// Authors: L. Di Matteo                                                                  (Aug 2013)
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
#endif

using namespace std;
using namespace mithep;

TString getEnv(const char* name);
void makeTable(TString myVar, TString myCut, TString myName, TString myAxisNameX, TString myAxisNameY, 
               vector<const Sample*>& listOfSamples, vector<const Sample*> listOfDatasets, TString inFileName,
               bool isBlind, bool isLog, 
               int nBins, float xLow, float xHigh,
               float* xlowVec = 0);

//==================================================================================================
void makeTables()
{
  // setup graphics stuff before starting
  MitStyle::Init();

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
  // define sample
  TString inFileName = "../monoph-2013-Dec13_reducedTest_Signal.root";
  std::cout << "inFileName " << inFileName << std::endl;
  
  //Make the tables
  // -- EG cuts --
  TString thisCut  = "*(phoEt > 145 && metCor > 140 && nlep == 0 && nalljets < 2 && phoMetDeltaPhi > 0.5)";
  thisCut += "*(phoWorstIso < 1.5 && metMin > 125)";
  
  // -- CiC cuts --
  //TString thisCut = "*(phoEt > 160 && metCor > 140 && nlep == 0 && phoMetDeltaPhi > 0.5 && (nalljets == 0 || (nalljets == 1 && jetMetDeltaPhi > 0.5 && jetMetDeltaPhi < 2.65 && phoJetDeltaPhi < 2.9)))";

  // -- MET filters with beam halo tight filters --
  //thisCut += "*(metFilterWord == 1023 || metFilterWord == 511)";
  // -- MET filters w.o. beam halo filters --
  thisCut += "*(metFilterWord == 1023 || metFilterWord == 511 || metFilterWord == 255)";
  cout << "table in blinded control region" << endl;
  makeTable("nphotons", thisCut, "nphotons", "nphotons", "", 
               listOfSamples, listOfDatasets, inFileName,
               1, false, 
               1, -1, 10,
               0);
  // full region, still blinded
  cout << "\n" << endl;
  cout << "table in full region" << endl;
  makeTable("nphotons", thisCut, "nphotons", "nphotons", "", 
               listOfSamples, listOfDatasets, inFileName,
               0, false, 
               1, -1, 10,
               0);
        
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

void makeTable(TString myVar, TString myCut, TString myName, TString myAxisNameX, TString myAxisNameY, 
               vector<const Sample*>& listOfSamples, vector<const Sample*> listOfDatasets, TString inFileName,
               bool isBlind, bool isLog, 
               int nBins, float xLow, float xHigh,
               float* xlowVec)
{
  // prepare the input file
  TFile* infile = new TFile(inFileName, "READ"); 
  infile -> cd();
  
  // prepare the stack
  THStack *hs = new THStack("hs","");
  // prepare the histos pointers
  TH1F*   hist[20];
  // prepare the tree pointers
  TTree*  tree[20];
  // prepare the legend
  TLegend* leg = new TLegend(.7485,.7225,.9597,.9604);
  leg->SetFillColor(0);
  // prepare the colors
  Int_t col[20] = {46,2,12,5,3,4,9,7,47,49,49,50,51,52,53,54,55,56,57,58};
  // prepare the cut
  if (isBlind) myCut += "*(phoMetDeltaPhi < 2.9)";     
  else  myCut += "*(phoMetDeltaPhi > 2)";
  // prepare the Y axis lable
  if (xlowVec != 0) myAxisNameY = "Events/" + myAxisNameY;
  else {
    float binWidth = (xHigh-xLow)/nBins;
    TString tempString;
    tempString.Form("%.2f ",binWidth); 
    myAxisNameY = "Events/" + tempString + myAxisNameY;
  }
  // prepare the legend strings
  vector<TString> theLegends;
  
  // loop through the datasets and produce the plots
  TH1F* hdata;
  //variable bin histo
  if (xlowVec != 0) hdata = new TH1F("hdata","",nBins,xlowVec);
  //fixed bin histo
  else hdata = new TH1F("hdata","",nBins,xLow,xHigh);

  TTree*  treedata[20];
  for (UInt_t iDatas=0; iDatas < listOfDatasets.size(); iDatas++) {
    //get the tree
    treedata[iDatas] = (TTree*) infile -> Get(listOfDatasets.at(iDatas)->Name()->Data());
  

    //fill the histogram
    if ( iDatas == 0 ) treedata[iDatas] -> Draw(myVar + " >> hdata","hlt_weight*evt_weight*kf_weight*pu_weight" + myCut,"goff");
    else treedata[iDatas] -> Draw(myVar + " >>+ hdata","hlt_weight*evt_weight*kf_weight*pu_weight" + myCut,"goff");
    
    if ( iDatas == 0 ) leg -> AddEntry(hdata, "DATA (19.8 fb^{-1})", "pl");    
    
  }//end loop on datasets
  if (xlowVec != 0) {
    for (int iBin = 1; iBin <= nBins; iBin++) hdata->SetBinError  (iBin,hdata->GetBinError(iBin)/hdata->GetBinWidth(iBin));
    for (int iBin = 1; iBin <= nBins; iBin++) hdata->SetBinContent(iBin,hdata->GetBinContent(iBin)/hdata->GetBinWidth(iBin));
  }
  
  int theHistCounter = 0;
  // loop through the samples and produce the plots
  for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {

    //determine if the histo is first of the series
    bool isFirstOfSerie = (*listOfSamples.at(iSample)->Legend()).CompareTo(" ");
    bool isLastOfSerie = false;
    if (iSample == listOfSamples.size() - 1) isLastOfSerie = true;
    if (iSample < listOfSamples.size() - 1 && (*listOfSamples.at(iSample+1)->Legend()).CompareTo(" ") != 0) isLastOfSerie = true;
    
    //get the tree
    tree[iSample] = (TTree*) infile -> Get(listOfSamples.at(iSample)->Name()->Data());
    //if sample first of the list create a new histogram
    if (isFirstOfSerie) {
       TString thisHistName = "h_" + *(listOfSamples.at(iSample)->Name());
       //variable bin histo
       if (xlowVec != 0) hist[theHistCounter] = new TH1F(thisHistName,thisHistName,nBins,xlowVec);
       //fixed bin histo
       else hist[theHistCounter] = new TH1F(thisHistName,thisHistName,nBins,xLow,xHigh);
       hist[theHistCounter] -> Sumw2();
       hist[theHistCounter] -> SetFillColor(col[theHistCounter]);
       hist[theHistCounter] -> SetFillStyle(1001);
       theLegends.push_back(*listOfSamples.at(iSample)->Legend());
       cout << *listOfSamples.at(iSample)->Legend() << " ";
    }

    //fill the histogram
    TString thisScale = Form("%f *", *(listOfSamples.at(iSample)->Scale()));
    if (isFirstOfSerie) tree[iSample] -> Draw(myVar + " >> " + TString(hist[theHistCounter] -> GetName()),thisScale + "hlt_weight*evt_weight*kf_weight*pu_weight" + myCut,"goff");
    else tree[iSample] -> Draw(myVar + " >>+ " + TString(hist[theHistCounter] -> GetName()),thisScale + "hlt_weight*evt_weight*kf_weight*pu_weight" + myCut,"goff");
    
    //add the systematic uncertainty
    for (int iBin = 1; iBin <= nBins; iBin++) {
      float thisSyst = *(listOfSamples.at(iSample)->Syst());
      float thisError = sqrt(hist[theHistCounter]->GetBinError(iBin)*hist[theHistCounter]->GetBinError(iBin) + 
      thisSyst*thisSyst*hist[theHistCounter]->GetBinContent(iBin)*hist[theHistCounter]->GetBinContent(iBin));
      hist[theHistCounter]->SetBinError  (iBin,thisError);
    }

    //add the histogram to the stack if the last of the series:
    //either last sample or ~ sample followed by non ~ sample
    if (isLastOfSerie) {
       if (xlowVec != 0) {
         for (int iBin = 1; iBin <= nBins; iBin++) hist[theHistCounter]->SetBinError  (iBin,hist[theHistCounter]->GetBinError(iBin)/hist[theHistCounter]->GetBinWidth(iBin));
         for (int iBin = 1; iBin <= nBins; iBin++) hist[theHistCounter]->SetBinContent(iBin,hist[theHistCounter]->GetBinContent(iBin)/hist[theHistCounter]->GetBinWidth(iBin));
       }
       hs -> Add(hist[theHistCounter]);
       cout << hist[theHistCounter] -> GetBinContent(1) << " +/- " <<  hist[theHistCounter] -> GetBinError(1) << endl;
       theHistCounter++;
    }
    
  }//end loop on samples

  //Fix the legend
  for (int iHisto = theHistCounter-1; iHisto >= 0; iHisto--) {
    leg -> AddEntry(hist[iHisto], theLegends[iHisto], "f");   
  }
  
  //get the maximum to properly set the frame
  float theMax = hdata -> GetBinContent(hdata -> GetMaximumBin()) + hdata -> GetBinError(hdata -> GetMaximumBin());
  TH1* theMCSum = (TH1*) hs->GetStack()->Last();
  cout << "total MC " << theMCSum -> GetBinContent(1) << " +/- " << theMCSum -> GetBinError(1) << endl;
  float theMaxMC = theMCSum->GetBinContent(theMCSum->GetMaximumBin()) + theMCSum->GetBinError(theMCSum->GetMaximumBin());
  if (theMaxMC > theMax) theMax = theMaxMC;
  
    
  if (isBlind) cout << "total DATA " << hdata -> Integral() << " +/- " <<  hdata -> GetBinError(1) << endl;
  
  //cleanup the memory allocation
  delete theMCSum;
  delete hs;
  delete leg;
  delete hdata;
  infile -> Close();
  delete infile;
  
  return;
}
