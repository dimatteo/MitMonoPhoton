//--------------------------------------------------------------------------------------------------
// Make the final plots from the reduced trees
//
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
void makeStack(TString myVar, TString myCut, TString myName, TString myAxisNameX, TString myAxisNameY, 
               vector<const Sample*>& listOfSamples, vector<const Sample*> listOfDatasets, TString inFileName,
               bool isBlind, bool isLog, 
               int nBins, float xLow, float xHigh,
               float* xlowVec = 0);

//==================================================================================================
void makePlots(bool isBlind=1)
{
  // setup graphics stuff before starting
  MitStyle::Init();

  // read all environment variables
  TString here   = getEnv("PWD");
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_ANA_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  cout << hstDir.Data() << endl;
  samples->ReadFile((mitHgg + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  
  vector<const Sample*> listOfDatasets;
  for (UInt_t iSample=0; iSample < *samples->NDataSamples(); iSample++) listOfDatasets.push_back(samples->GetDataSample(iSample));
  // define infolder
  TString inFileName = here + "/monoph-2013-July9_reduced.root";
  std::cout << "inFileName " << inFileName << std::endl;
  
  //Make the histos
  
  // MET: variable binning
  float xlowMET[7] = {140,160,190,250,400,700,1000};        
  makeStack("met", "*(1.)", "met", "MET [GeV]", "GeV", 
               listOfSamples, listOfDatasets, inFileName,
               isBlind, true, 
               6, -1, -1,
               xlowMET);
  // photon ET: variable binning
  float xlowPhET[6] = {160,190,250,400,700,1000};        
  makeStack("phoEt", "*(1.)", "phoEt", "#gamma_{ET} [GeV]", "GeV", 
               listOfSamples, listOfDatasets, inFileName,
               isBlind, true, 
               5, -1, -1,
               xlowPhET);
  // MET/photonET: fixed binning
  makeStack("met/phoEt", "*(1.)", "metOverPhET", "MET/#gamma_{ET}", "", 
               listOfSamples, listOfDatasets, inFileName,
               isBlind, false, 
               30, 0, 2,
               0);
  // nVtx: fixed binning
  makeStack("nvtx", "*(1.)", "nvtx", "N_{vtx}", "", 
               listOfSamples, listOfDatasets, inFileName,
               isBlind, false, 
               20, 0, 40,
               0);
  // jetEt: fixed binning
  makeStack("jetEt", "*(1.)", "jetEt", "Lead Jet_{ET}", "GeV", 
               listOfSamples, listOfDatasets, inFileName,
               isBlind, false, 
               40, 0, 200,
               0);
  // jetEta: fixed binning
  makeStack("jetEta", "*(1.)", "jetEta", "Lead Jet_{#eta}", "", 
               listOfSamples, listOfDatasets, inFileName,
               isBlind, false, 
               40, -5, 5,
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

void makeStack(TString myVar, TString myCut, TString myName, TString myAxisNameX, TString myAxisNameY, 
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
  if (isBlind) myCut += "*( phoMetDeltaPhi < 2.89)";        
  // prepare the Y axis lable
  if (xlowVec != 0) myAxisNameY = "Events/" + myAxisNameY;
  else {
    float binWidth = (xHigh-xLow)/nBins;
    TString tempString;
    tempString.Form("%.2f ",binWidth); 
    myAxisNameY = "Events/" + tempString + myAxisNameY;
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
       hist[theHistCounter] -> SetFillColor(col[theHistCounter]);
       hist[theHistCounter] -> SetFillStyle(1001);
       leg -> AddEntry(hist[theHistCounter], *listOfSamples.at(iSample)->Legend(), "f");   
    }

    //fill the histogram
    if (isFirstOfSerie) tree[iSample] -> Draw(myVar + " >> " + TString(hist[theHistCounter] -> GetName()),"evt_weight*kf_weight*pu_weight" + myCut);
    else tree[iSample] -> Draw(myVar + " >>+ " + TString(hist[theHistCounter] -> GetName()),"evt_weight*kf_weight*pu_weight" + myCut);
    
    //add the histogram to the stack if the last of the series:
    //either last sample or ~ sample followed by non ~ sample
    if (isLastOfSerie) {
       if (xlowVec != 0) {
         for (int iBin = 1; iBin < nBins; iBin++) hist[theHistCounter]->SetBinError  (iBin,hist[theHistCounter]->GetBinError(iBin)/hist[theHistCounter]->GetBinWidth(iBin));
         for (int iBin = 1; iBin < nBins; iBin++) hist[theHistCounter]->SetBinContent(iBin,hist[theHistCounter]->GetBinContent(iBin)/hist[theHistCounter]->GetBinWidth(iBin));
       }
       hs -> Add(hist[theHistCounter]);
       theHistCounter++;
    }
    
  }//end loop on samples

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
    if ( iDatas == 0 ) treedata[iDatas] -> Draw(myVar + " >> hdata","evt_weight*kf_weight*pu_weight" + myCut);
    else treedata[iDatas] -> Draw(myVar + " >>+ hdata","evt_weight*kf_weight*pu_weight" + myCut);
    
    if ( iDatas == 0 ) leg -> AddEntry(hdata, "DATA", "pl");    
    
  }//end loop on datasets
  //cout << "data " << hdata -> Integral() << endl;
  if (xlowVec != 0) {
    for (int iBin = 1; iBin < nBins; iBin++) hdata->SetBinError  (iBin,hdata->GetBinError(iBin)/hdata->GetBinWidth(iBin));
    for (int iBin = 1; iBin < nBins; iBin++) hdata->SetBinContent(iBin,hdata->GetBinContent(iBin)/hdata->GetBinWidth(iBin));
  }
  
  //get the maximum to properly set the frame
  float theMax = hdata -> GetBinContent(hdata -> GetMaximumBin()) + hdata -> GetBinError(hdata -> GetMaximumBin());
  TH1* theMCSum = (TH1*) hs->GetStack()->Last();
  float theMaxMC = theMCSum->GetBinContent(theMCSum->GetMaximumBin()) + theMCSum->GetBinError(theMCSum->GetMaximumBin());
  if (theMaxMC > theMax) theMax = theMaxMC;
  
    
  TCanvas* can = new TCanvas();
  can -> SetLogy(isLog);
  hs->Draw("hist");
  hdata->Draw("same,pe");
  leg->Draw("same");
  hs->GetXaxis()->SetTitle(myAxisNameX);
  hs->GetYaxis()->SetTitle(myAxisNameY);
  hs->GetXaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetXaxis()->SetLabelOffset(0.025);
  hs->GetYaxis()->SetLabelOffset(0.025);
  hs->GetXaxis()->SetTitleOffset(1.1);
  hs->GetYaxis()->SetTitleOffset(1.3);
  hs->SetMaximum(theMax);
  if (isLog) hs->SetMinimum(0.01);
  can->Modified();
  can -> SaveAs(myName + ".pdf","pdf");
  
  //cleanup the memory allocation
  delete theMCSum;
  delete hs;
  delete leg;
  delete hdata;
  delete can;
  infile -> Close();
  delete infile;
  
  return;
}
