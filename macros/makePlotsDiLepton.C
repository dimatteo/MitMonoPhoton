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
#include <TLine.h>
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Style/interface/MitStyle.h"
#endif

using namespace std;
using namespace mithep;

TString getEnv(const char* name);
void makeStack(TString myVar, TString myCut, TString myName, TString myAxisNameX, TString myAxisNameY, 
               vector<const Sample*>& listOfSignals, vector<const Sample*>& listOfSamples, vector<const Sample*> listOfDatasets, 
               TString inFileName,
               bool isBlind, bool isLog, bool drawSignal, bool drawLegend,
               int nBins, float xLow, float xHigh,
               float* xlowVec = 0);

TH1* makeRatioBand(TH1* hist);
TH1* makeRatioPlot(TH1* num, TH1* den);

//==================================================================================================
void makePlotsDiLepton(bool isBlind=0, bool drawSignal=0)
{
  // setup graphics stuff before starting
  MitStyle::Init();

  // read all environment variables
  TString here   = getEnv("PWD");
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_PLOTZ_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((mitHgg + TString("/config")).Data());
  vector<const Sample*> listOfSignals;
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) {
    if (samples->GetSample(iSample)->Name()->Contains("addmpho") ||
        samples->GetSample(iSample)->Name()->Contains("dmmpho") )
      listOfSignals.push_back(samples->GetSample(iSample));  
    else
      listOfSamples.push_back(samples->GetSample(iSample));  
  }
  vector<const Sample*> listOfDatasets;
  for (UInt_t iSample=0; iSample < *samples->NDataSamples(); iSample++) listOfDatasets.push_back(samples->GetDataSample(iSample));
  // define infolder
  TString inFileName = here + "/monoph-2013-Oct1_reduced_DiLepton.root";
  std::cout << "inFileName " << inFileName << std::endl;
  
  //Make the histos
  //==Default cut==
  //TString thisCut = "*(phoEt > 160 && metBosCor > 140 && nalljets < 2 && ((nlep > 0 && lep1Mass < 0.05 && lepAsPho > 0.5 && bosonMass > 10) || (nlep > 1 && lep1Mass > 0.05 && lep1Mass > 0.05)))";
  //==Ele only check ele veto==
  TString thisCut = "*(phoEt > 160 && metBosCor > 140 && nalljets < 2 && ((nlep > 0 && lep1Mass < 0.05 && bosonMass > 60 && bosonMass < 120) ))";
  //TString thisCut = "*(phoEt > 160 && nalljets < 2 && nlep > 1 && bosonPt > 140 && lep1Mass > 0.05 && lep1Mass > 0.05)";
  //TString thisCut = "*(phoEt > 160 && metBosCor > 140 && nalljets < 2 && ((nlep > 0 && lep1Mass < 0.05 && lepAsPho > 0.5 && bosonMass > 10) || (nlep > 1 && lep1Mass > 0.05 && lep1Mass > 0.05)))";
  //TString thisCut = "*(phoEt > 160 && nalljets < 2 && bosonPt > 140 && ((nlep > 0 && lep1Mass < 0.05 && lepAsPho > 0.5 && bosonMass > 10) || (nlep > 1 && lep1Mass > 0.05 && lep1Mass > 0.05)))";
  //TString thisCut = "*(phoEt > 160 && nalljets < 2 && bosonPt > 140 && nlep > 1 && abs(lep1Id) == 11 && abs(lep1Id) == 11)";
  //TString thisCut = "*(phoEt > 160 && nlep > 0 && abs(lep1Id) == 13 && abs(lep2Id) == 13)";

  //TString thisCut = "*(phoEt > 160 && nlep > 1)";

  //// One bin histo (only for counting)
  //makeStack("nphotons", thisCut, "nphotons", "nphotons", "", 
               //listOfSignals, listOfSamples, listOfDatasets, 
               //inFileName,
               //isBlind, false, drawSignal, true,
               //1, -1, 10,
               //0);
  // lepMetDeltaPhi: fixed binning
  makeStack("abs(TMath::ACos(TMath::Cos(lep1Phi-metCorPhi)))", thisCut, "lepMetDeltaPhi", "#Delta#phi(lep-MET)", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, false,
               30, 0, 3.15,
               0);
  // phoMetDeltaPhi: fixed binning
  makeStack("phoMetDeltaPhi", thisCut, "phoMetDeltaPhi", "#Delta#phi(#gamma-MET)", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, false,
               30, 0, 3.15,
               0);
  makeStack("jetMetDeltaPhi", thisCut, "jetMetDeltaPhi", "#Delta#phi(jet-MET)", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, false,
               30, 0, 3.15,
               0);
  makeStack("phoJetDeltaPhi", thisCut, "phoJetDeltaPhi", "#Delta#phi(#gamma-jet)", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, false,
               30, 0, 3.15,
               0);
  // MET: variable binning
  float xlowMET[6] = {140,150,165,195,260,500};        
  makeStack("metBosCor", thisCut, "met", "MET [GeV]", "GeV", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               5, -1, -1,
               xlowMET);
  // photon ET: variable binning
  float xlowPhET[5] = {160,175,195,260,500};        
  makeStack("phoEt", thisCut, "phoEt", "#gamma_{ET} [GeV]", "GeV", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, true, drawSignal, true,
               4, -1, -1,
               xlowPhET);
  // MET/photonET: fixed binning
  makeStack("metBosCor/phoEt", thisCut, "metOverPhET", "MET/#gamma_{ET}", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               30, 0, 2,
               0);
  // nVtx: fixed binning
  makeStack("nvtx", thisCut, "nvtx", "N_{vtx}", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               20, 0, 40,
               0);
  // jetEt: fixed binning
  makeStack("jet1Pt", thisCut, "jetPt", "Lead Jet_{PT}", "GeV", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               40, 0, 200,
               0);
  // nJets: fixed binning
  makeStack("nalljets", thisCut, "njets", "N_{jets}", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               5, 0, 5,
               0);
  // metSig: fixed binning
  makeStack("metSig", thisCut, "metSig", "MET Significance", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               50, 0, 100,
               0);
  // bosonMass: fixed binning
  makeStack("bosonMass", thisCut, "bosonMass", "m(Z)", "GeV", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               80, 0, 400,
               0);
  // bosonPt: fixed binning
  makeStack("bosonPt", thisCut, "bosonPt", "p_{T}(Z)", "GeV", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               40, 0, 400,
               0);
  // bosonEta: fixed binning
  makeStack("bosonEta", thisCut, "bosonEta", "#eta(Z)", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               40, -4, 4,
               0);
  // bosonMass: fixed binning
  makeStack("bosonPhoMass", thisCut, "bosonPhoMass", "m(Z#gamma)", "GeV", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               40, 0, 1000,
               0);
  // ele->pho fake rate: fixed binning
  makeStack("lepAsPho", thisCut, "lepAsPho", "ele is pho", "", 
               listOfSignals, listOfSamples, listOfDatasets, 
               inFileName,
               isBlind, false, drawSignal, true,
               3, -1, 2,
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
               vector<const Sample*>& listOfSignals, vector<const Sample*>& listOfSamples, vector<const Sample*> listOfDatasets, 
               TString inFileName,
               bool isBlind, bool isLog, bool drawSignal, bool drawLegend,
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
  TH1F* hsignal;
  //prepare data and signal histos
  if (xlowVec != 0) hdata   = new TH1F("hdata","",nBins,xlowVec);
  else hdata = new TH1F("hdata","",nBins,xLow,xHigh);
  if (xlowVec != 0) hsignal = new TH1F("hsignal","",nBins,xlowVec);
  else hsignal = new TH1F("hsignal","",nBins,xLow,xHigh);

  TTree*  treedata[20];
  for (UInt_t iDatas=0; iDatas < listOfDatasets.size(); iDatas++) {
    //get the tree
    treedata[iDatas] = (TTree*) infile -> Get(listOfDatasets.at(iDatas)->Name()->Data());

    //fill the histogram
    if ( iDatas == 0 ) treedata[iDatas] -> Draw(myVar + " >> hdata","evt_weight*kf_weight*pu_weight" + myCut);
    else treedata[iDatas] -> Draw(myVar + " >>+ hdata","evt_weight*kf_weight*pu_weight" + myCut);
    
    if ( isBlind && iDatas == 0 ) leg -> AddEntry(hdata, "DATA (19.8 fb^{-1})", "pl");    
    
  }//end loop on datasets
  if (xlowVec != 0) {
    for (int iBin = 1; iBin <= nBins; iBin++) hdata->SetBinError  (iBin,hdata->GetBinError(iBin)/hdata->GetBinWidth(iBin));
    for (int iBin = 1; iBin <= nBins; iBin++) hdata->SetBinContent(iBin,hdata->GetBinContent(iBin)/hdata->GetBinWidth(iBin));
  }

  TTree*  treesignal[20];
  for (UInt_t iSignal=0; iSignal < listOfSignals.size(); iSignal++) {
    //get the tree
    treesignal[iSignal] = (TTree*) infile -> Get(listOfSignals.at(iSignal)->Name()->Data());

    //fill the histogram
    TString thisScale = Form("%f *", *(listOfSignals.at(iSignal)->Scale()));
    if ( iSignal == 0 ) treesignal[iSignal] -> Draw(myVar + " >> hsignal",thisScale + "evt_weight*kf_weight*pu_weight" + myCut);
    else treesignal[iSignal] -> Draw(myVar + " >>+ hsignal",thisScale + "evt_weight*kf_weight*pu_weight" + myCut);
    
    if ( drawSignal && iSignal == 0 ) leg -> AddEntry(hsignal, "Signal", "l");    
    
  }//end loop on signals
  if (xlowVec != 0) {
    for (int iBin = 1; iBin <= nBins; iBin++) hsignal->SetBinError  (iBin,hsignal->GetBinError(iBin)/hsignal->GetBinWidth(iBin));
    for (int iBin = 1; iBin <= nBins; iBin++) hsignal->SetBinContent(iBin,hsignal->GetBinContent(iBin)/hsignal->GetBinWidth(iBin));
  }
  hsignal -> SetLineColor(49);
  hsignal -> SetLineWidth(4.0);
       
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
    }

    //fill the histogram
    TString thisScale = Form("%f *", *(listOfSamples.at(iSample)->Scale()));
    if (isFirstOfSerie) tree[iSample] -> Draw(myVar + " >> " + TString(hist[theHistCounter] -> GetName()),thisScale + "evt_weight*kf_weight*pu_weight" + myCut);
    else tree[iSample] -> Draw(myVar + " >>+ " + TString(hist[theHistCounter] -> GetName()),thisScale + "evt_weight*kf_weight*pu_weight" + myCut);
    
    //add the histogram to the stack if the last of the series:
    //either last sample or ~ sample followed by non ~ sample
    if (isLastOfSerie) {
       if (xlowVec != 0) {
         for (int iBin = 1; iBin <= nBins; iBin++) hist[theHistCounter]->SetBinError  (iBin,hist[theHistCounter]->GetBinError(iBin)/hist[theHistCounter]->GetBinWidth(iBin));
         for (int iBin = 1; iBin <= nBins; iBin++) hist[theHistCounter]->SetBinContent(iBin,hist[theHistCounter]->GetBinContent(iBin)/hist[theHistCounter]->GetBinWidth(iBin));
       }
       hs -> Add(hist[theHistCounter]);
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
  float theMaxMC = theMCSum->GetBinContent(theMCSum->GetMaximumBin()) + theMCSum->GetBinError(theMCSum->GetMaximumBin());
  if (theMaxMC > theMax) theMax = theMaxMC;
  
  //prepare the ratio band and plot
  TH1* theMCRatioBand = makeRatioBand(theMCSum);
  TH1* theRatioPlot = makeRatioPlot(hdata,theMCSum);
    
  TCanvas* can = new TCanvas();
  can -> SetLogy(isLog);
  
  TPad *pad1 = new TPad("pad1","top pad",0,0.30,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.13);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","bottom pad",0,0.0,1,0.30);
  pad2->SetTopMargin(0.02);
  pad2->SetLeftMargin(0.13);
  pad2->SetBottomMargin(0.4);
  pad2->SetGridy();
  pad2->Draw();
  
  pad1->cd();
  hs->Draw("hist");
  hdata->Draw("same,pe");
  if (drawSignal) hsignal->Draw("same,hist");
  if (drawLegend) leg->Draw("same");
  //hs->GetXaxis()->SetTitle(myAxisNameX);
  hs->GetYaxis()->SetTitle(myAxisNameY);
  hs->GetXaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetXaxis()->SetLabelOffset(0.025);
  hs->GetYaxis()->SetLabelOffset(0.035);
  //hs->GetXaxis()->SetTitleOffset(1.1);
  hs->GetYaxis()->SetTitleOffset(1.1);
  hs->SetMaximum(theMax);
  if (isLog) hs->SetMinimum(0.01);
  
  pad2->cd();
  theMCRatioBand->GetXaxis()->SetTitle(myAxisNameX);
  theMCRatioBand->GetXaxis()->SetTitleSize(0.16);
  theMCRatioBand->GetXaxis()->SetTitleOffset(1.1);
  theMCRatioBand->GetXaxis()->SetLabelSize(0.12);
  theMCRatioBand->GetXaxis()->SetLabelOffset(0.07);
  theMCRatioBand->GetYaxis()->SetTitle("Data/MC");
  theMCRatioBand->GetYaxis()->SetTitleSize(0.10);
  theMCRatioBand->GetYaxis()->SetTitleOffset(0.6);
  theMCRatioBand->GetYaxis()->SetLabelSize(0.06);
  theMCRatioBand->GetYaxis()->SetLabelOffset(0.03);
  theMCRatioBand->SetFillStyle(3001);
  theMCRatioBand->SetFillColor(kBlue);
  theMCRatioBand->SetLineWidth(1);
  theMCRatioBand->SetLineColor(kBlack);
  theMCRatioBand->SetMarkerSize(0.1);
  theMCRatioBand->SetMaximum(4.);
  theMCRatioBand->SetMinimum(0.);
  theMCRatioBand->Draw("E2");
  TLine *line = new TLine(xLow,1,xHigh,1);
  line->SetLineColor(kBlack);
  line->Draw("same");
  theRatioPlot->Draw("same,pe");
  
  can->cd();
  can->Modified();
  can -> SaveAs(myName + ".pdf","pdf");
  
  //cleanup the memory allocation
  delete theMCSum;
  delete hs;
  delete leg;
  delete hdata;
  delete pad1;
  delete pad2;
  delete can;
  delete theMCRatioBand;
  delete theRatioPlot;
  infile -> Close();
  delete infile;
  
  return;
}

TH1* makeRatioBand(TH1* hist)
{
  TH1* theResult = (TH1*) hist -> Clone("ratioBand");
  for (int iBin = 1; iBin <= hist -> GetNbinsX(); iBin++) {
    float theError = 1.;
    if (hist->GetBinContent(iBin) > 0) theError = sqrt(pow(hist->GetBinError(iBin)/hist->GetBinContent(iBin),2)+0.07*0.07);
    theResult -> SetBinError(iBin, theError);
    theResult -> SetBinContent(iBin, 1.);    
  }
   
  theResult->SetBit(TH1::kNoTitle); 
  return theResult;  
}

TH1* makeRatioPlot(TH1* num, TH1* den)
{
  TH1* theResult = (TH1*) num -> Clone("ratioPlot");
  for (int iBin = 1; iBin <= num -> GetNbinsX(); iBin++) {
    float theVal = 0.;
    float theError = 0.;
    if (den->GetBinContent(iBin) > 0) {
      theVal = num->GetBinContent(iBin)/den->GetBinContent(iBin);
      theError = num->GetBinError(iBin)/den->GetBinContent(iBin);
    }
    theResult -> SetBinError(iBin, theError);
    theResult -> SetBinContent(iBin, theVal);    
  }
   
  theResult->SetBit(TH1::kNoTitle); 
  return theResult;      
}


