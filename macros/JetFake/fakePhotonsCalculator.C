#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TPad.h"
#include <sstream>
#include <iostream>

const char* preFitFile = "sysBin.root";

const char* concat(const char* a,const char* b,const char* c="",const char* d="") {
  std::stringstream s;
  s<<a<<b<<c<<d;
  const char* e = s.str().c_str();
  return e;
}

const char* cut(const char* word,Int_t start,Int_t finish,Bool_t end=kFALSE) {
  std::stringstream ss;
  ss<<word;
  std::string S = ss.str();
  if ( end ) 
    finish = strlen(word);
  std::string s = S.substr(start,finish);
  const char* c = s.c_str();
  return c;
}

void fractionFit( Int_t whichBin , Bool_t print=kFALSE, const char* fitPrint="V") {

  TFile *open = TFile::Open(preFitFile);
  const Int_t nBins = 6;
  const Double_t xBins[nBins+1] = {
    145,
    160,
    190,
    250,
    400,
    700,
    1000,
  };
  const char* binNums[nBins] = {
    "145_160",
    "160_190",
    "190_250",
    "250_400",
    "400_700",
    "700_1000",
  };
  const Int_t nWords = 3;
  const char* binWords[nWords] = {
    "data",
    "fake",
    "mc",
  };
  const char* binNames[nBins*nWords];
  TObjArray *histArray = new TObjArray();
  Double_t scaleSB[nBins]; Double_t errorSB[nBins];
  Double_t scaleMC[nBins]; Double_t errorMC[nBins];
  Int_t hist = 0;
  Int_t data = 0;
  Int_t fake = 1;
  Int_t mc   = 2;

  Double_t correction[nBins];
  Double_t initialFakeRatio[nBins];
  Double_t fakeRatio[nBins];
  Double_t initialFakePhotons[nBins];
  Double_t fakePhotons[nBins];
  Double_t totalInitialFakePhotons = 0;
  Double_t totalFakePhotons = 0;

  TH1D *numerators     = (TH1D*)open->Get("numerators");
  TH1D *denominators   = (TH1D*)open->Get("denominators");
  TH1D *looseSelection = (TH1D*)open->Get("looseSelection");
  TH1D *fakeRate       = new TH1D("fakeRate","Corrected Fake Rate",nBins,xBins);
  TH1D *corrNums       = new TH1D("corrNums","",nBins,xBins);

  TH1D *fitFake;
  TH1D *fitMC;
  TH1D *fitData;
  TH1D *normFake;
  TH1D *normMC;
  TH1D *normData;

  // TCanvas *canvas = new TCanvas("canvas",binNums[whichBin],1200,400); canvas->Divide(2,1);
  
  TFile *save = new TFile("test.root","RECREATE");

  for ( Int_t i=0; i<nBins; i++ ) {
    for ( Int_t j=0; j<nWords; j++ ) {
      binNames[hist] = concat(binWords[j],binNums[i]);    
      histArray->Add((TH1D*)open->Get(binNames[hist]));
      hist++;
    }
    TObjArray *forFit = new TObjArray();
    forFit->Add((TH1D*)histArray->At(fake)); forFit->Add((TH1D*)histArray->At(mc));

    TFractionFitter *fit = new TFractionFitter((TH1D*)histArray->At(data),forFit,fitPrint);
    fit->Fit();
    fit->GetResult(0,scaleSB[i],errorSB[i]);
    fit->GetResult(1,scaleMC[i],errorMC[i]);

    // Set scale manually if necessary
    // scaleSB[5] = 5.73320e-02; scaleMC[5] = 9.42680e-01;

    TH1D *dataHist = (TH1D*)histArray->At(data);
    TH1D *fakeHist = (TH1D*)histArray->At(fake);
    TH1D *mcHist   = (TH1D*)histArray->At(mc);

    if ( i == whichBin ) {
      normFake = (TH1D*)fakeHist->Clone("normFake");
      normMC   = (TH1D*)mcHist->Clone("normMC");
      normData = (TH1D*)dataHist->Clone("normData");
    }

    TH1D *scaledFake = (TH1D*)fakeHist->Clone("scaledFake");
    TH1D *scaledMC   = (TH1D*)mcHist->Clone("scaledMC");
 
    scaledFake->Scale(TMath::Power(scaledFake->Integral(1,200),-1));
    scaledMC  ->Scale(TMath::Power(scaledMC->Integral(1,200),-1));
    scaledFake->Scale(dataHist->Integral(1,200));
    scaledMC  ->Scale(dataHist->Integral(1,200));
    scaledFake->Scale(scaleSB[i]);
    scaledMC  ->Scale(scaleMC[i]);

    correction[i] = scaledFake->Integral(1,65)/(scaledFake->Integral(1,65) + scaledMC->Integral(1,65));

    if ( i == whichBin ) {
      fitFake = (TH1D*)scaledFake->Clone("fitFake");
      fitData = (TH1D*)dataHist->Clone("fitData");      
      fitMC   = (TH1D*)scaledMC->Clone("fitMC");
    }

    data+=3; fake+=3; mc+=3;

    initialFakeRatio[i]      = (numerators->GetBinContent(i+1))/(denominators->GetBinContent(i+1));
    if ( denominators->GetBinContent(i+1) == 0 )
      initialFakeRatio[i]    = 0;
    initialFakePhotons[i]    = (looseSelection->GetBinContent(i+1))*initialFakeRatio[i];
    fakeRatio[i]             = correction[i]*initialFakeRatio[i];
    fakePhotons[i]           = initialFakePhotons[i]*correction[i];
    totalFakePhotons        += fakePhotons[i];
    totalInitialFakePhotons += initialFakePhotons[i];

    if ( i > 0 )
      fakeRate->Fill(xBins[i], fakeRatio[i]);
    corrNums->Fill(xBins[i],correction[i]*(Double_t)numerators->GetBinContent(i+1));

    if ( print == kTRUE) {
    std::cout<<"Bin "<<i<<": "<<"Initial Fake Ratio     = "<<initialFakeRatio[i]<<std::endl;
    std::cout<<"       Initial Fake Photons   = "          <<initialFakePhotons[i]<<std::endl;
    std::cout<<"       Corrected Fake Ratio   = "          <<fakeRatio[i]<<std::endl;
    std::cout<<"       Corrected Photon Count = "          <<fakePhotons[i]<<std::endl;
    std::cout<<                                              std::endl;
    }
  }
  
  TFile *saveFakeRate = new TFile("fakePhotonsOutput.root","RECREATE");
  Double_t valHist[nBins];
  Double_t errorHist[nBins];
  Double_t errorComb[nBins];

  for ( Int_t i=0; i<nBins; i++ ) {
    valHist[i]   = (Double_t)corrNums->GetBinContent(xBins[i]);
    errorHist[i] = (Double_t)corrNums->GetBinError(xBins[i]);

    errorComb[i] = TMath::Sqrt(valHist[i]*valHist[i]*TMath::Power(1/(scaleSB[i]+scaleMC[i])*TMath::Sqrt(errorSB[i]*errorSB[i] + (errorSB[i]*errorSB[i]+errorMC[i]*errorMC[i])*scaleSB[i]*scaleSB[i]/((scaleSB[i]+scaleMC[i])*(scaleSB[i]+scaleMC[i]))),2) + errorHist[i]*errorHist[i]*scaleSB[i]*scaleSB[i]/((scaleSB[i]+scaleMC[i])*(scaleSB[i]+scaleMC[i])));

    corrNums->SetBinError(xBins[i],errorComb[i]);
  }
  fakeRate    ->Write();
  corrNums    ->Write();
  denominators->Write();
  saveFakeRate->Close();

  if ( print == kTRUE ) {
    std::cout<<"Total, corrected fake photons: "<<totalFakePhotons<<std::endl;
  }

  // canvas->cd(1);
  fitFake->SetLineColor(kRed); fitFake->SetFillColor(kRed); fitFake->SetFillStyle(3354);
  fitMC  ->SetLineColor(kBlue); fitMC->SetFillColor(kBlue); fitMC->SetFillStyle(3345);
  fitData->SetLineColor(kBlack); fitData->SetMarkerStyle(20);
  // fitFake->SetTitle(concat("Fit P^{#gamma}_{T} Bin: ",binNums[whichBin]," GeV"));
  fitFake->SetTitle("");
  fitFake->GetXaxis()->SetRange(0,125);
  fitFake->GetXaxis()->SetTitle("#sigma_{i #eta i #eta}"); fitFake->GetXaxis()->SetTitleSize(0.05);
  fitFake->GetYaxis()->SetTitle("Number of Events"); fitFake->GetXaxis()->SetTitleSize(0.05);
  fitFake->Add(fitMC);
  fitFake->Draw("");
  fitMC  ->Draw("same");
  fitData->Draw("same P E1 ][");
  gPad   ->SetLogy();
  gStyle ->SetOptStat("");
  TLegend *legend = new TLegend(0.7,0.7,0.89,0.89);
  legend->SetHeader(concat(cut(binNums[whichBin],0,3)," #leq P^{#gamma}_{T} < ",cut(binNums[whichBin],4,4,kTRUE)," [GeV]"));
  legend->AddEntry("fitData","Numerator (Data)","p");
  legend->AddEntry("fitFake","QCD Sideband (Data)","f");
  legend->AddEntry("fitMC","#gamma + Jet (MC)","f");
  legend->SetFillColor(0); legend->SetLineColor(0);
  legend->Draw();
  ((TLegendEntry*)legend->GetListOfPrimitives()->First())->SetTextSize(0.03);



  // canvas->cd(2);
  // normFake->SetLineColor(kRed); normFake->SetFillColor(kRed); normFake->SetFillStyle(3354);
  // normMC->SetLineColor(kBlue); normFake->SetFillColor(kBlue); normMC->SetFillStyle(3345);
  // normData->SetLineColor(kBlack); normData->SetMarkerStyle(20);  
  // normMC->SetTitle(concat("Normalized ",binNums[whichBin]));
  // normMC->GetXaxis()->SetTitle("sig_ietaieta");
  // normMC->GetYaxis()->SetTitle("Number of Events");
  // normMC->DrawNormalized();
  // normFake->DrawNormalized("same");
  // normData->DrawNormalized("same P E1 ][");
  
  open->Close();
}

void fakeRate() {

  const Int_t nBins  = 6;
  const Double_t xBins[nBins+1] = {
    130,
    140,
    150,
    160,
    200,
    300,
    1000,
  };
  TFile       *openRate  = TFile::Open("sysBin.root");
  TH1D        *ratioHist = (TH1D*)openRate->Get("fakeRate");
  TH1D        *corrNums  = (TH1D*)openRate->Get("corrNums");
  TH1D        *dens      = (TH1D*)openRate->Get("denominators");
  TF1         *ratioFunc = new TF1("ratioFunc","pol1",0,1000);
  TEfficiency *eff = new TEfficiency(*corrNums,*dens);

  // ratioHist->Fit("ratioFunc");
  // eff->Fit(ratioFunc);
  ratioFunc->SetParameter(0,7.99101e-03);// - 2.72194e-04);
  ratioFunc->SetParameter(1,-9.71701e-06);// - 1.12266e-06);
  
  TH1 *effHist = (TH1*)eff->GetPassedHistogram();
  effHist->Divide((TH1*)eff->GetTotalHistogram());
  effHist->SetTitle("");
  effHist->GetXaxis()->SetTitle("P^{#gamma}_{T}");
  effHist->GetYaxis()->SetTitle("Fake Rate");
  // effHist->Fit("ratioFunc");
  effHist->SetMinimum(0.001);
  effHist->SetMaximum(0.01);
  effHist->SetLineColor(kBlack);

  effHist->Draw("E1");
  ratioFunc->Draw("same");

  gStyle->SetOptStat("");
  // ratioHist->Draw("same");

  TFile *openData = TFile::Open(preFitFile);
  TH1D  *phoPt    = (TH1D*)openData->Get("loosePhotons"); 
  Int_t  ptBins   = phoPt->GetNbinsX();
  Int_t nLoosePho = 0;
  Double_t fakes  = 0;
  for ( Int_t i=1; i<ptBins+1; i++ ) {
    if ( phoPt->GetBinContent(i) == 0 ) continue;
    nLoosePho += phoPt->GetBinContent(i);
    fakes += phoPt->GetBinContent(i)*ratioFunc->Eval(phoPt->GetBinCenter(i));
  }
  std::cout<<std::endl;
  std::cout<<"Estimated fake photons = "<<fakes<<std::endl;
  std::cout<<"Out of "<<nLoosePho<<" loose photons"<<std::endl;
  std::cout<<std::endl;
}

void numDen() {
  TFile *open = TFile::Open(preFitFile);

  TH1D *numerators     = (TH1D*)open->Get("numerators");
  TH1D *denominators   = (TH1D*)open->Get("denominators");

  numerators  ->SetLineColor(kRed);
  denominators->SetLineColor(kBlue);
  denominators->SetTitle("");
  denominators->GetYaxis()->SetTitle("Number of Events"); denominators->GetYaxis()->SetTitleSize(0.04);
  denominators->GetXaxis()->SetTitle("P^{#gamma}_{T}");   denominators->GetXaxis()->SetTitleSize(0.04);
  denominators->Draw("][");
  numerators  ->Draw("same ][");
  // gPad        ->SetLogy();
  gStyle      ->SetOptStat("");
  TLegend *ifLegend = new TLegend(0.65,0.7,0.89,0.89);
  ifLegend->SetHeader("Uncorrected Fake Rate");
  ifLegend->AddEntry("denominators","Denominator","l");
  ifLegend->AddEntry("numerators","Numerator","l");
  ifLegend->SetTextSize(0.02);
  ifLegend->SetFillColor(0); ifLegend->SetLineColor(0);
  ifLegend->Draw();
  ((TLegendEntry*)ifLegend->GetListOfPrimitives()->First())->SetTextSize(0.03);
}

void fakePhotonsCalculator(){
  std::cout<<"Try fakeRate() or fractionFit()"<<std::endl;
}
