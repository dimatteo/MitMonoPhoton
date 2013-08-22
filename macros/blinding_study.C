// blinding_study.C
//
// This macro studies the blinding cut on deltaPhi(photon-met) for one signal
// sample, defined by signal_num. It uses data files from my scratch area.
// With the signal sample determined, the macro produces a pseudodataset
// using scaled contributions from the available background sources and the
// given signal sample. It then makes a histogram of the signal significance
// as a function of the cut on deltaPhi(photon-met), and prints the deltaPhi
// cut needed to reduce the signal significance by a factor of 4.
//
// Christian Ferko

#include "TMath.h"
#include "TTree.h"
#include "MitGPTree.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
using namespace std;

void blinding_study() 
{
	gSystem->CompileMacro("MitGPTree.h");
	
// First we define MIT Style for the plots.

	TStyle *MitStyle = gStyle;
	//gStyle = MitStyle;

	// Canvas
	MitStyle->SetCanvasColor     (0);
	MitStyle->SetCanvasBorderSize(10);
	MitStyle->SetCanvasBorderMode(0);
	MitStyle->SetCanvasDefH      (700);
	MitStyle->SetCanvasDefW      (700);
	MitStyle->SetCanvasDefX      (100);
	MitStyle->SetCanvasDefY      (100);

	// Pads
	MitStyle->SetPadColor       (0);
	MitStyle->SetPadBorderSize  (10);
	MitStyle->SetPadBorderMode  (0);
	MitStyle->SetPadBottomMargin(0.13);
	MitStyle->SetPadTopMargin   (0.04);
	MitStyle->SetPadLeftMargin  (0.18);
	MitStyle->SetPadRightMargin (0.04);
	MitStyle->SetPadGridX       (0);
	MitStyle->SetPadGridY       (0);
	MitStyle->SetPadTickX       (0);
	MitStyle->SetPadTickY       (0);

	// Frames
	MitStyle->SetFrameFillStyle ( 0);
	MitStyle->SetFrameFillColor ( 0);
	MitStyle->SetFrameLineColor ( 1);
	MitStyle->SetFrameLineStyle ( 0);
	MitStyle->SetFrameLineWidth ( 1);
	MitStyle->SetFrameBorderSize(10);
	MitStyle->SetFrameBorderMode( 0);

	// Histograms
	MitStyle->SetHistFillColor(2);
	MitStyle->SetHistFillStyle(0);
	MitStyle->SetHistLineColor(1);
	MitStyle->SetHistLineStyle(0);
	MitStyle->SetHistLineWidth(2);
	MitStyle->SetNdivisions(505);

	// Functions
	MitStyle->SetFuncColor(1);
	MitStyle->SetFuncStyle(0);
	MitStyle->SetFuncWidth(2);

	// Various
	MitStyle->SetMarkerStyle(20);
	MitStyle->SetMarkerColor(kBlack);
	MitStyle->SetMarkerSize (1.2);

	MitStyle->SetTitleSize  (0.055,"X");
	MitStyle->SetTitleOffset(1.200,"X");
	MitStyle->SetLabelOffset(0.005,"X");
	MitStyle->SetLabelSize  (0.050,"X");
	MitStyle->SetLabelFont  (42   ,"X");
	MitStyle->SetTickLength (-0.03,"X");

	MitStyle->SetStripDecimals(kFALSE);

	MitStyle->SetTitleSize  (0.055,"Y");
	MitStyle->SetTitleOffset(1.800,"Y");
	MitStyle->SetLabelOffset(0.010,"Y");
	MitStyle->SetLabelSize  (0.050,"Y");
	MitStyle->SetLabelFont  (42   ,"Y");
	MitStyle->SetTickLength (-0.03,"Y");

	MitStyle->SetTextSize   (0.055);
	MitStyle->SetTextFont   (42);

	MitStyle->SetStatFont   (42);
	MitStyle->SetTitleFont  (42);
	MitStyle->SetTitleFont  (42,"X");
	MitStyle->SetTitleFont  (42,"Y");

	MitStyle->SetOptStat    (0);
  
// Here the style section ends and the macro begins.
	
	string sig_samples[] = 
	{
		"s12-dmmpho-v_m1-v7a",
		"s12-dmmpho-av_m1-v7a",
		"s12-dmmpho-v_m10-v7a",
		"s12-dmmpho-av_m10-v7a",
		"s12-dmmpho-v_m100-v7a",
		"s12-dmmpho-av_m100-v7a",
		"s12-dmmpho-v_m200-v7a",
		"s12-dmmpho-av_m200-v7a",
		"s12-dmmpho-av_m300-v7a",	
		"s12-dmmpho-v_m500-v7a",
		"s12-dmmpho-av_m500-v7a",		
		"s12-dmmpho-v_m1000-v7a",
		"s12-dmmpho-av_m1000-v7a",
		"s12-addmpho-md1_d2-v7a",
		"s12-addmpho-md1_d3-v7a",
		"s12-addmpho-md1_d4-v7a",
		"s12-addmpho-md1_d5-v7a",		
		"s12-addmpho-md1_d6-v7a",
		"s12-addmpho-md2_d2-v7a",
		"s12-addmpho-md2_d3-v7a",
		"s12-addmpho-md2_d5-v7a",
		"s12-addmpho-md2_d6-v7a",
		"s12-addmpho-md3_d2-v7a",
		"s12-addmpho-md3_d3-v7a",
		"s12-addmpho-md3_d4-v7a",
		"s12-addmpho-md3_d5-v7a",
		"s12-addmpho-md3_d6-v7a"
	}
	
	Int_t signal_num=0; 
	
// This macro considers only one signal sample at a time, so the variable signal_num defines
// which signal sample to work on (it is the index of the sample in sig_samples).
	
	double sig_weights[] =
	{
		4.81E-07,
		4.79E-07,
		4.80E-07,
		4.81E-07,
		4.77E-07,
		4.26E-07,
		4.23E-07,
		3.18E-07,
		2.19E-07,
		2.01E-07,
		9.52E-08,
		3.00E-08,
		8.22E-09,
		6.48E-02,
		1.73E-01,
		6.93E-02,
		9.17E-02,
		9.63E-01,
		5.87E-03,
		4.91E-03,
		4.37E-03,
		4.26E-03,
		1.45E-03,
		7.93E-04,
		5.53E-04,
		4.26E-04,
		3.24E-04
	};
	
// The signal weights are given by (sigma_MC * lumi)/(N_processed). Background weights are defined similarly below.

	string bg_samples[] = {"s12-zgptg130-v7c", "s12-wjets-ptw100-v7a", "s12-wgptg130-v7a", "s12-qcdht100-250-v7a", "s12-qcdht250-500-v7a", "s12-qcdht500-1000-v7a", 
	"s12-qcdht1000-v7a", "s12-2pibo10_25-v7a", "s12-2pibo25_250-v7a", "s12-2pibo250-v7a", "s12-2pibx10_25-v7a", "s12-2pibx25_250-v7a", "s12-2pibx250-v7a", 
	"s12-pj50_80-v7a", "s12-pj80_120-v7a", "s12-pj120_170-v7a", "s12-pj170_300-v7a", "s12-pj300_470-v7a", "s12-pj470_800-v7a", "s12-pj800_1400-v7a", 
	"s12-pj1400_1800-v7a", "s12-pj1800-v7a", "s12-zgllgptg130-v7a", "s12-zgllg-v7a"};

	double bg_weights[]=
	{
		5.28E-03,
		1.14E-01,
		1.38E-02,
		5.63E+03,
		2.34E+02,
		9.71E+00,
		3.09E-01,
		9.21E+00,
		9.90E-01,
		4.21E-04,
		1.66E+01,
		6.43E-01,
		4.91E-05,
		3.30E+01,
		5.49E+00,
		1.09E+00,
		3.01E-01,
		2.15E-02,
		2.10E-03,
		7.11E-05,
		4.45E-07,
		1.88E-08,
		3.18E-03,
		5.10E-01
	};
	
	string sigLine = sig_samples[signal_num];
	Double_t sigWeight = sig_weights[signal_num];
	
	TH1* sig_hist = new TH1F("sig_hist", TString("Phi Between Met and Photon for Signal Events [") + TString(sigLine) + TString("]"), 70, 0, 3.5); 
	TH1* bg_hist = new TH1F("bg_hist", TString("Phi Between Met and Photon for Background Events [") + TString(sigLine) + TString("]"), 70, 0, 3.5);
	
	TH1* signifhist = new TH1F("signifhist", TString(" "), 24, 2, 3.2);

	cout<<"The selected signal sample is: "<<TString(sigLine)<<" "<<endl;
	
	TString sigFilename = TString("monoph-2013-July9_") + TString(sigLine) + TString("_noskim.root");
	
	MitGPTree sigEvent;
	sigEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+sigFilename, 0);
	sigEvent.InitTree(0);
	
	int nDataSig=sigEvent.tree_->GetEntries();
	
	for (int evt=0; evt<nDataSig; ++evt) 
	{
		sigEvent.tree_->GetEntry(evt);
		Double_t pho1Pt = sigEvent.pho1_.Pt();
		Double_t jet1Pt = sigEvent.jet1_.Pt();
		Double_t met = sigEvent.met_;
		Double_t metPhi = sigEvent.metPhi_;
		Double_t nphotons = sigEvent.nphotons_;
		Double_t ncosmics = sigEvent.ncosmics_;
		Double_t phoPassEleVeto = sigEvent.phoPassEleVeto_a1_;
		Double_t jet1Pt = sigEvent.jet1_.Pt();
		Double_t nlep = sigEvent.nlep_;
		Double_t pho1Eta = sigEvent.pho1_.Eta();
		Double_t pho1Phi = sigEvent.pho1_.Phi();
		Double_t phoIsTrigger = sigEvent.phoIsTrigger_a1_;
		Double_t phoLeadTimeSpan = sigEvent.phoLeadTimeSpan_a1_;
		Double_t phoCoviEtaiEta = sigEvent.phoCoviEtaiEta_a1_;
		Double_t phoCoviPhiiPhi = sigEvent.phoCoviPhiiPhi_a1_;
		Double_t phoMipIsHalo = sigEvent.phoMipIsHalo_a1_;
		Double_t phoSeedTime = sigEvent.phoSeedTime_a1_;
		
		Double_t deltaPhi = TMath::ACos(TMath::Cos(pho1Phi - metPhi));
		
		if (TMath::Abs(pho1Eta)<1.479 && met>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==0 && ncosmics==0 && jet1Pt<100) sig_hist->Fill(deltaPhi, sigWeight);			
	}

	for (int i=0; i<24; i++)
	{
		string bgLine = bg_samples[i];
		Double_t bgWeight = bg_weights[i];

		TString bgFilename = TString("monoph-2013-July9_") + TString(bgLine) + TString("_noskim.root");

		MitGPTree bgEvent;
		bgEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/") + bgFilename, 0);
		bgEvent.InitTree(0);
		
		int nDataBg=bgEvent.tree_->GetEntries();
		
		for (int evt=0; evt<nDataBg; ++evt) 
		{
			bgEvent.tree_->GetEntry(evt);
			Double_t pho1Pt = bgEvent.pho1_.Pt();
			Double_t jet1Pt = bgEvent.jet1_.Pt();
			Double_t met = bgEvent.met_;
			Double_t metPhi = bgEvent.metPhi_;
			Double_t nphotons = bgEvent.nphotons_;
			Double_t ncosmics = bgEvent.ncosmics_;
			Double_t phoPassEleVeto = bgEvent.phoPassEleVeto_a1_;
			Double_t jet1Pt = bgEvent.jet1_.Pt();
			Double_t nlep = bgEvent.nlep_;
			Double_t pho1Eta = bgEvent.pho1_.Eta();
			Double_t pho1Phi = bgEvent.pho1_.Phi();
			Double_t phoIsTrigger = bgEvent.phoIsTrigger_a1_;
			Double_t phoLeadTimeSpan = bgEvent.phoLeadTimeSpan_a1_;
			Double_t phoCoviEtaiEta = bgEvent.phoCoviEtaiEta_a1_;
			Double_t phoCoviPhiiPhi = bgEvent.phoCoviPhiiPhi_a1_;
			Double_t phoMipIsHalo = bgEvent.phoMipIsHalo_a1_;
			Double_t phoSeedTime = bgEvent.phoSeedTime_a1_;
			
			Double_t deltaPhi = TMath::ACos(TMath::Cos(pho1Phi - metPhi));
			
			if (TMath::Abs(pho1Eta)<1.479 && met>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==0 && ncosmics==0 && jet1Pt<100) bg_hist->Fill(deltaPhi, bgWeight);			
		}	
	}
	
	Double_t sigTot = sig_hist->Integral();
	Double_t bgTot = bg_hist->Integral();
	
	Double_t signifTot = sigTot/TMath::Sqrt(bgTot);
	Double_t signifTarget = 0.25*signifTot;
	
	Int_t reduced=0;

// The strange limits in the for loop are a messy hack to get the integral to work out. Since hist->Integral(a,b)
// integrates from BIN a to BIN b (bin numbers instead of values of the x-axis), I choose to work from bin 64
// (an angle cut of 3.2) to bin 40 (an angle cut of 2).

	for (int i=64; i>=40; i--)
	{
		Double_t sigcount = sig_hist->Integral(0, i);
		Double_t bgcount = bg_hist->Integral(0, i);
		if (bgcount>0)
		{
			Double_t significance = sigcount/TMath::Sqrt(bgcount);
//			cout<<"At a phi cut of "<<sig_hist->GetBinCenter(i)+0.025<<" the significance is "<<significance<<", which is a fraction "<<significance/signifTot<<" of max."<<endl;
			if (significance<signifTarget && reduced==0)
			{
				cout<<"The signal significance is reduced to one-fourth of its uncut value at a phi cut of "<<sig_hist->GetBinCenter(i)+0.025<<endl;
				reduced=1;
			}
			signifhist->Fill(sig_hist->GetBinCenter(i), significance);
		}
	}
	
	signifhist->SetStats(kFALSE);
	signifhist->GetXaxis()->SetTitle("Maximum #Delta#phi");
	signifhist->GetYaxis()->SetTitle("Signal Significance");
	signifhist->Draw();
	c1->Print(TString("Blinding_") + TString(sigLine) + TString(".pdf"));
}