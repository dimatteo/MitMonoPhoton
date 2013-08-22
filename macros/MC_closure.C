// MC_closure.C
//
// This macro performs a closure test of the data-driven method used to
// estimate the expected yield of Z->nu nu events, which is implemented
// in Znunu_estimation.C. Here we generate 5000 pseudoexperiments, each
// with the expected number of Z->mu mu events with invariant masses 
// distributed according to the Monte Carlo, and in each we count
// how many events pass the invariant mass cut. In each pseudoexperiment
// we then compute the estimated Z->nu nu yield and calculate the fractional
// difference from the MC estimate. Finally we histogram all these fractional
// differences to find the overall bias of this method.
//
// Christian Ferko

#include "TMath.h"
#include "TTree.h"
#include "MitGPTree.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void MC_closure()
{

	TStyle *MitStyle = gStyle;// new TStyle("MIT-Style","The Perfect Style for Plots ;-)");
	MitStyle->SetCanvasColor     (0);
	MitStyle->SetPadColor       (0);
	MitStyle->SetOptStat    (0);
   
	//Now begins the actual macro

	gSystem->CompileMacro("MitGPTree.h");

    TH1* Zll_dimuon_mass = new TH1F("Zll_dimuon_mass", "Dimuon Mass in Zll Events Passing Selection", 2000, 0, 200);	
	TH1* frac_error_hist_tight = new TH1F("frac_error_hist_tight", "Fractional Difference between Data-Driven and Monte Carlo (75<M_{#mu#mu}<105)", 20, -0.3, 0.3);
	TH1* frac_error_hist_loose = new TH1F("frac_error_hist_loose", "Fractional Difference between Data-Driven and Monte Carlo (60<M_{#mu#mu}<120)", 20, -0.3, 0.3);
	
	TH2* ThreeBody_vs_TwoBody = new TH2F("ThreeBody_vs_TwoBody", "3-Body Mass vs. 2-Body Mass for Events Passing Cuts in Zgllg MC", 1000, 0, 1000, 1500, 0, 1500);

	// First apply selection to Zll to find the expected number of events and invariant mass distribution.
	
	TString filename = TString("monoph-2013-July9_s12-zgllgptg130-v7a_noskim.root");
	TFile ZllFile(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, "OLD");
	TH1F *h1 = new TH1F("h1", "Events", 1, -0.5, 0.5);
	h1=(TH1F*)ZllFile.FindObjectAny("hDEvents");
	Double_t ZllAllEvents = h1->GetBinContent(1); // all processed events
	ZllFile.Close();	
	
	MitGPTree ZllEvent;
	ZllEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, 1);
	ZllEvent.InitTree(0);
	
	Double_t ZllCtr_mm = 0;
	Double_t ZllCtr_tight=0; // tight mass range is 75<m(mu mu)<105
	Double_t ZllCtr_loose=0; // loose mass range is 60<m(mu mu)<120
	
	int nData=ZllEvent.tree_->GetEntries();
	Double_t Zll_efficiency = Double_t(nData)/ZllAllEvents;
	
	for (int evt=0; evt<nData; ++evt) 
	{
		ZllEvent.tree_->GetEntry(evt);
		Double_t pho1Pt = ZllEvent.pho1_.Pt();
		Double_t jet1Pt = ZllEvent.jet1_.Pt();
		Double_t met = ZllEvent.met_;
		Double_t metPhi = ZllEvent.metPhi_;
		Double_t nphotons = ZllEvent.nphotons_;
		Double_t ncosmics = ZllEvent.ncosmics_;
		Double_t phoPassEleVeto = ZllEvent.phoPassEleVeto_a1_;
		Double_t nlep = ZllEvent.nlep_;
		Double_t pho1Eta = ZllEvent.pho1_.Eta();
		Double_t pho1Phi = ZllEvent.pho1_.Phi();
		Double_t phoIsTrigger = ZllEvent.phoIsTrigger_a1_;
		Double_t phoLeadTimeSpan = ZllEvent.phoLeadTimeSpan_a1_;
		Double_t phoCoviEtaiEta = ZllEvent.phoCoviEtaiEta_a1_;
		Double_t phoCoviPhiiPhi = ZllEvent.phoCoviPhiiPhi_a1_;
		Double_t phoMipIsHalo = ZllEvent.phoMipIsHalo_a1_;
		Double_t phoSeedTime = ZllEvent.phoSeedTime_a1_;
		Double_t lep1M = ZllEvent.lep1_.M();
		Double_t lep1E = ZllEvent.lep1_.E();
		Double_t lep1Px = ZllEvent.lep1_.Px();
		Double_t lep1Py = ZllEvent.lep1_.Py();
		Double_t lep1Pz = ZllEvent.lep1_.Pz();	
		Double_t lep2M = ZllEvent.lep2_.M();
		Double_t lep2E = ZllEvent.lep2_.E();
		Double_t lep2Px = ZllEvent.lep2_.Px();
		Double_t lep2Py = ZllEvent.lep2_.Py();
		Double_t lep2Pz = ZllEvent.lep2_.Pz();
		Double_t pho1E = ZllEvent.pho1_.E();
		Double_t pho1Px = ZllEvent.pho1_.Px();
		Double_t pho1Py = ZllEvent.pho1_.Py();
		Double_t pho1Pz = ZllEvent.pho1_.Pz();	
		
		Double_t newMetX = ZllEvent.lep1_.Px() + ZllEvent.lep2_.Px() + met*TMath::Cos(metPhi);
		Double_t newMetY = ZllEvent.lep1_.Py() + ZllEvent.lep2_.Py() + met*TMath::Sin(metPhi);
		Double_t correctedMet = TMath::Sqrt(newMetX**2 + newMetY**2);
		Double_t correctedMetPhi = TMath::ATan2(newMetY, newMetX);
		
		Double_t ll_mass = TMath::Sqrt( (lep1E+lep2E)**2 - (lep1Px+lep2Px)**2 - (lep1Py+lep2Py)**2-(lep1Pz+lep2Pz)**2);
		
		Double_t llg_mass = TMath::Sqrt( (lep1E+lep2E+pho1E)**2 - (lep1Px+lep2Px+pho1Px)**2 - (lep1Py+lep2Py+pho1Py)**2-(lep1Pz+lep2Pz+pho1Pz)**2);

		if (TMath::Abs(pho1Eta)<1.479 && correctedMet>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==2 && lep1M>0.05 && lep2M>0.05 && jet1Pt<100)
		{
			ZllCtr_mm++;
            Zll_dimuon_mass->Fill(ll_mass);
			if (ll_mass>60 && ll_mass<120 ) ZllCtr_loose++;
			if (ll_mass>75 && ll_mass<105 ) ZllCtr_tight++;
			ThreeBody_vs_TwoBody->Fill(ll_mass, llg_mass);
		}
	}	
	
	Double_t Zll_Prediction = (ZllCtr_mm*4.78E-02*19500)/(ZllAllEvents); // this uses the LO cross section from prep

	Double_t Zmm_acceptance_loose = ZllCtr_loose/(Double_t(nData));	
	Double_t Zmm_acceptance_tight = ZllCtr_tight/(Double_t(nData));
	
	// Next get the acceptance and efficiency for Z->nu nu, along with the MC estimate.
	
	TString filename = TString("monoph-2013-July9_s12-zgptg130-v7c_noskim.root");
	TFile ZnunuFile(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, "OLD");
	TH1F *h1 = new TH1F("h1", "Events", 1, -0.5, 0.5);
	h1=(TH1F*)ZnunuFile.FindObjectAny("hDEvents");
	Double_t ZnunuAllEvents = h1->GetBinContent(1); // all processed events
	ZnunuFile.Close();

	MitGPTree ZnunuEvent;
	ZnunuEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, 0);
	ZnunuEvent.InitTree(0);

	Double_t ZnunuCtr=0;

	int nData=ZnunuEvent.tree_->GetEntries();
	
	Double_t Znunu_efficiency = Double_t(nData)/ZnunuAllEvents;
	
	for (int evt=0; evt<nData; ++evt) 
	{
		ZnunuEvent.tree_->GetEntry(evt);
		Double_t pho1Pt = ZnunuEvent.pho1_.Pt();
		Double_t jet1Pt = ZnunuEvent.jet1_.Pt();
		Double_t met = ZnunuEvent.met_;
		Double_t metPhi = ZnunuEvent.metPhi_;
		Double_t nphotons = ZnunuEvent.nphotons_;
		Double_t ncosmics = ZnunuEvent.ncosmics_;
		Double_t phoPassEleVeto = ZnunuEvent.phoPassEleVeto_a1_;
		Double_t nlep = ZnunuEvent.nlep_;
		Double_t pho1Eta = ZnunuEvent.pho1_.Eta();
		Double_t pho1Phi = ZnunuEvent.pho1_.Phi();
		Double_t phoIsTrigger = ZnunuEvent.phoIsTrigger_a1_;
		Double_t phoLeadTimeSpan = ZnunuEvent.phoLeadTimeSpan_a1_;
		Double_t phoCoviEtaiEta = ZnunuEvent.phoCoviEtaiEta_a1_;
		Double_t phoCoviPhiiPhi = ZnunuEvent.phoCoviPhiiPhi_a1_;
		Double_t phoMipIsHalo = ZnunuEvent.phoMipIsHalo_a1_;
		Double_t phoSeedTime = ZnunuEvent.phoSeedTime_a1_;
		
		if (TMath::Abs(pho1Eta)<1.479 && met>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==0 && jet1Pt<100 && ncosmics==0 && TMath::Cos(pho1Phi-metPhi)>-0.97) 
		{
			ZnunuCtr++;
		}
	}
	
	Double_t Znunu_acceptance = ZnunuCtr/Double_t(nData);
	Double_t ZnunuPrediction_MC = (Znunu_acceptance * Znunu_efficiency * 7.40E-02 * 19500); // cross section from prep
	
	// Now we generate 5000 pseudodatasets and find the DD estimate in each case.
	
	for (int i=0; i<5000; i++)
	{
		Double_t loose_data_count=0;
		Double_t tight_data_count=0;

		for (int n=0; n<TMath::Nint(Zll_Prediction); n++)
		{
			Double_t pseudo_dimuon_mass = Zll_dimuon_mass->GetRandom();
			if (pseudo_dimuon_mass>60 && pseudo_dimuon_mass<120) loose_data_count++;
			if (pseudo_dimuon_mass>75 && pseudo_dimuon_mass<105) tight_data_count++;
		}
		
		// Here are the branching ratios from pdg
		
		Double_t ZnunuBranch = 20; 
		Double_t ZeeBranch = 3.363;
		Double_t ZmmBranch = 3.366;
		Double_t ZttBranch = 3.370;
		Double_t ZllBranch = ZeeBranch+ZmmBranch+ZttBranch;
		
		// Before Guillelmo's suggest, we used (ZnunuBranch)/(ZllBranch) instead of guillelmo_ratio
		
		Double_t guillelmo_ratio = 1.548;
		
		// Finally, perform the Z nu nu estimation and find error with MC method

		Double_t ZnunuPrediction_DD_loose = (loose_data_count)*( (Znunu_efficiency*Znunu_acceptance)/(Zll_efficiency*Zmm_acceptance_loose) )*(guillelmo_ratio);
		Double_t ZnunuPrediction_DD_tight = (tight_data_count)*( (Znunu_efficiency*Znunu_acceptance)/(Zll_efficiency*Zmm_acceptance_tight) )*(guillelmo_ratio);
	
		frac_error_hist_loose->Fill( (ZnunuPrediction_DD_loose - ZnunuPrediction_MC)/ZnunuPrediction_MC );
		frac_error_hist_tight->Fill( (ZnunuPrediction_DD_tight - ZnunuPrediction_MC)/ZnunuPrediction_MC );
	}

	frac_error_hist_loose->GetXaxis()->SetTitle("(N_{DD} - N_{MC})/ N_{MC}");
	frac_error_hist_loose->GetYaxis()->SetTitle("Pseudoexperiments");
	frac_error_hist_loose->Draw();
	frac_error_hist_loose->Fit("gaus");
	c1->Print("frac_error_hist_loose.pdf");
	
	cout<<"\n"<<endl;
	
	frac_error_hist_tight->GetXaxis()->SetTitle("(N_{DD} - N_{MC})/ N_{MC}");
	frac_error_hist_tight->GetYaxis()->SetTitle("Pseudoexperiments");
	frac_error_hist_tight->Draw();
	frac_error_hist_tight->Fit("gaus");
	c1->Print("frac_error_hist_tight.pdf");
	
	ThreeBody_vs_TwoBody->GetXaxis()->SetTitle("M_{#mu #mu}");
	ThreeBody_vs_TwoBody->GetYaxis()->SetTitle("M_{#mu #mu #gamma}");
	ThreeBody_vs_TwoBody->Draw();
	c1->Print("ThreeBody_vs_TwoBody.pdf");
	
}
