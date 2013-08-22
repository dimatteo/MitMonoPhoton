// Znunu_estimation.C
//
// This macro performs the leading order data-driven estimate of the number of Z->Nu nu events in the control region,
// printing this result as well as the leading order Monte Carlo estimate. It uses the finished data files in my
// scratch area and also requires MitGPTree.h to be in the same directory as this file. At the end, this macro
// also saves histograms of the dimuon and dielectron invariant masses for events passing cuts in data and MC.
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

void Znunu_estimation()
{
	gSystem->CompileMacro("MitGPTree.h");
	
	TH1* dielectron_mass_MC = new TH1F("dielectron_mass_MC", "Mass of Dielectron Parent in Zgllg Mc (After cuts)", 200, 0, 200);
	TH1* dimuon_mass_MC = new TH1F("dimuon_mass_MC", "Mass of Dimuon Parent in Zgllg Mc (After cuts)", 200, 0, 200);
	TH1* dielectron_mass_data = new TH1F("dielectron_mass_data", "Mass of Dielectron Parent in Data (After Cuts)", 20, 0, 200);
	TH1* dimuon_mass_data = new TH1F("dimuon_mass_data", "Mass of Dimuon Parent in Data (After Cuts)", 200, 0, 200);

	// First we find the acceptance times efficiency for Z->nu nu using the Monte Carlo.
	
	TString filename = TString("monoph-2013-July9_s12-zgptg130-v7c_noskim.root");
	TFile ZnunuFile(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, "OLD");
	TH1F *h1 = new TH1F("h1", "Events", 1, -0.5, 0.5);
	h1=(TH1F*)ZnunuFile.FindObjectAny("hDEvents");
	Double_t ZnunuAllEvents = h1->GetBinContent(1); // total processed events
	ZnunuFile.Close();
	
	MitGPTree ZnunuEvent;
	ZnunuEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, 0);
	ZnunuEvent.InitTree(0);

	Double_t ZnunuCtr=0;

	int nData=ZnunuEvent.tree_->GetEntries();
	
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
		
		if (TMath::Abs(pho1Eta)<1.479 && met>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==0 && ncosmics==0 && jet1Pt<100 && TMath::Cos(pho1Phi - metPhi)>-0.97 ) ZnunuCtr++;
	}
	
	Double_t ZnunuAccxEff = ZnunuCtr/ZnunuAllEvents;
	
	cout<<"\n"<<endl;
	cout<<"At leading order, the MC Znunu estimate is "<<0.074*19500*ZnunuAccxEff<<" events."<<endl;
	cout<<"\n"<<endl;
	
	// Next we get the Zll acceptance times efficiency.
	
	TString filename = TString("monoph-2013-July9_s12-zgllgptg130-v7a_noskim.root");
	TFile ZllFile(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, "OLD");
	TH1F *h1 = new TH1F("h1", "Events", 1, -0.5, 0.5);
	h1=(TH1F*)ZllFile.FindObjectAny("hDEvents");
	Double_t ZllAllEvents = h1->GetBinContent(1); // total processed events
	ZllFile.Close();	
	
	MitGPTree ZllEvent;
	ZllEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, 1);
	ZllEvent.InitTree(0);
	
	Double_t ZllCtr_ee=0;
	Double_t ZllCtr_mm=0;
	
	int nData=ZllEvent.tree_->GetEntries();
	
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
		
		Double_t newMetX = ZllEvent.lep1_.Px() + ZllEvent.lep2_.Px() + met*TMath::Cos(metPhi);
		Double_t newMetY = ZllEvent.lep1_.Py() + ZllEvent.lep2_.Py() + met*TMath::Sin(metPhi);
		Double_t correctedMet = TMath::Sqrt(newMetX**2 + newMetY**2);
		Double_t correctedMetPhi = TMath::ATan2(newMetY, newMetX);
		
		Double_t ll_mass = TMath::Sqrt( (lep1E+lep2E)**2 - (lep1Px+lep2Px)**2 - (lep1Py+lep2Py)**2-(lep1Pz+lep2Pz)**2);

		if (TMath::Abs(pho1Eta)<1.479 && correctedMet>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==2)
		{
			if (lep1M<0.05 && lep2M<0.05)
			{
				dielectron_mass_MC->Fill(ll_mass);
				if (ll_mass>75 && ll_mass<105) ZllCtr_ee++;
			}
			else if (lep1M>0.05 && lep2M>0.05)
			{
				dimuon_mass_MC->Fill(ll_mass);
				if (ll_mass>75 && ll_mass<105) ZllCtr_mm++;
			}
		}
	}	
	
	Double_t ZmmAccxEff = (ZllCtr_mm)/ZllAllEvents;
	Double_t ZeeAccxEff = (ZllCtr_ee)/ZllAllEvents;
	Double_t ZllAccxEff = (ZllCtr_mm + ZllCtr_ee)/ZllAllEvents;
	
	// Now we count dilepton events in data which pass all our cuts.
	
	string data_samples[] = {"monoph-2013-July9_r12a-phh-j22-v1_noskim.root", "monoph-2013-July9_r12a-pho-j22-v1_noskim.root", "monoph-2013-July9_r12b-phh-j22-v1_noskim.root", "monoph-2013-July9_r12b-sph-j22-v1_noskim.root", "monoph-2013-July9_r12c-phh-j22-v1_noskim.root", "monoph-2013-July9_r12c-sph-j22-v1_noskim.root", "monoph-2013-July9_r12d-phh-j22-v1_noskim.root", "monoph-2013-July9_r12d-sph-j22-v1_noskim.root"};
	Double_t dataCount_ee=0;
	Double_t dataCount_mm=0;

	for (int i=0; i<8; i++)
	{
		TString filename = data_samples[i];

		MitGPTree dataEvent;
		
		dataEvent.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/")+ filename, 1);
		dataEvent.InitTree(0);
		
		int nData=dataEvent.tree_->GetEntries(); 
		
		for (int evt=0; evt<nData; ++evt) 
		{
			dataEvent.tree_->GetEntry(evt);
			Double_t pho1Pt = dataEvent.pho1_.Et();
			Double_t jet1Pt = dataEvent.jet1_.Pt();
			Double_t met = dataEvent.met_;
			Double_t metPhi = dataEvent.metPhi_;
			Double_t nphotons = dataEvent.nphotons_;
			Double_t ncosmics = dataEvent.ncosmics_;
			Double_t phoPassEleVeto = dataEvent.phoPassEleVeto_a1_;
			Double_t nlep = dataEvent.nlep_;
			Double_t pho1Eta = dataEvent.pho1_.Eta();
			Double_t pho1Phi = dataEvent.pho1_.Phi();
			Double_t phoIsTrigger = dataEvent.phoIsTrigger_a1_;
			Double_t phoLeadTimeSpan = dataEvent.phoLeadTimeSpan_a1_;
			Double_t phoCoviEtaiEta = dataEvent.phoCoviEtaiEta_a1_;
			Double_t phoCoviPhiiPhi = dataEvent.phoCoviPhiiPhi_a1_;
			Double_t phoMipIsHalo = dataEvent.phoMipIsHalo_a1_;
			Double_t phoSeedTime = dataEvent.phoSeedTime_a1_;
			Double_t lep1M = dataEvent.lep1_.M();
			Double_t lep1E = dataEvent.lep1_.E();
			Double_t lep1Px = dataEvent.lep1_.Px();
			Double_t lep1Py = dataEvent.lep1_.Py();
			Double_t lep1Pz = dataEvent.lep1_.Pz();	
			Double_t lep2M = dataEvent.lep2_.M();
			Double_t lep2E = dataEvent.lep2_.E();
			Double_t lep2Px = dataEvent.lep2_.Px();
			Double_t lep2Py = dataEvent.lep2_.Py();
			Double_t lep2Pz = dataEvent.lep2_.Pz();
			
			Double_t newMetX = dataEvent.lep1_.Px() + dataEvent.lep2_.Px() + met*TMath::Cos(metPhi);
			Double_t newMetY = dataEvent.lep1_.Py() + dataEvent.lep2_.Py() + met*TMath::Sin(metPhi);
			Double_t correctedMet = TMath::Sqrt(newMetX**2 + newMetY**2);
			Double_t correctedMetPhi = TMath::ATan2(newMetY, newMetX);

			Double_t ll_mass = TMath::Sqrt( (lep1E+lep2E)**2 - (lep1Px+lep2Px)**2 - (lep1Py+lep2Py)**2-(lep1Pz+lep2Pz)**2);
			
			if (TMath::Abs(pho1Eta)<1.479 && correctedMet>140 && pho1Pt>160 && nphotons>0 && phoPassEleVeto>0 && phoIsTrigger==1 && TMath::Abs(phoLeadTimeSpan) < 8. && phoCoviEtaiEta > 0.001 && phoCoviPhiiPhi > 0.001 && phoMipIsHalo == 0 && phoSeedTime > -1.5 && nlep==2)
			{
				if (lep1M<0.05 && lep2M<0.05){
					if (ll_mass>75 && ll_mass<105) dataCount_ee++;
					dielectron_mass_data->Fill(ll_mass);
				}
				else if (lep1M>0.05 && lep2M>0.05)
				{
					if (ll_mass>75 && ll_mass<105) dataCount_mm++;
					dimuon_mass_data->Fill(ll_mass);
				}
			}
		}
	}
	
	cout<<"Data has in total "<<dataCount_mm<<" dimuon events and "<<dataCount_ee<<" dielectron events passing cuts."<<endl;
	
	Double_t muon_uncertainty = Double_t(1)/TMath::Sqrt(dataCount_mm);
	Double_t electron_uncertainty = Double_t(1)/TMath::Sqrt(dataCount_ee);
	Double_t lepton_uncertainty = Double_t(1)/TMath::Sqrt(dataCount_ee+dataCount_mm);
	
	// Here are the branching ratios from pdg, and the ratio of LO cross-sections as suggested by Guillelmo.
	
	Double_t ZnunuBranch = 20; 
	Double_t ZeeBranch = 3.363;
	Double_t ZmmBranch = 3.366;
	Double_t ZttBranch = 3.370;
	Double_t ZllBranch = ZeeBranch+ZmmBranch+ZttBranch;

	// Note: Before Guillelmo's suggestion, we used ZnunuBranch/ZllBranch instead of guillelmo_ratio, leading
	// to an estimate that was about 28% higher.
	
	Double_t guillelmo_ratio = 1.548;
	
	//Finally, here we estimate the Z->nu nu number using the data counts and appropriate ratios.

	Double_t ZnunuPrediction_mm = (dataCount_mm)*(ZnunuAccxEff/ZmmAccxEff)*(guillelmo_ratio);
	Double_t ZnunuPrediction_ee = (dataCount_ee)*(ZnunuAccxEff/ZeeAccxEff)*(guillelmo_ratio);
	Double_t ZnunuPrediction_ll = (dataCount_mm + dataCount_ee)*(ZnunuAccxEff/ZllAccxEff)*(guillelmo_ratio);
	
	cout<<"\n"<<endl;
	cout<<"The data-driven Z->nu nu prediction is: "<<endl;
	cout<<"     "<<ZnunuPrediction_mm<<" +- "<<muon_uncertainty*ZnunuPrediction_mm<<" from just muons, "<<endl;
	cout<<"     "<<ZnunuPrediction_ee<<" +- "<<electron_uncertainty*ZnunuPrediction_ee<<" from just electrons, and "<<endl;
	cout<<"     "<<ZnunuPrediction_ll<<" +- "<<lepton_uncertainty*ZnunuPrediction_ll<<" from all leptons. "<<endl;
	cout<<"Here the uncertainties are only statistical from counting. "<<endl;
	cout<<"\n"<<endl;
	
	dielectron_mass_MC->Draw();
	c1->Print("dielectron_mass_MC.pdf");
	
	dimuon_mass_MC->Draw();
	c1->Print("dimuon_mass_MC.pdf");

	dielectron_mass_data->Draw();
	c1->Print("dielectron_mass_data.pdf");

	dimuon_mass_data->Draw();
	c1->Print("dimuon_mass_data.pdf");
	
}
