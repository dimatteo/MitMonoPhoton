#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "MitMonoPhoton/Core/MitGPTree.h"
#include "TLorentzVector.h"
#include "TFractionFitter.h"
#include "TH1D.h"
#include "THStack.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TPad.h"
#include "TMath.h"
#include "math.h"

Double_t min(Double_t num1, Double_t num2) {

  Double_t lesser;
  if ( num1 > num2 ) lesser = num2;
  if ( num1 < num2 ) lesser = num1;
  if ( num1 == num2 ) lesser = num1;
  return lesser;
}

Double_t deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = TMath::Abs(phi1 - phi2);
  if ( dPhi > 3.14159 ) dPhi = 2*3.14159 - TMath::Abs(phi1 - phi2);
  return dPhi;
}

void fakeMonoPhotons() {

  //*************************************************************************************************
  // Loose Photon Definition
  //*************************************************************************************************
  // Photon energy > x
  Double_t x_phoE = 160;

  // Low Missing energy < x
  Double_t x_met = 30; 

  // High MET defined as > x
  Double_t x_highMet = 140;

  // |Photon Eta| < x
  Double_t x_eta = 1.479; 
  
  // E_Had / E_Em < x
  Double_t x_hadEm = 0.05;
  
  //  x < sigma_ietai_eta < y
  Double_t x_covEta = 0.001;
  Double_t y_covEta = 0.013;

  // x < sigma_iphii_phi
  Double_t x_covPhi = 0.001;

  // |Lead time span| < x
  Double_t x_leadTime = 8.0;

  // No Pixel Seed

  // Is HLT Trigger Object

  // Not Electron

  // Not Beam Halo Tagged

  // Cosmic Veto

  // Seed Time > x
  Double_t x_seedTime = -1.5;

  // Jet Veto
  Double_t x_jetPt = 1000000.;

  // Iso < Min( x*(den Iso), y*(photon pt) ) 
  Double_t x_looseIso = 5.0; 
  // Double_t y_looseIso = 0.2;

  // QCD Sideband Iso x*x_iso3 < Iso3 < y*x_iso3
  Double_t x_sideband = 0.5;
  Double_t y_sideband = 2;

  //*************************************************************************************************
  // Numerator Definition (loose photon + ..)
  //*************************************************************************************************
  // Iso < x, if Photon R9 >= y
  Double_t x_iso1high = 6.0;
  Double_t x_iso2high = 10.0;
  Double_t x_iso3high = 3.8;
  Double_t y_isoR9    = 0.94;

  // Iso < x, if Photon R9 < y
  Double_t x_iso1low  = 4.7;
  Double_t x_iso2low  = 6.5;
  Double_t x_iso3low  = 2.5;

  Double_t x_iso1;
  Double_t x_iso2;
  Double_t x_iso3;
    
  //*************************************************************************************************
  // Denominator Definition (loose photon + ..)
  //*************************************************************************************************
  // One inverted isolation requirement (Iso > x + y*(photon pt))

  //*************************************************************************************************
  // Output objects
  //*************************************************************************************************
  const Int_t nBins = 6;
  Double_t xBins[nBins+1] = {
    145,
    160,
    190,
    250,
    400,
    700,
    1000,
  };
  TH1D    *numerators    = new TH1D("numerators","Numerators",nBins,xBins);
  TH1D    *denominators  = new TH1D("denominators","Denominators",nBins,xBins);
  TH1D    *looseSelection= new TH1D("looseSelection","Loose Selection",nBins,xBins);
  TH1D    *loosePhotons  = new TH1D("loosePhotons","Loose Photons",1000,0,1000);
  TH1D    *data145_160   = new TH1D("data145_160","Data",200,0,0.04);
  TH1D    *data160_190   = new TH1D("data160_190","Data",200,0,0.04);
  TH1D    *data190_250   = new TH1D("data190_250","Data",200,0,0.04);
  TH1D    *data250_400   = new TH1D("data250_400","Data",200,0,0.04);
  TH1D    *data400_700   = new TH1D("data400_700","Data",200,0,0.04);
  TH1D    *data700_1000  = new TH1D("data700_1000","Data",200,0,0.04);
  TH1D    *fake145_160   = new TH1D("fake145_160","Fake",200,0,0.04);
  TH1D    *fake160_190   = new TH1D("fake160_190","Fake",200,0,0.04);
  TH1D    *fake190_250   = new TH1D("fake190_250","Fake",200,0,0.04);
  TH1D    *fake250_400   = new TH1D("fake250_400","Fake",200,0,0.04);
  TH1D    *fake400_700   = new TH1D("fake400_700","Fake",200,0,0.04);
  TH1D    *fake700_1000  = new TH1D("fake700_1000","Fake",200,0,0.04);
  TH1D    *mc145_160     = new TH1D("mc145_160","MC",200,0,0.04);
  TH1D    *mc160_190     = new TH1D("mc160_190","MC",200,0,0.04);
  TH1D    *mc190_250     = new TH1D("mc190_250","MC",200,0,0.04);
  TH1D    *mc250_400     = new TH1D("mc250_400","MC",200,0,0.04);
  TH1D    *mc400_700     = new TH1D("mc400_700","MC",200,0,0.04);
  TH1D    *mc700_1000    = new TH1D("mc700_1000","MC",200,0,0.04);
  TH1D    *dataHist;
  TH1D    *fakeHist;
  TH1D    *mcHist;

  Double_t loosePhMet[nBins];
  Double_t num[nBins];
  Double_t den[nBins];
  Int_t    whichBin;

  for ( Int_t i=0; i<nBins; i++ ) {
    loosePhMet[i] = 0;
    num[i]        = 0;
    den[i]        = 0;
  }

  //*************************************************************************************************
  // Loop
  //*************************************************************************************************
  // Get list of files
  std::string dataString("r12");
  std::string mcString("s12");
  std::string buff;
  std::ifstream file ("files.txt");
  while ( file.good() ) {

    // Load each file name
    std::getline (file,buff);

    // Stop if end of file
    if( buff == "" ) break;

    // Initialize file
    MitGPTree event;
    TString filename = TString("monoph-2013-July9_") + TString(buff) + TString("_noskim.root");
    event.LoadTree(TString("/scratch/cferko/hist/monoph-2013-July9/merged/") + filename,2);
    event.InitTree(0);

    // Loop over events
    Int_t nEntries = event.tree_->GetEntries();
    for ( Int_t i=0; i<nEntries; i++ ) {

      event.tree_->GetEntry(i);

      // Loop over photons of each event 
      for ( Int_t j=1; j<5; j++ ) {
     
	//*******************************************************************************************
	// Set Variables
	//*******************************************************************************************
	// Define
	Double_t phoPt;
	Double_t phoE;
	Double_t phoEta;
	Double_t phoPhi;
	Double_t hadEm;
	Double_t covEta;
	Double_t covPhi;
	Double_t leadTime;
	Double_t phoR9;
	Bool_t   pxlSeed;
	Bool_t   hltMatch;
	Bool_t   notEle;
	Bool_t   isHalo;
	Double_t seedTime;
	Int_t    nleps;
	Int_t    ncosmics;
	Double_t jetPt;
	Double_t iso1;
	Double_t iso2;
	Double_t iso3;
	Double_t dPhi;

	if ( j == 1 ) {
	  phoPt    = event.pho1_.Pt();
	  phoE     = event.pho1_.E();
	  phoEta   = event.pho1_.Eta();
	  phoPhi   = event.pho1_.Phi();
	  hadEm    = event.phoHadOverEm_a1_;
	  covEta   = event.phoCoviEtaiEta_a1_;
	  covPhi   = event.phoCoviPhiiPhi_a1_;
	  leadTime = event.phoLeadTimeSpan_a1_;
	  phoR9    = event.phoR9_a1_;
	  pxlSeed  = event.phoHasPixelSeed_a1_;
	  hltMatch = event.phoIsTrigger_a1_;
	  notEle   = event.phoPassEleVeto_a1_;
	  isHalo   = event.phoMipIsHalo_a1_;
	  seedTime = event.phoSeedTime_a1_;
	  iso1     = event.phoCombIso1_a1_;
	  iso2     = event.phoCombIso2_a1_;
	  iso3     = event.phoCombIso3_a1_;
	}     
	if ( j == 2 ) {
	  phoPt    = event.pho2_.Pt();
	  phoE     = event.pho2_.E();
	  phoEta   = event.pho2_.Eta();
	  phoPhi   = event.pho2_.Phi();
	  hadEm    = event.phoHadOverEm_a2_;
	  covEta   = event.phoCoviEtaiEta_a2_;
	  covPhi   = event.phoCoviPhiiPhi_a2_;
	  leadTime = event.phoLeadTimeSpan_a2_;
	  phoR9    = event.phoR9_a2_;
	  pxlSeed  = event.phoHasPixelSeed_a2_;
	  hltMatch = event.phoIsTrigger_a2_;	  
	  notEle   = event.phoPassEleVeto_a2_;
	  isHalo   = event.phoMipIsHalo_a2_;
	  seedTime = event.phoSeedTime_a2_;
	  iso1     = event.phoCombIso1_a2_;
	  iso2     = event.phoCombIso2_a2_;
	  iso3     = event.phoCombIso3_a2_;
	}
	if ( j == 3 ) {
	  phoPt    = event.pho3_.Pt();
	  phoE     = event.pho3_.E();
	  phoEta   = event.pho3_.Eta();
	  phoPhi   = event.pho3_.Phi();
	  hadEm    = event.phoHadOverEm_a3_;
	  covEta   = event.phoCoviEtaiEta_a3_;
	  covPhi   = event.phoCoviPhiiPhi_a3_;
	  leadTime = event.phoLeadTimeSpan_a3_;
	  phoR9    = event.phoR9_a3_;
	  pxlSeed  = event.phoHasPixelSeed_a3_;
	  hltMatch = event.phoIsTrigger_a3_;
	  notEle   = event.phoPassEleVeto_a3_;
	  isHalo   = event.phoMipIsHalo_a3_;
	  seedTime = event.phoSeedTime_a3_;
	  iso1     = event.phoCombIso1_a3_;
	  iso2     = event.phoCombIso2_a3_;
	  iso3     = event.phoCombIso3_a3_;
	}
	if ( j == 4 ) {
	  phoPt    = event.pho4_.Pt();
	  phoE     = event.pho4_.E();
	  phoEta   = event.pho4_.Eta();
	  phoPhi   = event.pho4_.Phi();
	  hadEm    = event.phoHadOverEm_a4_;
	  covEta   = event.phoCoviEtaiEta_a4_;
	  covPhi   = event.phoCoviPhiiPhi_a4_;
	  leadTime = event.phoLeadTimeSpan_a4_;
	  phoR9    = event.phoR9_a4_;
	  pxlSeed  = event.phoHasPixelSeed_a4_;
	  hltMatch = event.phoIsTrigger_a4_;
	  notEle   = event.phoPassEleVeto_a4_;
	  isHalo   = event.phoMipIsHalo_a4_;
	  seedTime = event.phoSeedTime_a4_;
	  iso1     = event.phoCombIso1_a4_;
	  iso2     = event.phoCombIso2_a4_;
	  iso3     = event.phoCombIso3_a4_;
	}

	nleps      = event.nlep_;
	ncosmics   = event.ncosmics_;
	jetPt      = event.jet1_.Pt();
	dPhi       = deltaPhi(phoPhi,event.metPhi_);

	// Determine correct bin
	if ( phoPt < xBins[1] ) {
	  whichBin = 0;
	  dataHist = data145_160;
	  fakeHist = fake145_160;
	  mcHist   = mc145_160;
	}
	if ( xBins[1] <= phoPt && phoPt < xBins[2] ) {
	  whichBin = 1;
	  dataHist = data160_190;
	  fakeHist = fake160_190;
	  mcHist   = mc160_190;
	}
	if ( xBins[2] <= phoPt && phoPt < xBins[3] ) {
	  whichBin = 2;
	  dataHist = data190_250;
	  fakeHist = fake190_250;
	  mcHist   = mc190_250;
	}
	if ( xBins[3] <= phoPt && phoPt < xBins[4] ) {
	  whichBin = 3;
	  dataHist = data250_400;
	  fakeHist = fake250_400;
	  mcHist   = mc250_400;
	}
	if ( xBins[4] <= phoPt && phoPt < xBins[5] ) {
	  whichBin = 4;
	  dataHist = data400_700;
	  fakeHist = fake400_700;
	  mcHist   = mc400_700;
	}
	if ( xBins[5] <= phoPt ) {
	  whichBin = 5;
	  dataHist = data700_1000;
	  fakeHist = fake700_1000;
	  mcHist   = mc700_1000;
	}

	// Which isolation requirements
	if ( phoR9 >= y_isoR9 ) {
	  x_iso1 = x_iso1high;
	  x_iso2 = x_iso2high;
	  x_iso3 = x_iso3high;
	}
	if ( phoR9 < y_isoR9 ) {
	  x_iso1 = x_iso1low;
	  x_iso2 = x_iso2low;
	  x_iso3 = x_iso3low;
	}

	//*******************************************************************************************
	// Cuts
	//*******************************************************************************************
	// Regional cuts
	if ( phoE < x_phoE )                     continue;
	if ( TMath::Abs(phoEta) > x_eta )        continue;
	if ( hadEm > x_hadEm )                   continue;
	if ( covEta < x_covEta )                 continue;
	if ( covPhi < x_covPhi )                 continue;
	if ( TMath::Abs(leadTime) > x_leadTime ) continue;
	if ( !hltMatch )                         continue;
	if ( !notEle )                           continue;
	if ( isHalo )                            continue;
	if ( seedTime < x_seedTime )             continue;
	if ( nleps != 0 )                        continue;
	if ( ncosmics != 0 )                     continue;


	// Check whether reading data
	if ( buff.compare(0, dataString.length(), dataString) == 0 ) {

	  // Check whether loose photon + high met candidate	  
	  if ( event.met_ > x_met ) {
	    if ( event.met_ > x_highMet )
	      if ( iso1 > x_iso1 || iso2 > x_iso2 || iso3 > x_iso3 )
		if ( iso1 < x_looseIso*x_iso1 )
		  if ( iso2 < x_looseIso*x_iso2 )
		    if ( iso3 < x_looseIso*x_iso3 )
		      if ( covEta > x_covEta && covEta < y_covEta ) {
			loosePhMet[whichBin]++;
			loosePhotons->Fill(phoPt);
		      }
	    continue;
	  }
	  // Numerators
	  if ( iso1 < x_iso1 )
	    if ( iso2 < x_iso2 ) {
	      if ( iso3 < x_iso3 )
		// Jet veto
		if ( jetPt < x_jetPt ) {
		  // // BLIND
		  // if ( TMath::Cos(dPhi) > -0.97 ) {
		    dataHist->Fill(TMath::Abs(covEta));
		    if ( covEta > x_covEta && covEta < y_covEta )
		      num[whichBin]++;
		  // }
		}
	      // Sideband
	      if ( iso3 > x_sideband*x_iso3 && iso3 < y_sideband*x_iso3 )
		fakeHist->Fill(TMath::Abs(covEta));
	    }
	  // Denominators
	  if ( iso1 > x_iso1 || iso2 > x_iso2 || iso3 > x_iso3 )
	    if ( iso1 < x_looseIso*x_iso1 )
	      if ( iso2 < x_looseIso*x_iso2 )
		if ( iso3 < x_looseIso*x_iso3 )
		  if ( covEta > x_covEta && covEta < y_covEta )
		    den[whichBin]++;
	}
	// Check whether reading monte carlo
	if ( buff.compare(0, mcString.length(), mcString) == 0 ) {

	  if ( event.met_ > x_met ) continue;
	  if ( iso1 < x_iso1 )
	    if ( iso2 < x_iso2 )
	      if ( iso3 < x_iso3 )
		mcHist->Fill(TMath::Abs(covEta));
	}
      }
    }
  }

  //*************************************************************************************************
  // Plots for fakePhotonsCalculator.C
  //*************************************************************************************************
  for ( Int_t i=0; i<nBins; i++ ) {
    denominators->Fill(xBins[i],den[i]);
    numerators->Fill(xBins[i],num[i]);
    looseSelection->Fill(xBins[i],loosePhMet[i]);
  }

  TFile *preFit = new TFile("sysBin.root","RECREATE");
  data145_160->Write();  fake145_160->Write();  mc145_160->Write();
  data160_190->Write();  fake160_190->Write();  mc160_190->Write();
  data190_250->Write();  fake190_250->Write();  mc190_250->Write();
  data250_400->Write();  fake250_400->Write();  mc250_400->Write();
  data400_700->Write();  fake400_700->Write();  mc400_700->Write();
  data700_1000->Write(); fake700_1000->Write(); mc700_1000->Write();
  numerators->Write();   denominators->Write(); looseSelection->Write();
  loosePhotons->Write();
  preFit->Close();
}
