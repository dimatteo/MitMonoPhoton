 // $Id $

#include "MitMonoPhoton/SelMods/interface/MonoPhotonAnalysisMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include <sstream>

using namespace mithep;
ClassImp(mithep::MonoPhotonAnalysisMod)

//--------------------------------------------------------------------------------------------------
MonoPhotonAnalysisMod::MonoPhotonAnalysisMod(const char *name, const char *title) : 
  BaseMod(name,title),
  // Define all the Branches to load
  // define all the Branches to load
  fMetBranchName                 ("PFMet"),
  fPhotonsBranchName             (Names::gkPhotonBrn),
  fElectronsBranchName           (Names::gkElectronBrn),
  fMuonsBranchName               (Names::gkMuonBrn),
  fLeptonsName             (ModNames::gkMergedLeptonsName),
  fMetFromBranch                 (kTRUE),
  fPhotonsFromBranch             (kTRUE),
  fElectronsFromBranch           (kTRUE),
  fMuonsFromBranch               (kTRUE),

  // ----------------------------------------
  // Collections....
  fPhotons                       (0),
  fElectrons                     (0),
  fMuons                         (0),
  fMet                           (0),
  fMinNumPhotons                 (1),
  fMinNumLeptons                 (2),
  fMinPhotonEt                   (30),
  fMaxPhotonEta                  (2.4),
  fMinMetEt                      (30),
  // Counters....
  fNEventsSelected(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonAnalysisMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonAnalysisMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  // Load Branches
  ReqEventObject(fMetBranchName,    fMet,      fMetFromBranch);
  ReqEventObject(fPhotonsBranchName, fPhotons,  fPhotonsFromBranch); 
  ReqEventObject(fElectronsBranchName, fElectrons, fElectronsFromBranch);
  ReqEventObject(fMuonsBranchName, fMuons, fMuonsFromBranch);

  // Create your histograms here
  //*************************************************************************************************
  // Selection Histograms
  //*************************************************************************************************
  AddTH1(fMonoPhotonSelection,"hMonoPhotonSelection", ";Cuts;Number of Events",             5, 0, 5);
  
  // Create const char* labels with cut values
  std::ostringstream os1;
  std::string s1;
  os1<<"N Photons >= "<<fMinNumPhotons;
  s1 = os1.str();
  const char* labelCut1 = s1.c_str();
  std::ostringstream os2;
  std::string s2;
  os2<<"Photon Et >= "<<fMinPhotonEt;
  s2 = os2.str();
  const char* labelCut2 = s2.c_str();
  std::ostringstream os3;
  std::string s3;
  os3<<"Photon Eta <= "<<fMaxPhotonEta;
  s3 = os3.str();
  const char* labelCut3 = s3.c_str();
  std::ostringstream os4;
  std::string s4;
  os4<<"Met >= "<<fMinMetEt;
  s4 = os4.str();
  const char* labelCut4 = s4.c_str();
  std::ostringstream os5;
  std::string s5;
  os5<<"N Leptons >= "<<fMinNumLeptons;
  s5 = os5.str();
  const char* labelCut5 = s5.c_str();
  
  // Set selection histogram bin labels
  fMonoPhotonSelection->GetXaxis()->TAxis::SetBinLabel(1, "All Events");
  fMonoPhotonSelection->GetXaxis()->TAxis::SetBinLabel(2, labelCut1);
  fMonoPhotonSelection->GetXaxis()->TAxis::SetBinLabel(3, labelCut2);
  fMonoPhotonSelection->GetXaxis()->TAxis::SetBinLabel(4, labelCut3);
  fMonoPhotonSelection->GetXaxis()->TAxis::SetBinLabel(5, labelCut4);
  fMonoPhotonSelection->GetXaxis()->TAxis::SetBinLabel(6, labelCut5);

  //***********************************************************************************************
  // Histograms after preselection
  //***********************************************************************************************
  AddTH1(fPhotonEt           ,"hPhotonEt",";PhotonEt;Number of Events",400,0.,400.);
  AddTH1(fMetEt              ,"hMetEt",";MetEt;Number of Events",400,0.,400.);
  
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonAnalysisMod::Process()
{
  // Process entries of the tree.
  LoadEventObject(fPhotonsBranchName,   fPhotons);
  assert(fPhotons);
  LoadEventObject(fElectronsBranchName,   fElectrons);
  assert(fElectrons);
  LoadEventObject(fMuonsBranchName,    fMuons);
  assert(fMuons);
  LoadEventObject(fMetBranchName,   fMet);
  assert(fMet);
  ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(fLeptonsName);

  //*********************************************************************************************
  // Define Cuts
  //*********************************************************************************************
  const int nCuts = 5;
  bool passCut[nCuts] = {
	  false, 
	  false, 
	  false,
	  false,
	  false,
	  };
	  
  //***********************************************************************************************
  // Discard events with no identified photons
  //***********************************************************************************************
  if (fPhotons->GetEntries() >= fMinNumPhotons)  passCut[0] = true;

  //***********************************************************************************************
  // Discard events with no hard photons
  //***********************************************************************************************
  int    nHardPh = 0;
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
	  const Photon *ph = fPhotons->At(i);
	  if ( ph->Et() <= fMinPhotonEt ) continue; 
	  nHardPh ++;
  }
  if ( nHardPh > 0 ) passCut[1] = true;

  //***********************************************************************************************
  // Discard events with no hard photons inside eta range
  //***********************************************************************************************
  int    nHardEtaPh = 0;
  vector<int> theGoodPhotons;
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
	  const Photon *ph = fPhotons->At(i);
	  if ( ph->Et() <= fMinPhotonEt ) continue;
	  if ( TMath::Abs(ph->Eta()) >= fMaxPhotonEta ) continue;
	  nHardEtaPh ++;
    theGoodPhotons.push_back((int) i);
  }
  if ( nHardEtaPh > 0 ) passCut[2] = true;

  //***********************************************************************************************
  // Discard events with soft met
  //***********************************************************************************************
  const Met *stdMet = fMet->At(0);
  if ( stdMet->Pt() >= fMinMetEt )  passCut[3] = true;

  //***********************************************************************************************
  // Discard events with too few leptons
  //***********************************************************************************************
  if ((leptons->GetEntries() >= fMinNumLeptons) || (fElectrons->GetEntries() >= fMinNumLeptons) || (fMuons->GetEntries() >= fMinNumLeptons)) passCut[4] = true;
	  
  //*********************************************************************************************
  // Make Selection Histograms. Number of events passing each level of cut
  //*********************************************************************************************  
  //Cut Selection Histograms
  int zero = 0;
  fMonoPhotonSelection->Fill(zero,1);

  for (int k=0;k<nCuts;k++) {
    bool pass = true;
    bool passPreviousCut = true;
    for (int p=0;p<=k;p++) {
      pass = (pass && passCut[p]);
      if (p<k)
        passPreviousCut = (passPreviousCut&& passCut[p]);
    }
    if (pass)
      fMonoPhotonSelection->Fill(k+1,1);
  }

  
  bool passAllCuts = true;
  for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
  if(passAllCuts) {
	  fNEventsSelected++;
	  
	  // Make Preselection Histograms
    for ( int iGoodPh=0; iGoodPh<(int) theGoodPhotons.size(); iGoodPh++ )
      fPhotonEt->Fill(fPhotons->At(theGoodPhotons[iGoodPh])->Et()); 

	  fMetEt->Fill(fMet->At(0)->Et()); 
  }
  else 
	  this->SkipEvent(); // skip the event if does not passes the cuts
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonAnalysisMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
  cout << "selected events on MonoPhotonAnalysisMod: " << fNEventsSelected << endl;

} 
//--------------------------------------------------------------------------------------------------
void MonoPhotonAnalysisMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}
