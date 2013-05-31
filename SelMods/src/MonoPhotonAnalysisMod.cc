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

using namespace mithep;
ClassImp(mithep::MonoPhotonAnalysisMod)

//--------------------------------------------------------------------------------------------------
MonoPhotonAnalysisMod::MonoPhotonAnalysisMod(const char *name, const char *title) : 
  BaseMod(name,title),
  // define all the Branches to load
  fPhotonBranchName              (Names::gkPhotonBrn),
  fMetBranchName              ("PFMet"),
  // ----------------------------------------
  // collections....
  fPhotons                       (0),
  fMet                       (0),
  // counters....
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
  ReqEventObject(fPhotonBranchName,   fPhotons,      true);
  ReqEventObject(fMetBranchName,   fMet,      true);

  //Create your histograms here

  //*************************************************************************************************
  // Selection Histograms
  //*************************************************************************************************
  AddTH1(fHWWSelection,"hHWWSelection", ";Cut Number;Number of Events",             17, -1.5, 15.5);

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
  LoadEventObject(fPhotonBranchName,   fPhotons);
  assert(fPhotons);
  LoadEventObject(fMetBranchName,   fMet);
  assert(fMet);

  //*********************************************************************************************
  //Define Cuts
  //*********************************************************************************************
  const int nCuts = 3;
  bool passCut[nCuts] = {
	  false, 
	  false, 
	  false,
	  };
	  
  //***********************************************************************************************
  //Discard events with no identified photons
  //***********************************************************************************************
  if (fPhotons->GetEntries() < 1)
	  return;
  passCut[0] = true;

  //***********************************************************************************************
  //Discard events with soft photons
  //***********************************************************************************************
  double minPhEt = 30;
  int    nHardPh = 0;
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
	  const Photon *ph = fPhotons->At(i);
	  if ( ph->Et() < 30 ) continue; 
	  nHardPh ++;
  }
  if ( nHardPh == 0 )
      return;
  passCut[1] = true;

  //***********************************************************************************************
  //Discard events with soft met
  //***********************************************************************************************
  const Met *stdMet = fMet->At(0);
  if ( stdMet->Pt() < 30 )
      return;
  passCut[2] = true;
	  
  //*********************************************************************************************
  //Make Selection Histograms. Number of events passing each level of cut
  //*********************************************************************************************  
  bool passAllCuts = true;
  for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
  if(passAllCuts) fNEventsSelected++;
  else this->SkipEvent(); //skip the event if does not passes the cuts

  //Cut Selection Histograms
  fHWWSelection->Fill(-1,1);

  for (int k=0;k<nCuts;k++) {
    bool pass = true;
    bool passPreviousCut = true;
    for (int p=0;p<=k;p++) {
      pass = (pass && passCut[p]);
      if (p<k)
        passPreviousCut = (passPreviousCut&& passCut[p]);
    }
    if (pass)
      fHWWSelection->Fill(k,1);
  }
  
  //*****************************************************************************************
  //Make Preselection Histograms  
  //*****************************************************************************************
  fPhotonEt->Fill(fPhotons->At(0)->Et()); 
  fMetEt->Fill(fMet->At(0)->Et()); 
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
