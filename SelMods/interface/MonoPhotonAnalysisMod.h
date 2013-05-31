//--------------------------------------------------------------------------------------------------
// $Id $
//
// MonoPhotonAnalysisMod
//
// A Module for Selecting gamma+MET events
// and produces some distributions.
//
//
// Authors: LDM
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOPHOTON_SELMODS_MONOPHOTONANALYSISMOD_H
#define MITMONOPHOTON_SELMODS_MONOPHOTONANALYSISMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"


class TH1D;
class TH2D;

namespace mithep 
{
  class MonoPhotonAnalysisMod : public BaseMod
  {
    public:
    MonoPhotonAnalysisMod(const char *name="MonoPhotonAnalysisMod", 
		 const char *title="Example analysis module with all branches");
      ~MonoPhotonAnalysisMod() {}

      // setting all the input Names
      void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
      void                SetInputMetName(const char *n){ fMetBranchName= n;        }

    protected:
      TString                  fPhotonBranchName;	     //name of input photon branch
      TString                  fMetBranchName;           //name of input met branch
      Int_t                    fNEventsSelected;         //selected events

      TH1D                    *fHWWSelection;            //histogram for cut flow monitoring

      TH1D                    *fPhotonEt;                //histogram of photon transverse energy spectrum
      TH1D                    *fMetEt;                   //histogram of met spectrum

      const PhotonCol              *fPhotons;
      const PFMetCol               *fMet;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(MonoPhotonAnalysisMod,1) // TAM example analysis module
  };
}
#endif
