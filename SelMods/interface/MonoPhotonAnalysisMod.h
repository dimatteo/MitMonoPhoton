//--------------------------------------------------------------------------------------------------
// $Id $
//
// MonoPhotonAnalysisMod
//
// A Module for Selecting gamma+MET events
// and produces some distributions.
//
//
// Authors: LDM, TRS
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
      void                SetInputPhotonsName(const char *n){ fPhotonsBranchName= n;        }
      void                SetInputMetName(const char *n){ fMetBranchName= n;        }
      void                SetPhotonsFromBranch(Bool_t b)    { fPhotonsFromBranch = b; }
      void                SetMetFromBranch(Bool_t b)    { fMetFromBranch = b; }

    protected:
      TString                  fPhotonsBranchName;	   //name of input photon branch
      TString                  fMetBranchName;           //name of input met branch
      Bool_t                   fPhotonsFromBranch;       //photons are loaded from a branch
      Bool_t                   fMetFromBranch;           //met is loaded from a branch
      
      TH1D                    *fMonoPhotonSelection;     //histogram for cut flow monitoring

      TH1D                    *fPhotonEt;                //histogram of photon transverse energy spectrum
      TH1D                    *fMetEt;                   //histogram of met spectrum

      const PhotonCol         *fPhotons;
      const PFMetCol          *fMet;

      void         SetMinNumPhotons(Int_t n)   { fMinNumPhotons = n; }
      void         SetMinPhotonEt(Double_t x)  { fMinPhotonEt = x;   }
      void         SetMaxPhotonEta(Double_t x) { fMaxPhotonEta = x;  }
      void         SetMinMetEt(Double_t x)     { fMinMetEt = x;      }

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      unsigned int fMinNumPhotons;
      Double_t     fMinPhotonEt;
      Double_t     fMaxPhotonEta;
      Double_t     fMinMetEt;

      Int_t                    fNEventsSelected;         //selected events

      ClassDef(MonoPhotonAnalysisMod,1) // TAM example analysis module
  };
}
#endif
