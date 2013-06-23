//--------------------------------------------------------------------------------------------------
// $Id: MonoPhotonTreeWriter.h,v 1.6 2013/06/22 20:58:23 ceballos Exp $
//
// MonoPhotonTreeWriter
//
// Authors: L. Di Matteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOPHOTON_MODS_MONOPHOTONTREEWRITER_H
#define MITMONOPHOTON_MODS_MONOPHOTONTREEWRITER_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

#include "MitMonoPhoton/Core/MitGPTree.h"

class TNtuple;
class TRandom3;

namespace mithep 
{

  class MonoPhotonTreeWriter : public BaseMod
  {
  public:
    MonoPhotonTreeWriter(const char *name ="MonoPhotonTreeWriter", 
		         const char *title="Selecting PhotonPairs");
    
    ~MonoPhotonTreeWriter();

    // setting all the input Names
    void                SetMetName(const char *n)         { fMetName= n;                 }
    void                SetPhotonsName(const char *n)     { fPhotonsName= n;             }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetElectronsName(const char *n)   { fElectronsName = n;          }
    void                SetElectronsFromBranch(bool b)    { fElectronsFromBranch = b;    }
    void                SetMuonsName(const char *n)       { fMuonsName = n;              }
    void                SetMuonsFromBranch(bool b)        { fMuonsFromBranch = b;        }
    void                SetJetsName(const char *n)        { fJetsName = n;               }
    void                SetJetsFromBranch(bool b)         { fJetsFromBranch = b;         }
    void                SetLeptonsName(const char *n)     { fLeptonsName = n;            }

    void                SetSuperClustersName(const char *n){ fSuperClustersName = n;     }
    void                SetTracksName(const char *n)      { fTracksName = n;             }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }

    // is Data Or Not?
    void                SetIsData (Bool_t b)              { fIsData = b; };

    void                SetTupleName(const char* c)          { fTupleName = c; }

  protected:
    void                Process();
    void                SlaveBegin();
    void                SlaveTerminate();
    // Private auxiliary methods...
    // Names of the input Collections
    TString             fMetName;
    TString             fPhotonsName;
    TString             fElectronsName;
    TString             fMuonsName;
    TString             fJetsName;
    TString             fLeptonsName;

    TString             fSuperClustersName;
    TString             fTracksName;
    TString             fPVName;
    TString             fPileUpDenName;    
    TString             fPileUpName;
    TString             fBeamspotName;
    TString             fMCEvInfoName;

    // is it Data or MC?
    Bool_t              fIsData;
    
    // there are not some PV pre-selection?
    Bool_t              fPhotonsFromBranch;
    Bool_t              fElectronsFromBranch;
    Bool_t              fMuonsFromBranch;
    Bool_t              fJetsFromBranch;
    Bool_t              fPVFromBranch;

    const PFMetCol                *fMet;
    const PhotonCol               *fPhotons;
    const ElectronCol             *fElectrons;
    const MuonCol                 *fMuons;
    const JetCol                  *fJets;

    const TrackCol                *fTracks;
    const VertexCol               *fPV;
    const BeamSpotCol             *fBeamspot;
    const MCEventInfo             *fMCEventInfo;
    const PileupInfoCol           *fPileUp;    
    const PileupEnergyDensityCol  *fPileUpDen;
    const SuperClusterCol         *fSuperClusters;   
    
    // --------------------------------
    TFile	                   *fOutputFile;
    TString	                   fTupleName;
    MitGPTree                      fMitGPTree;

    Int_t                          fNEventsSelected;

    ClassDef(MonoPhotonTreeWriter, 1) // Photon identification module
      };
}
#endif
