//--------------------------------------------------------------------------------------------------
// $Id: MonoPhotonTreeWriter.h,v 1.5 2013/06/22 01:16:12 dimatteo Exp $
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
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"


class TNtuple;
class TRandom3;

namespace mithep 
{
  
  class MonoPhotonEvent
  {
    public:  
      // ------------ PHOTON STUFF -------------------
      static const Int_t kMaxPh = 5;
      Int_t   nPhotons;
      //kin
      Float_t a_photonE[kMaxPh];
      Float_t a_photonEt[kMaxPh];
      Float_t a_photonEta[kMaxPh];
      Float_t a_photonPhi[kMaxPh];
      //iso
      Float_t a_photonHCALisoDR03[kMaxPh];
      Float_t a_photonECALisoDR03[kMaxPh];
      Float_t a_photonHollowConeTKisoDR03[kMaxPh];
      Float_t a_photonHCALisoDR04[kMaxPh];
      Float_t a_photonECALisoDR04[kMaxPh];
      Float_t a_photonHollowConeTKisoDR04[kMaxPh];
      //shape
      Float_t a_photonCoviEtaiEta[kMaxPh];
      Float_t a_photonR9[kMaxPh];
      //misc
      Float_t a_photonSeedTime[kMaxPh];
      Float_t a_photonHadOverEm[kMaxPh];

      // ------------ VERTEX STUFF -------------------
      Float_t bsX;
      Float_t bsY;
      Float_t bsZ;
      Float_t bsSigmaZ; 
      Float_t vtxX;
      Float_t vtxY;      
      Float_t vtxZ;
      Int_t   nVtx;
      Float_t numPU;
      Float_t numPUminus;
      Float_t numPUplus;
      UInt_t  evt;
      UInt_t  run;
      UInt_t  lumi;
      UChar_t evtcat;

      // ------------ MET STUFF -------------------
      Float_t pfmet;
      Float_t pfmetphi;
      Float_t pfmetSig;
      Float_t pfSumEt;

      // ------------ JET STUFF -------------------
      static const Int_t kMaxJet = 5;
      Int_t   nJets;
      Float_t a_jetE[kMaxJet];
      Float_t a_jetPt[kMaxJet];
      Float_t a_jetEta[kMaxJet];
      Float_t a_jetPhi[kMaxJet];
      Float_t a_jetMass[kMaxJet];

      // ------------ ELECTRON STUFF -------------------
      static const Int_t kMaxEle = 2;
      Int_t   nElectrons;
      Float_t a_eleE[kMaxEle];
      Float_t a_elePt[kMaxEle];
      Float_t a_eleEta[kMaxEle];
      Float_t a_elePhi[kMaxEle];

      // ------------ MUON STUFF -------------------
      static const Int_t kMaxMu = 2;
      Int_t   nMuons;
      Float_t a_muE[kMaxMu];
      Float_t a_muPt[kMaxMu];
      Float_t a_muEta[kMaxMu];
      Float_t a_muPhi[kMaxMu];

      // ------------ Tracks STUFF -------------------
      static const Int_t kMaxTrack = 5;
      Int_t   nTracks;
      Float_t a_trkPt[kMaxTrack];
      Float_t a_trkEta[kMaxTrack];
      Float_t a_trkPhi[kMaxTrack];

  };
  
  
  class MonoPhotonTreeWriter : public BaseMod
  {
  public:
    MonoPhotonTreeWriter(const char *name ="MonoPhotonTreeWriter", 
		         const char *title="Selecting PhotonPairs");
    
    ~MonoPhotonTreeWriter();

    // setting all the input Names
    void                SetMetName(const char *n)         { fMetName= n;                 }
    void                SetPhotonsName(const char *n)     { fPhotonsName= n;              }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetElectronsName(const char *n)   { fElectronsName = n;          }
    void                SetElectronsFromBranch(bool b)    { fElectronsFromBranch = b;    }
    void                SetMuonsName(const char *n)       { fMuonsName = n;              }
    void                SetMuonsFromBranch(bool b)        { fMuonsFromBranch = b;        }
    void                SetJetsName(const char *n)        { fJetsName = n;               }
    void                SetJetsFromBranch(bool b)         { fJetsFromBranch = b;         }

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

    TString             fSuperClustersName;
    TString             fTracksName;
    TString             fPVName;
    TString             fPileUpDenName;    
    TString             fPileUpName;
    TString             fBeamspotName;
    
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
    const PileupInfoCol           *fPileUp;    
    const PileupEnergyDensityCol  *fPileUpDen;
    const SuperClusterCol         *fSuperClusters;   
    
    // --------------------------------
    // ntuple Tuple
    TString                        fTupleName;
    MonoPhotonEvent*               fMonoPhotonEvent;
    TTree*                         hMonoPhotonTuple;
    
    Int_t                          fNEventsSelected;

    ClassDef(MonoPhotonTreeWriter, 1) // Photon identification module
      };
}
#endif
