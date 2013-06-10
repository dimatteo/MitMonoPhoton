//--------------------------------------------------------------------------------------------------
// $Id: MonoPhotonTreeWriter.h,v 1.1 2013/05/31 18:20:56 dimatteo Exp $
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
      static const Int_t kMaxPh = 10;
	  Int_t   nPhotons;
      Float_t a_photonE[kMaxPh];
      Float_t a_photonEt[kMaxPh];
      Float_t a_photonEta[kMaxPh];
      Float_t a_photonPhi[kMaxPh];
      Float_t a_photonHCALisoDR03[kMaxPh];
      Float_t a_photonECALisoDR03[kMaxPh];
      Float_t a_photonHollowConeTKisoDR03[kMaxPh];
      Float_t a_photonHCALisoDR04[kMaxPh];
      Float_t a_photonECALisoDR04[kMaxPh];
      Float_t a_photonHollowConeTKisoDR04[kMaxPh];
            
      // ------------ VERTEX STUFF -------------------
      Float_t bsX;
      Float_t bsY;
      Float_t bsZ;
      Float_t bsSigmaZ; 
      Float_t vtxX;
      Float_t vtxY;      
      Float_t vtxZ;
      Int_t   nVtx;
      Int_t   numPU;
      Int_t   numPUminus;
      Int_t   numPUplus;
      UInt_t  evt;
      UInt_t  run;
      UInt_t  lumi;
      UChar_t evtcat;
      UInt_t  nobj;

      // ------------ MET STUFF -------------------
      Float_t pfmet;
      Float_t pfmetphi;
      Float_t pfmetx;
      Float_t pfmety;

      // ------------ BTAG STUFF -------------------
      Float_t jet1pt;
      Float_t jet1eta;
      Float_t jet1phi;
      Float_t jet1mass;

  };
  
  
  class MonoPhotonTreeWriter : public BaseMod
  {
  public:
    MonoPhotonTreeWriter(const char *name ="MonoPhotonTreeWriter", 
		       const char *title="Selecting PhotonPairs");
    
    ~MonoPhotonTreeWriter();

    // setting all the input Names
    void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetTrackName(const char *n)       { fTrackBranchName = n;        }
    void                SetElectronName(const char *n)    { fElectronName = n;           }
    void                SetConversionName(const char *n)  { fConversionName = n;         }
    void                SetPUDensityName(const char *n)   { fPileUpDenName = n;          }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetMCParticle(const char *n)      { fMCParticleName = n;         }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }
    void                SetPFCandName(const char *n)      { fPFCandName = n;             }
    void                SetSuperClusterName(const char *n) { fSuperClusterName = n;      }
    void                SetPFJetName(const char *n)       { fPFJetName = n;              }
    void                SetGenJetName(const char *n)      { fGenJetName = n;             }
    void                SetuncorrPFJetName(const char *n) { funcorrPFJetName = n;        }
    void                SetPFNoPileUpName(const char *n)  { fPFNoPileUpName  = n;       } 
    void                SetPFPileUpName(const char *n)    { fPFPileUpName  = n;         }


    void                SetPFJetsFromBranch(Bool_t b)     { fPFJetsFromBranch = b;       }
    void                SetEnableJets(Bool_t b)           { fEnableJets = b;             }
    void                SetApplyJetId(Bool_t b)           { fApplyJetId = b;             }
    void                SetApplyLeptonTag(Bool_t b)       { fApplyLeptonTag = b;         }
    void                SetApplyVBFTag(Bool_t b)          { fApplyVBFTag = b;            }
    void                SetApplyBTag(Bool_t b)            { fApplyBTag = b;              }
    void                SetApplyPFMetCorr(Bool_t b)       { fApplyPFMetCorrections = b;  }
    void                SetPhFixDataFile(const char *n)   { fPhFixDataFile = n;          }




    // set basic Cut variables (FOR PRE-SELECTION)

    // is Data Or Not?
    void                SetIsData (Bool_t b)                 { fIsData = b; };
    

    void                SetApplyElectronVeto(Bool_t b)   { fApplyElectronVeto = b;     }          

    void                SetTupleName(const char* c)          { fTupleName = c; }
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name)    { fGoodElectronName = name; }
    void                SetWriteDiphotonTree(Bool_t b)       { fWriteDiphotonTree = b; }
    void                SetWriteSingleTree(Bool_t b)         { fWriteSingleTree = b; }
    void                SetLoopOnGoodElectrons(Bool_t b)     { fLoopOnGoodElectrons = b; }
    void                SetEnablePFPhotons(Bool_t b)         { fEnablePFPhotons = b; }
    void                SetExcludeSinglePrompt(Bool_t b)     { fExcludeSinglePrompt = b; }
    void                SetExcludeDoublePrompt(Bool_t b)     { fExcludeDoublePrompt = b; }

    void                SetLeptonTagElectronsName(TString name) { fLeptonTagElectronsName = name; }
    void                SetLeptonTagMuonsName    (TString name) { fLeptonTagMuonsName     = name; }
    void                SetFillClusterArrays(Bool_t b)          { fFillClusterArrays = b; }

    void                SetDo2012LepTag(Bool_t b) {                         fDo2012LepTag = b; }
    
    void                SetBeamspotWidth(Double_t x)            { fBeamspotWidth = x; }

      void                SetElectronMVAWeightsSubdet0Pt10To20(TString s)  
                          { fElectronMVAWeights_Subdet0Pt10To20  = s; }
      void                SetElectronMVAWeightsSubdet1Pt10To20(TString s)  
                          { fElectronMVAWeights_Subdet1Pt10To20  = s; }
      void                SetElectronMVAWeightsSubdet2Pt10To20(TString s)  
                          { fElectronMVAWeights_Subdet2Pt10To20  = s; }
      void                SetElectronMVAWeightsSubdet0Pt20ToInf(TString s) 
                          { fElectronMVAWeights_Subdet0Pt20ToInf = s; }
      void                SetElectronMVAWeightsSubdet1Pt20ToInf(TString s) 
                          { fElectronMVAWeights_Subdet1Pt20ToInf = s; }
      void                SetElectronMVAWeightsSubdet2Pt20ToInf(TString s) 
                          { fElectronMVAWeights_Subdet2Pt20ToInf = s; }

      void                SetDoSynching(bool b) {fDoSynching = b;}


    
  protected:
    void                Process();
    void                SlaveBegin();
    // Private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z, Float_t& mass);
    Float_t             GetEventCat    (PhotonTools::CiCBaseLineCats cat1,
					PhotonTools::CiCBaseLineCats cat2);

    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fPFPhotonName;
    TString             fElectronName;
    TString             fGoodElectronName;
    TString             fConversionName;
    TString             fPFConversionName;
    TString             fTrackBranchName;
    TString             fPileUpDenName;    
    TString             fPVName;
    TString             fBeamspotName;
    TString             fPFCandName;
    TString             fPFNoPileUpName;  //name of pfnpu collection
    TString             fPFPileUpName;    //name of pfpu collection

    TString             fMCParticleName;
    TString             fMCEventInfoName;
    
    TString             fPileUpName;
    TString             fSuperClusterName;
    TString             fPFMetName;
    TString             fPFJetName;


    TString             funcorrPFJetName;
    TString             fGenJetName;   //added to do pfmet correction 05/01/2012


    TString             fLeptonTagElectronsName;
    TString             fLeptonTagMuonsName;

    
    // is it Data or MC?
    Bool_t              fIsData;
    
    // in case there's some PV pre-selection
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t              fPFJetsFromBranch;

    Bool_t              fDoSynching;

    const PhotonCol               *fPhotons;
    const PhotonCol               *fPFPhotons;
    const ElectronCol             *fElectrons;
    const ElectronCol             *fGoodElectrons;    
    const DecayParticleCol        *fConversions;
    const DecayParticleCol        *fPFConversions;
    const TrackCol                *fTracks;
    const PileupEnergyDensityCol  *fPileUpDen;
    const VertexCol               *fPV;
    const BeamSpotCol             *fBeamspot;
    const PFCandidateCol          *fPFCands;
    const MCParticleCol           *fMCParticles;
    const MCEventInfo             *fMCEventInfo;
    const PileupInfoCol           *fPileUp;    
    const SuperClusterCol         *fSuperClusters;   
    const PFMetCol                *fPFMet;
    const JetCol                  *fPFJets;
    const GenJetCol               *fGenJets;
    const PFJetCol                *funcorrPFJets;
    
    const ElectronCol             *fLeptonTagElectrons;
    const MuonCol                 *fLeptonTagMuons;
    const PFCandidateCol          *fPFNoPileUpCands;  //!pfnpu collection
    const PFCandidateCol          *fPFPileUpCands;    //!pfpu collection

    // --------------------------------
    Bool_t                         fLoopOnGoodElectrons; //loop over good elecs instead of photons
    Bool_t                         fApplyElectronVeto;   //invert elec veto (for cic sel. only atm)
    Bool_t                         fWriteDiphotonTree;
    Bool_t                         fWriteSingleTree;

    Bool_t                         fEnablePFPhotons;
    
    Bool_t                         fExcludeSinglePrompt;
    Bool_t                         fExcludeDoublePrompt;
    
    Bool_t                         fEnableJets;
    Bool_t                         fApplyJetId;

    Bool_t                         fApplyLeptonTag;
    Bool_t                         fApplyVBFTag;
    Bool_t                         fApplyBTag;
    Bool_t                         fApplyPFMetCorrections;

    Bool_t                         fFillClusterArrays;
    Bool_t                         fFillVertexTree;
    
    Bool_t                         fDo2012LepTag;

    TString                        fPhFixDataFile;
    PhotonFix                      fPhfixph;
    PhotonFix                      fPhfixele;
    
    Double_t                       fBeamspotWidth;

    // --------------------------------
    // validation Tuple
    TString                        fTupleName;
    MonoPhotonEvent*               fMonoPhotonEvent;
    TTree*                         hMonoPhotonTuple;

    MVAMet                         fMVAMet;
    JetIDMVA                       fJetId;
    MVAVBF                         fMVAVBF;   
    
    VertexTools           fVtxTools;

    ElectronIDMVA                 *fElectronIDMVA;
    TString                   fElectronMVAWeights_Subdet0Pt10To20;
    TString                   fElectronMVAWeights_Subdet1Pt10To20;
    TString                   fElectronMVAWeights_Subdet2Pt10To20;
    TString                   fElectronMVAWeights_Subdet0Pt20ToInf;
    TString                   fElectronMVAWeights_Subdet1Pt20ToInf;
    TString                   fElectronMVAWeights_Subdet2Pt20ToInf;

    RhoUtilities::RhoType    fTheRhoType;

    ClassDef(MonoPhotonTreeWriter, 1) // Photon identification module
      };
}
#endif
