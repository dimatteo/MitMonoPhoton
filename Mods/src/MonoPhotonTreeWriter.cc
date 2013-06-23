#include "MitMonoPhoton/Mods/interface/MonoPhotonTreeWriter.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/GenericParticle.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "TDataMember.h"
#include "TFile.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>

using namespace mithep;

ClassImp(mithep::MonoPhotonTreeWriter)

//--------------------------------------------------------------------------------------------------
MonoPhotonTreeWriter::MonoPhotonTreeWriter(const char *name, const char *title) : 
  // Base Module...
  BaseMod                 (name,title),

  fMetName                ("PFMet"),
  fPhotonsName            (Names::gkPhotonBrn),
  fElectronsName          (Names::gkElectronBrn),
  fMuonsName              (Names::gkMuonBrn),
  fJetsName               (Names::gkPFJetBrn),
  fLeptonsName            (ModNames::gkMergedLeptonsName),

  fSuperClustersName      ("PFSuperClusters"),
  fTracksName             (Names::gkTrackBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fBeamspotName           (Names::gkBeamSpotBrn),
  fMCEvInfoName           (Names::gkMCEvtInfoBrn),

  fIsData                 (false),
  fPhotonsFromBranch      (kTRUE),  
  fElectronsFromBranch    (kTRUE),  
  fMuonsFromBranch        (kTRUE),  
  fJetsFromBranch         (kTRUE),
  fPVFromBranch           (kTRUE),

  // ----------------------------------------
  fPhotons                (0),
  fElectrons              (0),
  fMuons                  (0),
  fJets                   (0),
  fTracks                 (0),
  fPV                     (0),
  fBeamspot               (0),
  fMCEventInfo            (0),
  fPileUp                 (0),
  fPileUpDen              (0),
  fSuperClusters          (0),

  fDecay(0),
  fOutputFile(0),
  fTupleName("hMonoPhotonTree"),

  fNEventsSelected(0)

{
  // Constructor
}

MonoPhotonTreeWriter::~MonoPhotonTreeWriter()
{
  // Destructor
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::Process()
{

  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fMetName,           fMet,           true);
  LoadEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch); 
  LoadEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  LoadEventObject(fJetsName,          fJets,          fJetsFromBranch);

  LoadEventObject(fPVName,            fPV,            fPVFromBranch);    
  LoadEventObject(fBeamspotName,      fBeamspot);
  
  LoadEventObject(fSuperClustersName, fSuperClusters);
  LoadEventObject(fTracksName,        fTracks,        true);

  LoadEventObject(fPileUpDenName,     fPileUpDen,     true);
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
  }
  ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);

  fNEventsSelected++;

  // ------------------------------------------------------------  
  // load event based information
      
  if( !fIsData ) {
  LoadBranch(fPileUpName);
  } 
  if( !fIsData ) {
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0) fMitGPTree.npu_	    = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1) fMitGPTree.npuPlusOne_  = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1) fMitGPTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
  }
  fMitGPTree.InitVariables();

  fMitGPTree.run_   = GetEventHeader()->RunNum();
  fMitGPTree.lumi_  = GetEventHeader()->LumiSec();
  fMitGPTree.event_ = GetEventHeader()->EvtNum();
  fMitGPTree.nvtx_  = fPV->GetEntries();
  fMitGPTree.scale1fb_ = 1000.0;
  
  if(fDecay == 0) fMitGPTree.dstype_ = MitGPTree::data;
  else            fMitGPTree.dstype_ = MitGPTree::other;

  fMitGPTree.met_    = fMet->At(0)->Pt();
  fMitGPTree.metPhi_ = fMet->At(0)->Phi();
  fMitGPTree.sumEt_  = fMet->At(0)->SumEt();
  fMitGPTree.metSig_ = fMet->At(0)->PFMetSig();

  // LEPTONS
  fMitGPTree.nlep_ = leptons->GetEntries();
  if (leptons->GetEntries() >= 1) {
    const Particle *lep = leptons->At(0);
    fMitGPTree.lep1_ = lep->Mom();
    if     (lep->ObjType() == kMuon    ) fMitGPTree.lid1_ = 13;
    else if(lep->ObjType() == kElectron) fMitGPTree.lid1_ = 11;
    else                                 assert(0);
    if(lep->Charge() < 0) fMitGPTree.lid1_ = -1 * fMitGPTree.lid1_;
  }
  if (leptons->GetEntries() >= 2) {
    Particle *lep = leptons->At(1);
    fMitGPTree.lep2_ = lep->Mom();
    if     (lep->ObjType() == kMuon    ) fMitGPTree.lid2_ = 13;
    else if(lep->ObjType() == kElectron) fMitGPTree.lid2_ = 11;
    else                                 assert(0);
    if(lep->Charge() < 0) fMitGPTree.lid2_ = -1 * fMitGPTree.lid2_;
  }
  if (leptons->GetEntries() >= 3) {
    Particle *lep = leptons->At(2);
    fMitGPTree.lep3_ = lep->Mom();
    if     (lep->ObjType() == kMuon    ) fMitGPTree.lid3_ = 13;
    else if(lep->ObjType() == kElectron) fMitGPTree.lid3_ = 11;
    else                                 assert(0);
    if(lep->Charge() < 0) fMitGPTree.lid3_ = -1 * fMitGPTree.lid3_;
  }

  //PHOTONS  
  fMitGPTree.nphotons_ = fPhotons->GetEntries();
  if(fPhotons->GetEntries() >= 1) {
    const Photon *photon = fPhotons->At(0);
    fMitGPTree.pho1_			  = photon->Mom();
    fMitGPTree.phoHCALisoDR03_a1_	  = photon->HcalTowerSumEtDr03();
    fMitGPTree.phoECALisoDR03_a1_	  = photon->EcalRecHitIsoDr03();
    fMitGPTree.phoHollowConeTKisoDR03_a1_ = photon->HollowConeTrkIsoDr03();
    fMitGPTree.phoHCALisoDR04_a1_	  = photon->HcalTowerSumEtDr04();
    fMitGPTree.phoECALisoDR04_a1_	  = photon->EcalRecHitIsoDr04();
    fMitGPTree.phoHollowConeTKisoDR04_a1_ = photon->HollowConeTrkIsoDr04();
    fMitGPTree.phoCoviEtaiEta_a1_	  = photon->CoviEtaiEta();
    fMitGPTree.phoR9_a1_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a1_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a1_ 	  = photon->HadOverEm();
  }
  if(fPhotons->GetEntries() >= 2) {
    const Photon *photon = fPhotons->At(1);
    fMitGPTree.pho2_			  = photon->Mom();
    fMitGPTree.phoHCALisoDR03_a2_	  = photon->HcalTowerSumEtDr03();
    fMitGPTree.phoECALisoDR03_a2_	  = photon->EcalRecHitIsoDr03();
    fMitGPTree.phoHollowConeTKisoDR03_a2_ = photon->HollowConeTrkIsoDr03();
    fMitGPTree.phoHCALisoDR04_a2_	  = photon->HcalTowerSumEtDr04();
    fMitGPTree.phoECALisoDR04_a2_	  = photon->EcalRecHitIsoDr04();
    fMitGPTree.phoHollowConeTKisoDR04_a2_ = photon->HollowConeTrkIsoDr04();
    fMitGPTree.phoCoviEtaiEta_a2_	  = photon->CoviEtaiEta();
    fMitGPTree.phoR9_a2_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a2_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a2_ 	  = photon->HadOverEm();
  }
  if(fPhotons->GetEntries() >= 3) {
    const Photon *photon = fPhotons->At(2);
    fMitGPTree.pho3_			  = photon->Mom();
    fMitGPTree.phoHCALisoDR03_a3_	  = photon->HcalTowerSumEtDr03();
    fMitGPTree.phoECALisoDR03_a3_	  = photon->EcalRecHitIsoDr03();
    fMitGPTree.phoHollowConeTKisoDR03_a3_ = photon->HollowConeTrkIsoDr03();
    fMitGPTree.phoHCALisoDR04_a3_	  = photon->HcalTowerSumEtDr04();
    fMitGPTree.phoECALisoDR04_a3_	  = photon->EcalRecHitIsoDr04();
    fMitGPTree.phoHollowConeTKisoDR04_a3_ = photon->HollowConeTrkIsoDr04();
    fMitGPTree.phoCoviEtaiEta_a3_	  = photon->CoviEtaiEta();
    fMitGPTree.phoR9_a3_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a3_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a3_ 	  = photon->HadOverEm();
  }
  if(fPhotons->GetEntries() >= 4) {
    const Photon *photon = fPhotons->At(3);
    fMitGPTree.pho4_			  = photon->Mom();
    fMitGPTree.phoHCALisoDR03_a4_	  = photon->HcalTowerSumEtDr03();
    fMitGPTree.phoECALisoDR03_a4_	  = photon->EcalRecHitIsoDr03();
    fMitGPTree.phoHollowConeTKisoDR03_a4_ = photon->HollowConeTrkIsoDr03();
    fMitGPTree.phoHCALisoDR04_a4_	  = photon->HcalTowerSumEtDr04();
    fMitGPTree.phoECALisoDR04_a4_	  = photon->EcalRecHitIsoDr04();
    fMitGPTree.phoHollowConeTKisoDR04_a4_ = photon->HollowConeTrkIsoDr04();
    fMitGPTree.phoCoviEtaiEta_a4_	  = photon->CoviEtaiEta();
    fMitGPTree.phoR9_a4_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a4_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a4_ 	  = photon->HadOverEm();
  }

  //JETS
  fMitGPTree.njets_ = fJets->GetEntries();
  if (fJets->GetEntries() >= 1) {
    const Jet *jet = fJets->At(0);
    fMitGPTree.jet1_     = jet->Mom();
    fMitGPTree.jet1Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }
  if (fJets->GetEntries() >= 2) {
    const Jet *jet = fJets->At(1);
    fMitGPTree.jet2_     = jet->Mom();
    fMitGPTree.jet2Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }
  if (fJets->GetEntries() >= 3) {
    const Jet *jet = fJets->At(2);
    fMitGPTree.jet3_     = jet->Mom();
    fMitGPTree.jet3Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }
  if (fJets->GetEntries() >= 4) {
    const Jet *jet = fJets->At(3);
    fMitGPTree.jet4_     = jet->Mom();
    fMitGPTree.jet4Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }
        
  //TRACKS
  fMitGPTree.ntracks_ = 0;
  for(unsigned int i = 0; i < fTracks->GetEntries(); i++) {
    const mithep::Track* pTrack = fTracks->At(i);
    if(pTrack->Pt() <= 15) continue;
    Bool_t isLepton = kFALSE;
    for (unsigned int arrayIndex=0; arrayIndex<leptons->GetEntries(); arrayIndex++ ) {
       const Particle *lep = leptons->At(arrayIndex);
      if(MathUtils::DeltaR(pTrack->Mom(), lep->Mom()) < 0.05) {
        isLepton = kTRUE;
        break;
      }
    }
    if(isLepton == kTRUE) continue;
    GenericParticle *p = new GenericParticle(pTrack->Px(), pTrack->Py(), pTrack->Pz(), 
                                             pTrack->P(), pTrack->Charge());
    if(fMitGPTree.ntracks_ == 0) fMitGPTree.track1_ = p->Mom();
    if(fMitGPTree.ntracks_ == 1) fMitGPTree.track2_ = p->Mom();
    if(fMitGPTree.ntracks_ == 2) fMitGPTree.track3_ = p->Mom();
    delete p;
    fMitGPTree.ntracks_++;
  }

  Double_t Q	     = 0.0;
  Int_t    id1       = 0;
  Double_t x1	     = 0.0;
  Double_t pdf1      = 0.0;
  Int_t    id2       = 0;
  Double_t x2	     = 0.0;
  Double_t pdf2      = 0.0;
  Int_t    processId = 0;
  if(fIsData == kFALSE){
     LoadBranch(fMCEvInfoName);
     Q         = fMCEventInfo->Scale();
     id1       = fMCEventInfo->Id1();
     x1        = fMCEventInfo->X1();
     pdf1      = fMCEventInfo->Pdf1();
     id2       = fMCEventInfo->Id2();
     x2        = fMCEventInfo->X2();
     pdf2      = fMCEventInfo->Pdf2();
     processId = fMCEventInfo->ProcessId();
  }
  fMitGPTree.Q_         = Q;
  fMitGPTree.id1_       = id1;
  fMitGPTree.x1_        = x1;
  fMitGPTree.pdf1_      = pdf1;
  fMitGPTree.id2_       = id2;
  fMitGPTree.x2_        = x2;
  fMitGPTree.pdf2_      = pdf2;
  fMitGPTree.processId_ = processId;

  fMitGPTree.tree_->Fill();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.
  ReqEventObject(fMetName,           fMet,            true);
  ReqEventObject(fPhotonsName,       fPhotons,        fPhotonsFromBranch); 
  ReqEventObject(fElectronsName,     fElectrons,      fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,          fMuonsFromBranch);
  ReqEventObject(fJetsName,          fJets,           fJetsFromBranch);

  ReqEventObject(fSuperClustersName,  fSuperClusters, true);
  ReqEventObject(fTracksName,         fTracks,        true);

  ReqEventObject(fPVName,             fPV,            fPVFromBranch);    
  ReqEventObject(fBeamspotName,       fBeamspot,      true);

  ReqEventObject(fPileUpDenName,   fPileUpDen,    true);
  if (!fIsData) {
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCEvInfoName,       fMCEventInfo);
  }

  //***********************************************************************************************
  //Create Smurf Ntuple Tree  
  //***********************************************************************************************
  fOutputFile = new TFile(fTupleName.Data(), "RECREATE");
  fMitGPTree.CreateTree(-1);
  
}
//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::SlaveTerminate()
{
  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  cout << "Processed events on MonoPhotonTreeWriter: " << fNEventsSelected << endl;
}
