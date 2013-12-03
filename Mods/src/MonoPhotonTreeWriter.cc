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
  fMetCorName             ("PFMetT0T1Shift"),
  fPhotonsName            (Names::gkPhotonBrn),
  fElectronsName          (Names::gkElectronBrn),
  fMuonsName              (Names::gkMuonBrn),
  fJetsName               (Names::gkPFJetBrn),
  fCosmicsName            ("CosmicMuons"),
  fLeptonsName            (ModNames::gkMergedLeptonsName),

  fGenISRPhotonsName      (ModNames::gkMCISRPhotonsName),
  fGenRadPhotonsName      (ModNames::gkMCRadPhotonsName),
  fGenPhotonsName         (ModNames::gkMCPhotonsName),
  fGenAllLeptonsName      (ModNames::gkMCAllLeptonsName),
  fGenJetsName            ("AKT5GenJets"),

  fSuperClustersName      ("PFSuperClusters"),
  fTracksName             (Names::gkTrackBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fBeamspotName           (Names::gkBeamSpotBrn),
  fMCEvInfoName           (Names::gkMCEvtInfoBrn),
  fAllElectronsName       ("Electrons"),
  fConversionsName        ("MergedConversions"),
  fPfCandidatesName       ("PFCandidates"),
  fHltObjsName            (Names::gkHltObjBrn),
  fEvtSelDataName         ("EvtSelData"),

  fIsData                 (false),
  fPhotonsFromBranch      (kTRUE),  
  fElectronsFromBranch    (kTRUE),  
  fMuonsFromBranch        (kTRUE),  
  fJetsFromBranch         (kTRUE),
  fCosmicsFromBranch      (kTRUE),
  fPVFromBranch           (kTRUE),

  // ----------------------------------------
  fPhotons                (0),
  fElectrons              (0),
  fMuons                  (0),
  fJets                   (0),
  fCosmics                (0),
  fTracks                 (0),
  fPV                     (0),
  fBeamspot               (0),
  fMCEventInfo            (0),
  fPileUp                 (0),
  fPileUpDen              (0),
  fSuperClusters          (0),
  fAllElectrons           (0),
  fConversions            (0),
  fPfCandidates           (0),
  fHltObjs                (0),
  fEvtSelData             (0),

  fDecay(0),
  fOutputFile(0),
  fTupleName("hMonoPhotonTree"),
  fFillNtupleType(0),
  fIsCicPhotonId(kTRUE),

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
  LoadEventObject(fMetCorName,        fMetCor,        false);
  LoadEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch); 
  LoadEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  LoadEventObject(fJetsName,          fJets,          fJetsFromBranch);
  LoadEventObject(fCosmicsName,       fCosmics,       fCosmicsFromBranch);
  LoadEventObject(fPVName,            fPV,            fPVFromBranch);    
  LoadEventObject(fBeamspotName,      fBeamspot);
  
  LoadEventObject(fSuperClustersName, fSuperClusters);
  LoadEventObject(fTracksName,        fTracks,        true);
  LoadEventObject(fAllElectronsName,  fAllElectrons);
  LoadEventObject(fConversionsName,   fConversions);
  LoadEventObject(fPfCandidatesName,  fPfCandidates);

  LoadEventObject(fPileUpDenName,     fPileUpDen,     true);
  LoadEventObject(fEvtSelDataName,    fEvtSelData,    true);
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    LoadEventObject(fGenISRPhotonsName ,    fGenISRPhotons, false);
    LoadEventObject(fGenRadPhotonsName ,    fGenRadPhotons, false);
    LoadEventObject(fGenPhotonsName    ,    fGenPhotons, false);
    LoadEventObject(fGenAllLeptonsName ,    fGenAllLeptons, false);
    LoadEventObject(fGenJetsName       ,    fGenJets, true);
  }
  ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  const TriggerObjectCol *fHltObjs = GetHLTObjects(fHltObjsName);
  if (!fHltObjs){
    printf("HLTEvtSelMod::TriggerObjectCol not found\n");
    return;
  }
  
  fNEventsSelected++;

  // ------------------------------------------------------------  
  // load event based information
  fMitGPTree.InitVariables();
      
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

  fMitGPTree.cuts_ = 0;  

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
  fMitGPTree.metCor_    = fMetCor->At(0)->Pt();
  fMitGPTree.metCorPhi_ = fMetCor->At(0)->Phi();

  fMitGPTree.metFilterWord_ = fEvtSelData->metFiltersWord();

  // LEPTONS
  fMitGPTree.nlep_ = leptons->GetEntries();
  // auxiliary, just to dump muon info
  MuonOArr *muons = GetObjThisEvt<MuonOArr>("HggLeptonTagMuons");
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
  float rho; //needed for isolation variables
  if (fPileUpDen->At(0)->RhoKt6PFJets()>0.) rho = fPileUpDen->At(0)->RhoKt6PFJets();
  else rho = fPileUpDen->At(0)->Rho();

  fMitGPTree.nphotons_ = fPhotons->GetEntries();
  if(fPhotons->GetEntries() >= 1) {
    const Photon *photon = fPhotons->At(0);
    //Get the relevant quantities for this photon
    int matchType = -1.;
    float matchPt = -1.;
    if (!fIsData)
      GetPhoMCMatch(photon->Eta(), photon->Phi(), photon->Pt(), matchType, matchPt);
    float combIso1,combIso2,combIso3;
    unsigned int wVtxInd = 0;
    if (fIsCicPhotonId) {
      wVtxInd = 0;
      float trackIsoSel03   = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.0, fPfCandidates);
      float trackIsoWorst04 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.4, 0.00, fPfCandidates, &wVtxInd, fPV);
      float ecalIso3        = IsolationTools::PFGammaIsolation(photon, 0.3, 0.0, fPfCandidates);
      float ecalIso4        = IsolationTools::PFGammaIsolation(photon, 0.4, 0.0, fPfCandidates);
      combIso1 = (ecalIso3+trackIsoSel03   + 2.5 - 0.09*rho)*50./photon->Et();
      combIso2 = (ecalIso4+trackIsoWorst04 + 2.5 - 0.23*rho)*50./(photon->MomVtx(fPV->At(wVtxInd)->Position()).Pt());
      combIso3 = (trackIsoSel03)*50./photon->Et();
    }
    else {
      wVtxInd = 0;
      float CHIso03 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.0, fPfCandidates);
      float NHIso03 = IsolationTools::PFNeutralHadronIsolation(photon, 0.3, 0.0, fPfCandidates);
      float PHIso03 = IsolationTools::PFGammaIsolation(photon, 0.3, 0.0, fPfCandidates);
    
      combIso1 = TMath::Max(CHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),0),(float) 0.);
      combIso2 = TMath::Max(NHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),1),(float) 0.);
      combIso3 = TMath::Max(PHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),2),(float) 0.);
    }
    wVtxInd = 0;
    float phoWorstIso03 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.00, fPfCandidates, &wVtxInd, fPV);
    Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(photon, fAllElectrons, fConversions, fBeamspot->At(0));  
    //Fill the photons branches
    fMitGPTree.pho1_			  = photon->Mom();
    fMitGPTree.phoMatchType_a1_	  = matchType;
    fMitGPTree.phoMatchPt_a1_	    = matchPt;
    fMitGPTree.phoCombIso1_a1_	  = combIso1;
    fMitGPTree.phoCombIso2_a1_	  = combIso2;
    fMitGPTree.phoCombIso3_a1_ = combIso3;
    fMitGPTree.phoWorstIso_a1_ = phoWorstIso03;
    fMitGPTree.phoPassEleVeto_a1_	  = PassEleVetoRaw;
    fMitGPTree.phoHasPixelSeed_a1_ = photon->HasPixelSeed();
    fMitGPTree.phoCoviEtaiEta_a1_	  = photon->CoviEtaiEta();
    fMitGPTree.phoCoviPhiiPhi_a1_	  = TMath::Sqrt(TMath::Abs(photon->SCluster()->Seed()->CoviPhiiPhi()));
    fMitGPTree.phoR9_a1_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a1_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a1_ 	  = photon->HadOverEm();
    fMitGPTree.phoLeadTimeSpan_a1_  	  = photon->SCluster()->LeadTimeSpan();
    fMitGPTree.phoRoundness_a1_  	  = photon->SCluster()->Roundness();
    fMitGPTree.phoAngle_a1_  	  = photon->SCluster()->Angle();
    fMitGPTree.phoMipChi2_a1_ = photon->MipChi2();
    fMitGPTree.phoMipTotEnergy_a1_ = photon->MipTotEnergy();
    fMitGPTree.phoMipSlope_a1_ = photon->MipSlope();
    fMitGPTree.phoMipIntercept_a1_ = photon->MipIntercept();
    fMitGPTree.phoMipNhitCone_a1_ = photon->MipNhitCone();
    fMitGPTree.phoMipIsHalo_a1_ = photon->MipIsHalo();
    // get the HLT matched object information
    for (UInt_t i=0; i<fHltObjs->GetEntries(); i++) {
      const TriggerObject *trigobj = fHltObjs->At(i);
      if (trigobj->TriggerType()==TriggerObject::TriggerPhoton && MathUtils::DeltaR(photon->SCluster(),trigobj)<0.3)
        fMitGPTree.phoIsTrigger_a1_ = true;
    }
  }
  if(fPhotons->GetEntries() >= 2) {
    const Photon *photon = fPhotons->At(1);
    //Get the relevant quantities for this photon
    int matchType = -1.;
    float matchPt = -1.;
    if (!fIsData)
      GetPhoMCMatch(photon->Eta(), photon->Phi(), photon->Pt(), matchType, matchPt);
    float combIso1,combIso2,combIso3;
    unsigned int wVtxInd = 0;
    if (fIsCicPhotonId) {
      wVtxInd = 0;
      float trackIsoSel03   = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.0, fPfCandidates);
      float trackIsoWorst04 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.4, 0.00, fPfCandidates, &wVtxInd, fPV);
      float ecalIso3        = IsolationTools::PFGammaIsolation(photon, 0.3, 0.0, fPfCandidates);
      float ecalIso4        = IsolationTools::PFGammaIsolation(photon, 0.4, 0.0, fPfCandidates);
      combIso1 = (ecalIso3+trackIsoSel03   + 2.5 - 0.09*rho)*50./photon->Et();
      combIso2 = (ecalIso4+trackIsoWorst04 + 2.5 - 0.23*rho)*50./(photon->MomVtx(fPV->At(wVtxInd)->Position()).Pt());
      combIso3 = (trackIsoSel03)*50./photon->Et();
    }
    else {
      wVtxInd = 0;
      float CHIso03 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.0, fPfCandidates);
      float NHIso03 = IsolationTools::PFNeutralHadronIsolation(photon, 0.3, 0.0, fPfCandidates);
      float PHIso03 = IsolationTools::PFGammaIsolation(photon, 0.3, 0.0, fPfCandidates);
    
      combIso1 = TMath::Max(CHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),0),(float) 0.);
      combIso2 = TMath::Max(NHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),1),(float) 0.);
      combIso3 = TMath::Max(PHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),2),(float) 0.);
    }
    wVtxInd = 0;
    float phoWorstIso03 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.00, fPfCandidates, &wVtxInd, fPV);
    Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(photon, fAllElectrons, fConversions, fBeamspot->At(0));  
    //Fill the photons branches
    fMitGPTree.pho2_			  = photon->Mom();
    fMitGPTree.phoMatchType_a2_	  = matchType;
    fMitGPTree.phoMatchPt_a2_	    = matchPt;
    fMitGPTree.phoCombIso1_a2_	  = combIso1;
    fMitGPTree.phoCombIso2_a2_	  = combIso2;
    fMitGPTree.phoCombIso3_a2_ = combIso3;
    fMitGPTree.phoWorstIso_a2_ = phoWorstIso03;
    fMitGPTree.phoPassEleVeto_a2_	  = PassEleVetoRaw;
    fMitGPTree.phoHasPixelSeed_a2_ = photon->HasPixelSeed();
    fMitGPTree.phoCoviEtaiEta_a2_	  = photon->CoviEtaiEta();
    fMitGPTree.phoCoviPhiiPhi_a2_	  = TMath::Sqrt(TMath::Abs(photon->SCluster()->Seed()->CoviPhiiPhi()));
    fMitGPTree.phoR9_a2_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a2_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a2_ 	  = photon->HadOverEm();
    fMitGPTree.phoLeadTimeSpan_a2_  	  = photon->SCluster()->LeadTimeSpan();
    fMitGPTree.phoRoundness_a2_  	  = photon->SCluster()->Roundness();
    fMitGPTree.phoAngle_a2_  	  = photon->SCluster()->Angle();
    fMitGPTree.phoMipChi2_a2_ = photon->MipChi2();
    fMitGPTree.phoMipTotEnergy_a2_ = photon->MipTotEnergy();
    fMitGPTree.phoMipSlope_a2_ = photon->MipSlope();
    fMitGPTree.phoMipIntercept_a2_ = photon->MipIntercept();
    fMitGPTree.phoMipNhitCone_a2_ = photon->MipNhitCone();
    fMitGPTree.phoMipIsHalo_a2_ = photon->MipIsHalo();
    // get the HLT matched object information
    for (UInt_t i=0; i<fHltObjs->GetEntries(); i++) {
      const TriggerObject *trigobj = fHltObjs->At(i);
      if (trigobj->TriggerType()==TriggerObject::TriggerPhoton && MathUtils::DeltaR(photon->SCluster(),trigobj)<0.3)
        fMitGPTree.phoIsTrigger_a2_ = true;
    }
  }
  if(fPhotons->GetEntries() >= 3) {
    const Photon *photon = fPhotons->At(2);
    //Get the relevant quantities for this photon
    int matchType = -1.;
    float matchPt = -1.;
    if (!fIsData)
      GetPhoMCMatch(photon->Eta(), photon->Phi(), photon->Pt(), matchType, matchPt);
    float combIso1,combIso2,combIso3;
    unsigned int wVtxInd = 0;
    if (fIsCicPhotonId) {
      wVtxInd = 0;
      float trackIsoSel03   = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.0, fPfCandidates);
      float trackIsoWorst04 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.4, 0.00, fPfCandidates, &wVtxInd, fPV);
      float ecalIso3        = IsolationTools::PFGammaIsolation(photon, 0.3, 0.0, fPfCandidates);
      float ecalIso4        = IsolationTools::PFGammaIsolation(photon, 0.4, 0.0, fPfCandidates);
      combIso1 = (ecalIso3+trackIsoSel03   + 2.5 - 0.09*rho)*50./photon->Et();
      combIso2 = (ecalIso4+trackIsoWorst04 + 2.5 - 0.23*rho)*50./(photon->MomVtx(fPV->At(wVtxInd)->Position()).Pt());
      combIso3 = (trackIsoSel03)*50./photon->Et();
    }
    else {
      wVtxInd = 0;
      float CHIso03 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.0, fPfCandidates);
      float NHIso03 = IsolationTools::PFNeutralHadronIsolation(photon, 0.3, 0.0, fPfCandidates);
      float PHIso03 = IsolationTools::PFGammaIsolation(photon, 0.3, 0.0, fPfCandidates);
    
      combIso1 = TMath::Max(CHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),0),(float) 0.);
      combIso2 = TMath::Max(NHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),1),(float) 0.);
      combIso3 = TMath::Max(PHIso03 - rho*GetEA(photon->SCluster()->AbsEta(),2),(float) 0.);
    }
    wVtxInd = 0;
    float phoWorstIso03 = IsolationTools::PFChargedIsolation(photon, fPV->At(0), 0.3, 0.00, fPfCandidates, &wVtxInd, fPV);
    Bool_t PassEleVetoRaw = PhotonTools::PassElectronVetoConvRecovery(photon, fAllElectrons, fConversions, fBeamspot->At(0));  
    //Fill the photons branches
    fMitGPTree.pho3_              = photon->Mom();
    fMitGPTree.phoMatchType_a3_	  = matchType;
    fMitGPTree.phoMatchPt_a3_	    = matchPt;
    fMitGPTree.phoCombIso1_a3_	  = combIso1;
    fMitGPTree.phoCombIso2_a3_	  = combIso2;
    fMitGPTree.phoCombIso3_a3_ = combIso3;
    fMitGPTree.phoWorstIso_a3_ = phoWorstIso03;
    fMitGPTree.phoPassEleVeto_a3_	  = PassEleVetoRaw;
    fMitGPTree.phoHasPixelSeed_a3_ = photon->HasPixelSeed();
    fMitGPTree.phoCoviEtaiEta_a3_	  = photon->CoviEtaiEta();
    fMitGPTree.phoCoviPhiiPhi_a3_	  = TMath::Sqrt(TMath::Abs(photon->SCluster()->Seed()->CoviPhiiPhi()));
    fMitGPTree.phoR9_a3_		  = photon->SCluster()->R9();
    fMitGPTree.phoSeedTime_a3_  	  = photon->SCluster()->SeedTime();
    fMitGPTree.phoHadOverEm_a3_ 	  = photon->HadOverEm();
    fMitGPTree.phoLeadTimeSpan_a3_  	  = photon->SCluster()->LeadTimeSpan();
    fMitGPTree.phoRoundness_a3_  	  = photon->SCluster()->Roundness();
    fMitGPTree.phoAngle_a3_  	  = photon->SCluster()->Angle();
    fMitGPTree.phoMipChi2_a3_ = photon->MipChi2();
    fMitGPTree.phoMipTotEnergy_a3_ = photon->MipTotEnergy();
    fMitGPTree.phoMipSlope_a3_ = photon->MipSlope();
    fMitGPTree.phoMipIntercept_a3_ = photon->MipIntercept();
    fMitGPTree.phoMipNhitCone_a3_ = photon->MipNhitCone();
    fMitGPTree.phoMipIsHalo_a3_ = photon->MipIsHalo();
    // get the HLT matched object information
    for (UInt_t i=0; i<fHltObjs->GetEntries(); i++) {
      const TriggerObject *trigobj = fHltObjs->At(i);
      if (trigobj->TriggerType()==TriggerObject::TriggerPhoton && MathUtils::DeltaR(photon->SCluster(),trigobj)<0.3)
        fMitGPTree.phoIsTrigger_a3_ = true;
    }
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
       
  //COSMIC MUONS
  fMitGPTree.ncosmics_ = 0;
  for(unsigned int i = 0; i < fCosmics->GetEntries(); i++) {
   
    const mithep::Muon* cMuon = fCosmics->At(i);    
    //save only cosmics with an associated valid standalone track
    if(!cMuon->HasStandaloneTrk() || !cMuon->IsStandaloneMuon()) continue;

    if(fMitGPTree.ncosmics_ == 0) {
      fMitGPTree.cosmic1_ = cMuon->Mom();
    }
    if(fMitGPTree.ncosmics_ == 1) {
      fMitGPTree.cosmic2_ = cMuon->Mom();
    }
    if(fMitGPTree.ncosmics_ == 2) {
      fMitGPTree.cosmic3_ = cMuon->Mom();
    }
    fMitGPTree.ncosmics_++;
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
  fMitGPTree.rho_       = rho;

  fMitGPTree.tree_->Fill();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.
  ReqEventObject(fMetName,           fMet,            true);
  ReqEventObject(fMetCorName,        fMetCor,         false);
  ReqEventObject(fPhotonsName,       fPhotons,        fPhotonsFromBranch); 
  ReqEventObject(fElectronsName,     fElectrons,      fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,          fMuonsFromBranch);
  ReqEventObject(fJetsName,          fJets,           fJetsFromBranch);
  ReqEventObject(fCosmicsName,       fCosmics,        fCosmicsFromBranch);

  ReqEventObject(fSuperClustersName,  fSuperClusters, true);
  ReqEventObject(fTracksName,         fTracks,        true);

  ReqEventObject(fPVName,             fPV,            fPVFromBranch);    
  ReqEventObject(fBeamspotName,       fBeamspot,      true);
  ReqEventObject(fAllElectronsName,   fAllElectrons,  true);
  ReqEventObject(fConversionsName,    fConversions,   true);
  ReqEventObject(fPfCandidatesName,   fPfCandidates,  true);

  ReqEventObject(fPileUpDenName,   fPileUpDen,    true);
  ReqEventObject(fEvtSelDataName,  fEvtSelData,   true);  
  if (!fIsData) {
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCEvInfoName,       fMCEventInfo);
    ReqEventObject(fGenISRPhotonsName,  fGenISRPhotons, false);
    ReqEventObject(fGenRadPhotonsName,  fGenRadPhotons, false);
    ReqEventObject(fGenPhotonsName   ,  fGenPhotons, false);
    ReqEventObject(fGenAllLeptonsName,  fGenAllLeptons, false);
    ReqEventObject(fGenJetsName      ,  fGenJets, true);
  }

  //***********************************************************************************************
  //Create Smurf Ntuple Tree  
  //***********************************************************************************************
  fMitGPTree.CreateTree(fFillNtupleType);
  AddOutput(fMitGPTree.tree_);
}
//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::SlaveTerminate()
{
  //fOutputFile->cd();
  //fOutputFile->Write();
  //fOutputFile->Close();
  cout << "Processed events on MonoPhotonTreeWriter: " << fNEventsSelected << endl;
}

//--------------------------------------------------------------------------------------------------
float MonoPhotonTreeWriter::GetCorrDeltaPhi(float phi1, float phi2)
{
  float corrDeltaPhi = TMath::Abs(phi1 - phi2);
  if (corrDeltaPhi > TMath::Pi())
    corrDeltaPhi = TMath::TwoPi() - corrDeltaPhi;     
  return corrDeltaPhi;
}

//--------------------------------------------------------------------------------------------------
float MonoPhotonTreeWriter::GetEA(float absEta, int isoType)
{
  // prepare the values of the EA (cat 0/1/2 is CH/NH/Photons)
  float effective_area[] = { 
    0.012, 	0.030, 	0.148,    
    0.010, 	0.057, 	0.130,    
    0.014, 	0.039, 	0.112,      
    0.012, 	0.015, 	0.216,      
    0.016, 	0.024, 	0.262,  
    0.020, 	0.039, 	0.260,
    0.012, 	0.072, 	0.266};

  int _etaCat = 0;    
  if (absEta >= 1 && absEta < 1.479)  _etaCat = 1;
  if (absEta >= 1.479 && absEta < 2.0)_etaCat = 2;
  if (absEta >= 2.0 && absEta < 2.2)  _etaCat = 3;
  if (absEta >= 2.2 && absEta < 2.3)  _etaCat = 4;
  if (absEta >= 2.3 && absEta < 2.4)  _etaCat = 5;
  if (absEta >= 2.4)                  _etaCat = 6;

  return effective_area[isoType+_etaCat*3];
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::GetPhoMCMatch(float eta, float phi, float pt, int &matchType, float &matchPt)
{
  // matched particle type
  matchType = -1; // default: unmatched

  // the DR for matching, require also pt matching at 50% for e/g and 100% for jets
  float theMinDRSq = 1000.;
  float genPt = -1.;

  // loop on ISR photons: first choice
  if (fGenISRPhotons->GetEntries() >= 1)
    for (UInt_t i=0; i < fGenISRPhotons->GetEntries(); ++i) {
      if (fGenISRPhotons->At(i)->Pt() < 30.) continue;
      float genEta = fGenISRPhotons->At(i)->Eta();
      float genPhi = fGenISRPhotons->At(i)->Phi();
      genPt        = fGenISRPhotons->At(i)->Pt();
      float thisDRSq =  (eta-genEta)*(eta-genEta)
                     +  GetCorrDeltaPhi(phi,genPhi)*GetCorrDeltaPhi(phi,genPhi);
      if (thisDRSq < theMinDRSq && fabs((pt-genPt)/genPt) < 1.5) theMinDRSq = thisDRSq;
    }
  // return an ISR match if match cone small enough
  if (theMinDRSq < 0.1*0.1) {
    matchType = 0;
    matchPt = genPt;
    return;
  }
    
  // loop on FSR photons: second choice
  theMinDRSq = 1000.;
  if (fGenRadPhotons->GetEntries() >= 1) 
    for (UInt_t i=0; i < fGenRadPhotons->GetEntries(); ++i) {
      if (fGenRadPhotons->At(i)->Pt() < 30.) continue;
      float genEta = fGenRadPhotons->At(i)->Eta();
      float genPhi = fGenRadPhotons->At(i)->Phi();
      genPt        = fGenRadPhotons->At(i)->Pt();
      float thisDRSq =  (eta-genEta)*(eta-genEta)
                     +  GetCorrDeltaPhi(phi,genPhi)*GetCorrDeltaPhi(phi,genPhi);
      if (thisDRSq < theMinDRSq && fabs((pt-genPt)/genPt) < 1.5) theMinDRSq = thisDRSq;
    }
  // return an FSR match if match cone small enough
  if (theMinDRSq < 0.1*0.1) {
    matchType = 1;
    matchPt = genPt;
    return;
  }

  // loop on photons: third choice
  theMinDRSq = 1000.;
  if (fGenPhotons->GetEntries() >= 1)
    for (UInt_t i=0; i < fGenPhotons->GetEntries(); ++i) {
      if (fGenPhotons->At(i)->Pt() < 30.) continue;
      float genEta = fGenPhotons->At(i)->Eta();
      float genPhi = fGenPhotons->At(i)->Phi();
      genPt        = fGenPhotons->At(i)->Pt();
      float thisDRSq =  (eta-genEta)*(eta-genEta)
                     +  GetCorrDeltaPhi(phi,genPhi)*GetCorrDeltaPhi(phi,genPhi);
      if (thisDRSq < theMinDRSq && fabs((pt-genPt)/genPt) < 1.5) theMinDRSq = thisDRSq;
    }
  // return an photon match if match cone small enough
  if (theMinDRSq < 0.1*0.1) {
    matchType = 2;
    matchPt = genPt;
    return;
  }
    
  // loop on electrons: fourth choice
  theMinDRSq = 1000.;
  if (fGenAllLeptons->GetEntries() >= 1)  
    for (UInt_t i=0; i < fGenAllLeptons->GetEntries(); ++i) {
      if (fGenAllLeptons->At(i)->Pt() < 30.) continue;
      if (fabs(fGenAllLeptons->At(i)->PdgId()) != 11) continue;
      float genEta = fGenAllLeptons->At(i)->Eta();
      float genPhi = fGenAllLeptons->At(i)->Phi();
      genPt        = fGenAllLeptons->At(i)->Pt();
      float thisDRSq =  (eta-genEta)*(eta-genEta)
                     +  GetCorrDeltaPhi(phi,genPhi)*GetCorrDeltaPhi(phi,genPhi);
      if (thisDRSq < theMinDRSq && fabs((pt-genPt)/genPt) < 1.5) theMinDRSq = thisDRSq;
    }
  // return an electron match if match cone small enough
  if (theMinDRSq < 0.1*0.1) {
    matchType = 3;
    matchPt = genPt;
    return;
  }

  // loop on genJets: fifth choice
  theMinDRSq = 1000.;
  if (fGenJets->GetEntries() >= 1)  
    for (UInt_t i=0; i < fGenJets->GetEntries(); ++i) {
      if (fGenJets->At(i)->Pt() < 30.) continue;
      float genEta = fGenJets->At(i)->Eta();
      float genPhi = fGenJets->At(i)->Phi();
      genPt        = fGenJets->At(i)->Pt();
      float thisDRSq =  (eta-genEta)*(eta-genEta)
                     +  GetCorrDeltaPhi(phi,genPhi)*GetCorrDeltaPhi(phi,genPhi);
      if (thisDRSq < theMinDRSq && fabs((pt-genPt)/genPt) < 2) theMinDRSq = thisDRSq;
    }
  // return a jet match if match cone small enough
  if (theMinDRSq < 0.5*0.5) {
    matchType = 4;
    matchPt = genPt;
    return;
  }
  
  return;
                  
}
