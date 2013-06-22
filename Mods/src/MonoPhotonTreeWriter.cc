#include "MitMonoPhoton/Mods/interface/MonoPhotonTreeWriter.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
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

  fSuperClustersName      ("PFSuperClusters"),
  fTracksName             (Names::gkTrackBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fBeamspotName           (Names::gkBeamSpotBrn),

  fIsData                 (false),
  fPhotonsFromBranch      (kTRUE),  
  fElectronsFromBranch    (kTRUE),  
  fMuonsFromBranch        (kTRUE),  
  fPVFromBranch           (kTRUE),
  fJetsFromBranch         (kTRUE),

  // ----------------------------------------
  // collections....
  fPhotons                (0),
  fElectrons              (0),
  fTracks                 (0),
  fPileUpDen              (0),
  fPV                     (0),
  fBeamspot               (0),
  fPileUp                 (0),
  fSuperClusters          (0),
  fJets                   (0),

  fTupleName              ("hMonoPhotonTree")

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

  //if(GetEventHeader()->EvtNum()==9008 || GetEventHeader()->EvtNum()==9008 || GetEventHeader()->EvtNum()==9010){
  //  printf("check MonoPhotonTreeWriter 0\n");
  // }
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fMetName,           fMet,            true);
  LoadEventObject(fPhotonsName,       fPhotons,        fPhotonsFromBranch); 
  LoadEventObject(fElectronsName,     fElectrons,      fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,          fMuonsFromBranch);
  LoadEventObject(fJetsName,          fJets,           fJetsFromBranch);

  LoadEventObject(fPVName,            fPV);    
  LoadEventObject(fBeamspotName,      fBeamspot);
  
  LoadEventObject(fSuperClustersName, fSuperClusters);
  LoadEventObject(fTracksName,        fTracks,        true);

  LoadEventObject(fPileUpDenName,     fPileUpDen,    true);
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
  }

  // ------------------------------------------------------------  
  // load event based information
  Int_t _numPU      = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus  = -1.;        // some sensible default values....
      
  if( !fIsData ) {
  LoadBranch(fPileUpName);
  } 
  if( !fIsData ) {
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing()==0) _numPU = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }

  fMonoPhotonEvent->nVtx = fPV->GetEntries();
  fMonoPhotonEvent->bsX = fBeamspot->At(0)->X();
  fMonoPhotonEvent->bsY = fBeamspot->At(0)->Y();
  fMonoPhotonEvent->bsZ = fBeamspot->At(0)->Z();
  fMonoPhotonEvent->bsSigmaZ = fBeamspot->At(0)->SigmaZ();
  fMonoPhotonEvent->vtxX = (fMonoPhotonEvent->nVtx>0) ? fPV->At(0)->X() : -99.;
  fMonoPhotonEvent->vtxY = (fMonoPhotonEvent->nVtx>0) ? fPV->At(0)->Y() : -99.;  
  fMonoPhotonEvent->vtxZ = (fMonoPhotonEvent->nVtx>0) ? fPV->At(0)->Z() : -99.;
  fMonoPhotonEvent->numPU = _numPU;
  fMonoPhotonEvent->numPUminus = _numPUminus;
  fMonoPhotonEvent->numPUplus = _numPUplus;
  fMonoPhotonEvent->evt = GetEventHeader()->EvtNum();
  fMonoPhotonEvent->run = GetEventHeader()->RunNum();
  fMonoPhotonEvent->lumi = GetEventHeader()->LumiSec();
  fMonoPhotonEvent->pfmet = fMet->At(0)->Pt();
  fMonoPhotonEvent->pfmetphi = fMet->At(0)->Phi();
  fMonoPhotonEvent->pfmetx = fMet->At(0)->Px();
  fMonoPhotonEvent->pfmety = fMet->At(0)->Py();
  
  //JETS  
  fMonoPhotonEvent->nJets = fJets->GetEntries();
  for ( int arrayIndex=0; arrayIndex<fMonoPhotonEvent->kMaxJet; arrayIndex++ ) {
	  if ( fJets->GetEntries() > 0 && arrayIndex < (int) fJets->GetEntries() ) {
		  const Jet *jet = fJets->At(arrayIndex);
		  //kin
		  fMonoPhotonEvent->a_jetE[arrayIndex] = jet->E();
		  fMonoPhotonEvent->a_jetPt[arrayIndex] = jet->Pt();
		  fMonoPhotonEvent->a_jetEta[arrayIndex] = jet->Eta();
		  fMonoPhotonEvent->a_jetPhi[arrayIndex] = jet->Phi();
		  fMonoPhotonEvent->a_jetMass[arrayIndex] = jet->Mass();
	  }
	  else {
		  //kin
		  fMonoPhotonEvent->a_jetE[arrayIndex] = -1;
		  fMonoPhotonEvent->a_jetPt[arrayIndex] = -1;
		  fMonoPhotonEvent->a_jetEta[arrayIndex] = -100;
		  fMonoPhotonEvent->a_jetPhi[arrayIndex] = -100;
		  fMonoPhotonEvent->a_jetMass[arrayIndex] = -1;
	  }
  }

        
  //PHOTONS  
  fMonoPhotonEvent->nPhotons = fPhotons->GetEntries();
  for ( int arrayIndex=0; arrayIndex<fMonoPhotonEvent->kMaxPh; arrayIndex++ ) {
	  if ( fPhotons->GetEntries() > 0 && arrayIndex < (int) fPhotons->GetEntries() ) {
		  const Photon *photon = fPhotons->At(arrayIndex);
		  //kin
		  fMonoPhotonEvent->a_photonE[arrayIndex] = photon->E();
		  fMonoPhotonEvent->a_photonEt[arrayIndex] = photon->Et();
		  fMonoPhotonEvent->a_photonEta[arrayIndex] = photon->Eta();
		  fMonoPhotonEvent->a_photonPhi[arrayIndex] = photon->Phi();
      //iso
		  fMonoPhotonEvent->a_photonHCALisoDR03[arrayIndex] = photon->HcalTowerSumEtDr03();
		  fMonoPhotonEvent->a_photonECALisoDR03[arrayIndex] = photon->EcalRecHitIsoDr03();
		  fMonoPhotonEvent->a_photonHollowConeTKisoDR03[arrayIndex] = photon->HollowConeTrkIsoDr03();
		  fMonoPhotonEvent->a_photonHCALisoDR04[arrayIndex] = photon->HcalTowerSumEtDr04();
		  fMonoPhotonEvent->a_photonECALisoDR04[arrayIndex] = photon->EcalRecHitIsoDr04();
		  fMonoPhotonEvent->a_photonHollowConeTKisoDR04[arrayIndex] = photon->HollowConeTrkIsoDr04();
      //shape
		  fMonoPhotonEvent->a_photonCoviEtaiEta[arrayIndex] = photon->CoviEtaiEta();
		  fMonoPhotonEvent->a_photonR9[arrayIndex] = photon->SCluster()->R9();
      //misc
		  fMonoPhotonEvent->a_photonSeedTime[arrayIndex] = photon->SCluster()->SeedTime();
		  fMonoPhotonEvent->a_photonHadOverEm[arrayIndex] = photon->HadOverEm();
	  }
	  else {
		  //kin
		  fMonoPhotonEvent->a_photonE[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonEt[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonEta[arrayIndex] = -100;
		  fMonoPhotonEvent->a_photonPhi[arrayIndex] = -100;
		  //iso
		  fMonoPhotonEvent->a_photonHCALisoDR03[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonECALisoDR03[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonHollowConeTKisoDR03[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonHCALisoDR04[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonECALisoDR04[arrayIndex] = -1;
		  fMonoPhotonEvent->a_photonHollowConeTKisoDR04[arrayIndex] = -1;
      //shape
		  fMonoPhotonEvent->a_photonCoviEtaiEta[arrayIndex] = -100;
		  fMonoPhotonEvent->a_photonR9[arrayIndex] = -1;
      //misc
		  fMonoPhotonEvent->a_photonSeedTime[arrayIndex] = -10000;
		  fMonoPhotonEvent->a_photonHadOverEm[arrayIndex] = -1;
	  }
  }

  //ELECTRONS  
  fMonoPhotonEvent->nElectrons = fElectrons->GetEntries();
  for ( int arrayIndex=0; arrayIndex<fMonoPhotonEvent->kMaxEle; arrayIndex++ ) {
	  if ( fElectrons->GetEntries() > 0 && arrayIndex < (int) fElectrons->GetEntries() ) {
		  const Electron *ele = fElectrons->At(arrayIndex);
		  //kin
		  fMonoPhotonEvent->a_eleE[arrayIndex] = ele->E();
		  fMonoPhotonEvent->a_elePt[arrayIndex] = ele->Pt();
		  fMonoPhotonEvent->a_eleEta[arrayIndex] = ele->Eta();
		  fMonoPhotonEvent->a_elePhi[arrayIndex] = ele->Phi();
	  }
	  else {
		  //kin
		  fMonoPhotonEvent->a_eleE[arrayIndex] = -1;
		  fMonoPhotonEvent->a_elePt[arrayIndex] = -1;
		  fMonoPhotonEvent->a_eleEta[arrayIndex] = -100;
		  fMonoPhotonEvent->a_elePhi[arrayIndex] = -100;
	  }
  }

  //ELECTRONS  
  fMonoPhotonEvent->nMuons = fMuons->GetEntries();
  for ( int arrayIndex=0; arrayIndex<fMonoPhotonEvent->kMaxMu; arrayIndex++ ) {
	  if ( fMuons->GetEntries() > 0 && arrayIndex < (int) fMuons->GetEntries() ) {
		  const Muon *mu = fMuons->At(arrayIndex);
		  //kin
		  fMonoPhotonEvent->a_muE[arrayIndex] = mu->E();
		  fMonoPhotonEvent->a_muPt[arrayIndex] = mu->Pt();
		  fMonoPhotonEvent->a_muEta[arrayIndex] = mu->Eta();
		  fMonoPhotonEvent->a_muPhi[arrayIndex] = mu->Phi();
	  }
	  else {
		  //kin
		  fMonoPhotonEvent->a_muE[arrayIndex] = -1;
		  fMonoPhotonEvent->a_muPt[arrayIndex] = -1;
		  fMonoPhotonEvent->a_muEta[arrayIndex] = -100;
		  fMonoPhotonEvent->a_muPhi[arrayIndex] = -100;
	  }
  }
  
  hMonoPhotonTuple -> Fill();
  
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
  }

  fMonoPhotonEvent = new MonoPhotonEvent;
  
  TFile *ftmp = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  hMonoPhotonTuple = new TTree(fTupleName.Data(),fTupleName.Data());
    
  //make flattish tree from classes so we don't have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::MonoPhotonEvent");
  TList  *elist  = eclass->GetListOfDataMembers();

  for (int i=0; i<elist->GetEntries(); ++i) {
    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i));//ming
    if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
    if (TString(tdm->GetName()).BeginsWith("kMax")) continue;
    TString typestring;
    if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
    else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
    else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
    else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
    else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
    else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
    else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
    else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
    else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
    else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
    else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
    else continue;
    //determine if the data member is an array
    bool dataMemberIsArray = false;
    if (TString(tdm->GetName()).BeginsWith("a_")) dataMemberIsArray = true;
    //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
    Char_t *addr = (Char_t*)fMonoPhotonEvent;//ming:?
    assert(sizeof(Char_t)==1);
    if ( dataMemberIsArray ) hMonoPhotonTuple->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s[10]/%s",tdm->GetName(),typestring.Data()));
    else hMonoPhotonTuple->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
  }
  
  AddOutput(hMonoPhotonTuple);
  
}
