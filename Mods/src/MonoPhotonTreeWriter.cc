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
  // define all the Branches to load
  fPhotonBranchName       (Names::gkPhotonBrn),
  fPFPhotonName           ("PFPhotons"),
  fElectronName           (Names::gkElectronBrn),
  fGoodElectronName       (Names::gkElectronBrn),  
  fConversionName         (Names::gkMvfConversionBrn),  
  fPFConversionName              ("PFPhotonConversions"),  
  fTrackBranchName        (Names::gkTrackBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fBeamspotName           (Names::gkBeamSpotBrn),
  fPFCandName             (Names::gkPFCandidatesBrn),
  fPFNoPileUpName         ("PFNoPileUp"),
  fPFPileUpName           ("PFPileUp"),
  fMCParticleName         (Names::gkMCPartBrn),
  fMCEventInfoName        (Names::gkMCEvtInfoBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fSuperClusterName       ("PFSuperClusters"),
  fPFMetName              ("PFMet"),
  fPFJetName              (Names::gkPFJetBrn),
  funcorrPFJetName        ("AKt5PFJets"),
  fGenJetName             ("AKT5GenJets"),
  fLeptonTagElectronsName ("HggLeptonTagElectrons"),
  fLeptonTagMuonsName     ("HggLeptonTagMuons"),

  fIsData                 (false),
  fPhotonsFromBranch      (kTRUE),  
  fPVFromBranch           (kTRUE),
  fGoodElectronsFromBranch(kTRUE),
  fPFJetsFromBranch       (kTRUE),
  // ----------------------------------------
  // flag for synchronization, adds vertex variables
  // should be on for synching trees
  fDoSynching             (kFALSE),

  // ----------------------------------------
  // collections....
  fPhotons                (0),
  fPFPhotons              (0),
  fElectrons              (0),
  fConversions            (0),
  fPFConversions          (0),
  fTracks                 (0),
  fPileUpDen              (0),
  fPV                     (0),
  fBeamspot               (0),
  fPFCands                (0),
  fMCParticles            (0),
  fMCEventInfo            (0),
  fPileUp                 (0),
  fSuperClusters          (0),
  fPFJets                 (0),
  fGenJets                (0),
  funcorrPFJets           (0),

  fLeptonTagElectrons     (0),
  fLeptonTagMuons         (0),
  fPFNoPileUpCands        (0),
  fPFPileUpCands          (0),

  fLoopOnGoodElectrons    (kFALSE),
  fApplyElectronVeto      (kTRUE),  
  fWriteSingleTree        (kTRUE),
  fEnablePFPhotons        (kTRUE),
  fExcludeSinglePrompt    (kFALSE),
  fExcludeDoublePrompt    (kFALSE),
  fEnableJets             (kFALSE),
  fApplyJetId             (kFALSE),
  fApplyLeptonTag         (kFALSE),
  fApplyVBFTag            (kFALSE),
  fApplyBTag              (kFALSE),
  fApplyPFMetCorrections  (kFALSE),
  fFillClusterArrays      (kFALSE),
  fFillVertexTree         (kFALSE),
  fDo2012LepTag           (kFALSE),
  fPhFixDataFile          (gSystem->Getenv("CMSSW_BASE") +
		           TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat")),
  fBeamspotWidth          (5.8),

  fElectronIDMVA(0),
  fElectronMVAWeights_Subdet0Pt10To20(""),
  fElectronMVAWeights_Subdet1Pt10To20(""),
  fElectronMVAWeights_Subdet2Pt10To20(""),
  fElectronMVAWeights_Subdet0Pt20ToInf(""),
  fElectronMVAWeights_Subdet1Pt20ToInf(""),
  fElectronMVAWeights_Subdet2Pt20ToInf(""),
  fTheRhoType(RhoUtilities::DEFAULT),

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
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);

  // lepton tag collections
  if( fApplyLeptonTag ) {
    LoadEventObject(fLeptonTagElectronsName, fLeptonTagElectrons);
    LoadEventObject(fLeptonTagMuonsName,     fLeptonTagMuons);
  }
  
  if (fEnablePFPhotons) LoadEventObject(fPFPhotonName,   fPFPhotons);
  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fConversionName,     fConversions);
  if ( fDoSynching ) LoadEventObject(fPFConversionName,     fPFConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);    
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fSuperClusterName,   fSuperClusters);
  LoadEventObject(fPFMetName,          fPFMet);  
//   LoadEventObject(fPFNoPileUpName,     fPFNoPileUpCands);
//   LoadEventObject(fPFPileUpName,     fPFPileUpCands);

  if (fEnableJets){
    LoadEventObject(fPFJetName,        fPFJets);  
    //LoadEventObject(funcorrPFJetName,  funcorrPFJets);
    LoadBranch(funcorrPFJetName);
    //   if(!fIsData) LoadEventObject(fGenJetName,        fGenJets);
  }

  // ------------------------------------------------------------  
  // load event based information
  Int_t _numPU      = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus  = -1.;        // some sensible default values....
      
  if( !fIsData ) {
    LoadBranch(fMCParticleName);
    LoadBranch(fMCEventInfoName);
    LoadBranch(fPileUpName);
    if (fEnableJets) LoadEventObject(fGenJetName,        fGenJets);
  }  else fGenJets = NULL;
  
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
  fMonoPhotonEvent->nobj = fPhotons->GetEntries();
  fMonoPhotonEvent->pfmet = fPFMet->At(0)->Pt();
  fMonoPhotonEvent->pfmetphi = fPFMet->At(0)->Phi();
  fMonoPhotonEvent->pfmetx = fPFMet->At(0)->Px();
  fMonoPhotonEvent->pfmety = fPFMet->At(0)->Py();
  
  //JETS  
  fMonoPhotonEvent->nJets = fPFJets->GetEntries();
  for ( int arrayIndex=0; arrayIndex<fMonoPhotonEvent->kMaxJet; arrayIndex++ ) {
	  if ( fPFJets->GetEntries() > 0 && arrayIndex < (int) fPFJets->GetEntries() ) {
		  const Jet *jet = fPFJets->At(arrayIndex);
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
  
  hMonoPhotonTuple -> Fill();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MonoPhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  if( fApplyLeptonTag ) {
    ReqEventObject(fLeptonTagElectronsName,    fLeptonTagElectrons,    false);  
    ReqEventObject(fLeptonTagMuonsName,        fLeptonTagMuons,        false);  
  }

//   ReqEventObject(fPFNoPileUpName,     fPFNoPileUpCands,    false);
//   ReqEventObject(fPFPileUpName,     fPFPileUpCands,    false);

  ReqEventObject(fPhotonBranchName,fPhotons,      fPhotonsFromBranch);
  if (fEnablePFPhotons) ReqEventObject(fPFPhotonName,fPFPhotons,      true);
  ReqEventObject(fTrackBranchName, fTracks,       true);
  ReqEventObject(fElectronName,    fElectrons,    true);  
  ReqEventObject(fGoodElectronName,fGoodElectrons,fGoodElectronsFromBranch);  
  ReqEventObject(fPileUpDenName,   fPileUpDen,    true);
  ReqEventObject(fPVName,          fPV,           fPVFromBranch);
  ReqEventObject(fConversionName,  fConversions,  true);
  if ( fDoSynching ) ReqEventObject(fPFConversionName,     fPFConversions,  true);
  ReqEventObject(fBeamspotName,    fBeamspot,     true);
  ReqEventObject(fPFCandName,      fPFCands,      true);
  ReqEventObject(fSuperClusterName,fSuperClusters,true);
  ReqEventObject(fPFMetName,       fPFMet,        true);
  if (fEnableJets){
    ReqEventObject(fPFJetName,       fPFJets,       fPFJetsFromBranch);
    ReqBranch(funcorrPFJetName, funcorrPFJets);
    //   if (!fIsData) ReqEventObject(fGenJetName, fGenJets, true);
  }
  if (!fIsData) {
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCParticleName,     fMCParticles);
    ReqBranch(fMCEventInfoName,    fMCEventInfo);
    if (fEnableJets) ReqEventObject(fGenJetName, fGenJets, true);
  }
  if (fIsData) {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") +
      TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  }
  else {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") +
      TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");
  }

  fPhfixph.initialise("4_2",std::string(fPhFixDataFile));
  fPhfixele.initialise("4_2e",std::string(fPhFixDataFile));
  
//   fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
//                       TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2_42.root")))
//                       );  
  
  fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1cov_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2cov_52.root")))
                      );                      
                 
  fJetId.Initialize(JetIDMVA::kMedium,
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                          JetIDMVA::kCut,
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")))

                          );

  fMVAVBF.InitializeMVA();
                      

  if( fDoSynching ) {
    fVtxTools.InitP(2);
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kIDEGamma2012NonTrigV1,
			       fTheRhoType);
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
