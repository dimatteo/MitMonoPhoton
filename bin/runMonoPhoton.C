// $Id: runMonoPhoton.C,v 1.1 2013/05/31 18:20:57 dimatteo Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

#include "MitMonoPhoton/SelMods/interface/MonoPhotonAnalysisMod.h"
#include "MitMonoPhoton/Mods/interface/MonoPhotonTreeWriter.h"

#endif

//--------------------------------------------------------------------------------------------------
void runMonoPhoton(const char *fileset    = "0000",
				   const char *skim       = "noskim",
				   const char *dataset    = "r12a-pho-j13-v1", 
				   const char *book       = "t2mit/filefi/029",
				   const char *catalogDir = "/home/cmsprod/catalog",
				   const char *outputName = "hgg",
				   int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;
  
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    sprintf(json, "%s", "~");
    //printf(" JSON file was not properly defined. EXIT!\n");
    //return;
  } 
  
  TString jsonFile = TString("/home/mingyang/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/mingyang/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) );
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
     sprintf(overlap,"%s", "-1.0");
    //printf(" OVERLAP file was not properly defined. EXIT!\n");
    //return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  MCProcessSelectionMod *mcselmod = new MCProcessSelectionMod;
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/mingyang/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information : [TBF]
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",160000,161176);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v*",161216,165633);

  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",165970,173198);
  
  hltModP->AddTrigger("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v*",165970,166967);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",167039,173198);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",173236,178380);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,189999);

  hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v*",190000,999999); 
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  // SepPUMod->SetVertexName("OutVtxCiC");
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kFALSE);
  
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC(!isData);

  GoodPVFilterMod *testPVFilterMod = new GoodPVFilterMod;
  testPVFilterMod->SetMinVertexNTracks(0);
  testPVFilterMod->SetMinNDof         (4.0);
  testPVFilterMod->SetMaxAbsZ         (24.0);
  testPVFilterMod->SetMaxRho          (2.0);
  testPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  testPVFilterMod->SetIsMC(!isData);
  
  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetPtMin(20);  
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);
  eleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml")))); 
  eleIdMod -> SetWhichVertex(-1);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.2); //h  
  eleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");
  eleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  eleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  eleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  eleIdMod -> SetNExpectedHitsInnerCut(1);   
  eleIdMod -> SetApplyConversionFilterType1(kTRUE);
  eleIdMod -> SetPVName(Names::gkPVBeamSpotBrn);   

  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  muonIdMod -> SetIDType("Tight");
  muonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.2);
  muonIdMod -> SetDZCut(0.5);
  muonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdMod -> SetPFIsoCut(0.2); //h
  muonIdMod -> SetOutputName("HggLeptonTagMuons");
  muonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  muonIdMod -> SetPFPileUpName("pfpileupcands");
  muonIdMod -> SetPVName(Names::gkPVBeamSpotBrn); 
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L3Absolute_AK5PF.txt")).Data())); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    
  
  
  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
  }
  
  if (TString(dataset).Contains("-qcd2em") || TString(dataset).Contains("-qcd-2em")) {
    excludedoubleprompt = kTRUE;
  }
  
  Bool_t is25 = kFALSE;
  if (TString(book).Contains("025")) is25 = kTRUE;
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(3);
  photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv3ph_52x.root")).Data()));
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetIsData(isData);
  
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetOutputVtxName("OutVtxCiC");        
  photcic->SetPhotonSelType("CiCPFSelection");
  photcic->SetVertexSelType("CiCMVA2012Selection");
  photcic->SetUseSingleLegConversions(kFALSE);
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  //photcic->SetMCSmearFactors2012HCP(0.0099,0.0098,0.0113,0.0094,0.0195,0.0198,0.0291,0.0275,0.0325,0.0356);
  photcic->SetMCSmearFactors2012HCP(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photcic->SetMCSmearFactors2012HCPMVA(0.01,0.011,0.0106,0.0106,0.0186,0.0196,0.0283,0.0267,0.0343,0.0345);
  photcic->UseSpecialSmearForDPMVA(true);
  photcic->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photcic->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photcic->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photcic->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photcic->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photcic->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photcic->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photcic->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photcic->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photcic->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photcic->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photcic->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photcic->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photcic->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photcic->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photcic->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);

  //-----------------------------------------------------------------------------------------------------------------
  
  //  photcic->SetDoMCR9Scaling(kTRUE);
  //  photcic->SetMCR9Scale(1.0035, 1.0035);
  photcic->SetDoShowerShapeScaling(kTRUE); 
  photcic->SetShowerShapeType("2012ShowerShape");
  //photcic->SetDoMCErrScaling(kTRUE);
  //photcic->SetMCErrScale(1.07, 1.045); 
  //photcic->SetMCErrScale(1, 1); //ming:scale(sigE/E)
  photcic->SetJetsName(jetCorr->GetOutputName());    
  //photcic->SetRescaledBeamspotWidth(5.0);
  photcic->SetIsData(isData);
  photcic->SetApplyLeptonTag(kTRUE);
  photcic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcic->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcic->Set2012HCP(kTRUE);
        
  PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
  phottreecic->SetPhotonsFromBranch(kFALSE);
  phottreecic->SetInputPhotonsName(photcic->GetOutputName());
  phottreecic->SetEnableJets(kTRUE);
  phottreecic->SetApplyJetId(kTRUE);
  phottreecic->SetPFJetsFromBranch(kFALSE);
  phottreecic->SetPFJetName(jetCorr->GetOutputName());
  phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
  phottreecic->SetIsData(isData);
  if (is25) phottreecic->SetEnablePFPhotons(kFALSE);
  phottreecic->SetApplyLeptonTag(kTRUE);
  phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  phottreecic->SetLeptonTagMuonsName("HggLeptonTagMuons");
  phottreecic->SetApplyVBFTag(kTRUE);
  phottreecic->SetApplyPFMetCorr(kTRUE);
  phottreecic->SetBeamspotWidth(5.0);
   
  //------------------------------------------------------------------------------------------------
  // select events with photon+MET
  //------------------------------------------------------------------------------------------------
  MonoPhotonAnalysisMod         *phplusmet = new MonoPhotonAnalysisMod("MonoPhotonSelector");
  phplusmet->SetMinNumPhotons(1);
  phplusmet->SetMinPhotonEt(30);
  phplusmet->SetMaxPhotonEta(2.4);
  phplusmet->SetMinMetEt(30);
  //phplusmet->SetInputPhotonsName(photreg->GetOutputName()); //phplusmet instead photcic?
  //photcic->SetInputMetName(photreg->GetOutputName());

  MonoPhotonTreeWriter *phplusmettree = new MonoPhotonTreeWriter("MonoPhotonTreeWriter");
  phplusmettree->SetPhotonsFromBranch(kTRUE);
  //phplusmettree->SetInputPhotonsName(photreg->GetOutputName());
  phplusmettree->SetEnableJets(kTRUE);
  //phplusmettree->SetApplyJetId(kTRUE);
  phplusmettree->SetPFJetsFromBranch(kFALSE);
  phplusmettree->SetPFJetName(jetCorr->GetOutputName());
  phplusmettree->SetExcludeDoublePrompt(excludedoubleprompt);
  phplusmettree->SetIsData(isData);
  if (is25) phplusmettree->SetEnablePFPhotons(kFALSE);
  phplusmettree->SetApplyLeptonTag(kTRUE);
  phplusmettree->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  phplusmettree->SetLeptonTagMuonsName("HggLeptonTagMuons");
  phplusmettree->SetApplyVBFTag(kTRUE);
  phplusmettree->SetApplyPFMetCorr(kTRUE);
  phplusmettree->SetBeamspotWidth(5.0);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
    
  if (!TString(dataset).Contains("meridiani")) {      
  //mcselmod         ->Add(hltModP);
  //hltModP         ->Add(goodPVFilterMod);
    mcselmod->Add(goodPVFilterMod);
  }
  else {
    mcselmod->Add(goodPVFilterMod);
  }
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  
  // simple object id modules
  jetCorr          ->Add(SepPUMod); 
  SepPUMod         ->Add(eleIdMod);
  eleIdMod         ->Add(muonIdMod);
   
  // Test gamma+met ana
  muonIdMod        ->Add(phplusmet);
  phplusmet        ->Add(phplusmettree);
  //phplusmet        ->Add(testPVFilterMod);

  // Working version of Hgg ana
  //muonIdMod          ->Add(photcic);
  //photcic         ->Add(phottreecic);

  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  TFile::SetReadaheadSize(128*1024*1024);
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  if (TString(dataset).Contains("meridiani")) {    
    ana->SetUseHLT(kFALSE);
  }
  
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset);
  //ana->AddDataset(d);
  ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/s12-qcd-2em3040-v7a/3E54E8F4-D2ED-E111-98CE-0026189438B1.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  ana->SetCacheSize(0);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
