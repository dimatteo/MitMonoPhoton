#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TPad.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TBox.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TTree.h>
#include "MitMonoPhoton/Utils/interface/TreeReducer.h"
#include "MitMonoPhoton/Core/MitGPTree.h"
#include "MitMonoPhoton/Core/MitGPTreeReduced.h"

#include "TLorentzVector.h"

ClassImp(mithep::TreeReducer)

using namespace std;
using namespace mithep;

const TH1D *TreeReducer::sPUWeights = 0;
const TH1D *TreeReducer::sPUWeightsUp = 0;
const TH1D *TreeReducer::sPUWeightsDown = 0;

//--------------------------------------------------------------------------------------------------
TreeReducer::TreeReducer(const Sample *mySample) :
  fSample      (mySample),
  fPUTarget    (0),
  fPUTargetUp  (0),
  fPUTargetDown(0),
  fOutFile     (0),
  fVerbose     (0),
  fInputBaseDir(0),
  fLumi(0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
TreeReducer::~TreeReducer()
{
  // Destructor
}

//--------------------------------------------------------------------------------------------------
void TreeReducer::MakeTree()
{
  // say what we are doing
  if ( fVerbose ) printf("\n ==== Reducing tree for -- %s ====\n\n",fSample->Name()->Data());

  // some useful definitions
  TString *dirFwk  = new TString("AnaFwkMod");
  TString *allEvts = new TString("hDAllEvents");
  TString  slash  ("/");

  // make sure the sample file exists
  TString inFileName = fInputBaseDir+slash+*fSample->File();

  TFile *fif = new TFile(inFileName.Data());
  if (fif->IsOpen() == kFALSE) {
    printf(" WARNING -- sample  %s  does not have a histogram file. Continue without!\n",
           fSample->Name()->Data());
    return;
  }

  // determine the sample type
  bool isMC = false, isSignal = true, isEleFake = false;
  int treeType = 0;
  if (fInputSelection == "Lepton" || fInputSelection == "DiLepton" ) {
    treeType = 1;
    isSignal = false;
  }
  if (fInputSelection == "EleFake") {
    isEleFake = true;
  }
  float thisXsec   = *fSample->Xsec();
  float thisScale  = *fSample->Scale();
  if (thisXsec > 0)         isMC = true;
  if (thisXsec == -2) {treeType = 2;} //QCD
  if (thisXsec == -3) {isEleFake = true;} //EleFake
  if (thisXsec == -4 || fIsBHStudy) {treeType = 3;} //BeamHalo

  // determine if overlap needs to be removed (only for Zjets and Wjets)
  bool isOverlapSample = false;
  if (*fSample->Name() == "s12-wjets-ptw100-v7a"
    ||*fSample->Name() == "s12-zjets-ptz100-v7a" ) isOverlapSample = true;

  // Get the event scaling and the PU weights: MC only
  float baseEvtWeight = 1.;
  if ( isMC ) {
    TDirectory *dirTmp = (TDirectory*) gROOT->FindObject(dirFwk->Data());
    if (dirTmp)
      fif->cd(dirFwk->Data());
    TH1D *hAllEvts = (TH1D*) gROOT->FindObject(allEvts->Data());
    if (! hAllEvts) {
      printf(" WARNING -- sample  %s  does not have a framework file. Next sample!\n",
             fSample->Name()->Data());
      return;
    }
    float nGenEvts   = hAllEvts->GetEntries();
    baseEvtWeight *= fLumi*thisXsec/nGenEvts;
    
    //set pileup weights
    if (fPUTarget) {
      if (sPUWeights) {
        delete sPUWeights;
        delete sPUWeightsUp;
        delete sPUWeightsDown;
    }
      TH1D *pusource = (TH1D*)dirTmp->Get("hNPUTrue")->Clone();
      pusource->Rebin(10);
      pusource->Scale(1.0/pusource->GetSumOfWeights());
          
      sPUWeights = new TH1D( (*fPUTarget) / (*pusource) );
      sPUWeightsUp = new TH1D( (*fPUTargetUp) / (*pusource) );
      sPUWeightsDown = new TH1D( (*fPUTargetDown) / (*pusource) );
    }

  }

  // set up the output tree
  fOutFile   -> cd();
  MitGPTreeReduced outtree;
  TString outTreeName = (*fSample->Name());
  if (treeType == 2 || treeType == 3 || thisXsec == -3) outTreeName = (*fSample->Legend());
  outtree.CreateTree(outTreeName.Data());

  // set up the input tree
  MitGPTree intree;
  intree.LoadTree(inFileName,treeType);
  intree.InitTree(0);

  // set up the resolution functions
  SetResoFunctions(!isMC);
  
  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {

    outtree.InitVariables();
    intree.tree_-> GetEntry(iEntry);
    
    //cut on lepton if not signal region
    if (!isSignal && intree.nlep_ == 0) continue;
    
    //apply the selection and define the good photon         
    int theGoodPhoton = -1;
    if ( fIsEGSelection && !EventIsSelectedEG(intree, treeType, theGoodPhoton, isEleFake) ) continue;
    if ( !fIsEGSelection && !EventIsSelected(intree, treeType, theGoodPhoton, isEleFake) ) continue;

    //get this event weights, not photon Et dependant
    float thisPUWeight = 1.;
    float thisPUWeightUp = 1.;
    float thisPUWeightDown = 1.;
    if ( isMC ) thisPUWeight = PUWeight(intree.npu_);
    if ( isMC ) thisPUWeightUp = PUWeightUp(intree.npu_);
    if ( isMC ) thisPUWeightDown = PUWeightDown(intree.npu_);
    
    //get the relevant variables
    outtree.isData_ = !isMC;
    
    outtree.evt_weight_ = baseEvtWeight;
    outtree.pu_weight_  = thisPUWeight;
    outtree.pu_weightup_  = thisPUWeightUp;
    outtree.pu_weightdo_  = thisPUWeightDown;
  
    outtree.nvtx_  = intree.nvtx_;
    outtree.metCor_= intree.metCor_;
    outtree.metCorPhi_= intree.metCorPhi_;
    outtree.metSig_= intree.metSig_;
    outtree.metFilterWord_ = intree.metFilterWord_;
  
    outtree.nphotons_ = intree.nphotons_;
    outtree.phoR9_ = 0.;
    int theGoodMatchType = -1;
    float theGoodMatchPt = -1.;
    if ( theGoodPhoton == 0 ) {
      outtree.phoEt_ = intree.pho1_.Et();
      outtree.phoEta_ = intree.pho1_.Eta();
      outtree.phoPhi_ = intree.pho1_.Phi();
      outtree.phoCombIso1_ = intree.phoCombIso1_a1_;
      outtree.phoCombIso2_ = intree.phoCombIso2_a1_;
      outtree.phoCombIso3_ = intree.phoCombIso3_a1_;
      outtree.phoWorstIso_ = std::max(intree.phoWorstIso_a1_-intree.rho_*GetChHadEA(outtree.phoEta_), (float)0);
      outtree.phoHasPixelSeed_ = intree.phoHasPixelSeed_a1_;
      outtree.phoR9_ =  intree.phoR9_a1_;
      if (isMC) theGoodMatchType =  intree.phoMatchType_a1_;
      if (isMC) theGoodMatchPt =  intree.phoMatchPt_a1_;
      outtree.phoSeedTime_ = intree.phoSeedTime_a1_;
      outtree.phoCoviEtaiEta_ = intree.phoCoviEtaiEta_a1_;
      outtree.phoCoviPhiiPhi_ = intree.phoCoviPhiiPhi_a1_;
      outtree.phoMipTotEnergy_ = intree.phoMipTotEnergy_a1_;
      outtree.phoRoundness_ = intree.phoRoundness_a1_;
      outtree.phoAngle_ = intree.phoAngle_a1_;
    }
    else if ( theGoodPhoton == 1 ) {
      outtree.phoEt_ = intree.pho2_.Et();
      outtree.phoEta_ = intree.pho2_.Eta();
      outtree.phoPhi_ = intree.pho2_.Phi();
      outtree.phoCombIso1_ = intree.phoCombIso1_a2_;
      outtree.phoCombIso2_ = intree.phoCombIso2_a2_;
      outtree.phoCombIso3_ = intree.phoCombIso3_a2_;
      outtree.phoWorstIso_ = std::max(intree.phoWorstIso_a2_-intree.rho_*GetChHadEA(outtree.phoEta_), (float)0);
      outtree.phoHasPixelSeed_ = intree.phoHasPixelSeed_a2_;
      outtree.phoR9_ =  intree.phoR9_a2_;
      if (isMC) theGoodMatchType =  intree.phoMatchType_a2_;
      if (isMC) theGoodMatchPt =  intree.phoMatchPt_a2_;
      outtree.phoSeedTime_ = intree.phoSeedTime_a2_;
      outtree.phoCoviEtaiEta_ = intree.phoCoviEtaiEta_a2_;
      outtree.phoCoviPhiiPhi_ = intree.phoCoviPhiiPhi_a2_;
      outtree.phoMipTotEnergy_ = intree.phoMipTotEnergy_a2_;
      outtree.phoRoundness_ = intree.phoRoundness_a2_;
      outtree.phoAngle_ = intree.phoAngle_a2_;
    }
    else if ( theGoodPhoton == 2 ) {
      outtree.phoEt_ = intree.pho3_.Et();
      outtree.phoEta_ = intree.pho3_.Eta();
      outtree.phoPhi_ = intree.pho3_.Phi();
      outtree.phoCombIso1_ = intree.phoCombIso1_a3_;
      outtree.phoCombIso2_ = intree.phoCombIso2_a3_;
      outtree.phoCombIso3_ = intree.phoCombIso3_a3_;
      outtree.phoWorstIso_ = std::max(intree.phoWorstIso_a3_-intree.rho_*GetChHadEA(outtree.phoEta_), (float)0);
      outtree.phoHasPixelSeed_ = intree.phoHasPixelSeed_a3_;
      outtree.phoR9_ =  intree.phoR9_a3_;
      if (isMC) theGoodMatchType =  intree.phoMatchType_a3_;
      if (isMC) theGoodMatchPt =  intree.phoMatchPt_a3_;
      outtree.phoSeedTime_ = intree.phoSeedTime_a3_;
      outtree.phoCoviEtaiEta_ = intree.phoCoviEtaiEta_a3_;
      outtree.phoCoviPhiiPhi_ = intree.phoCoviPhiiPhi_a3_;
      outtree.phoMipTotEnergy_ = intree.phoMipTotEnergy_a3_;
      outtree.phoRoundness_ = intree.phoRoundness_a3_;
      outtree.phoAngle_ = intree.phoAngle_a3_;
    }
    else
      cout << "Error in the selection function! theGoodPhoton==" << theGoodPhoton  << endl; 
    
    //only for Wjets and Zjets discard events if photon ISR
    if (isOverlapSample && theGoodMatchType == 0) continue;

    //if photon does not come from electron/ISR discard it (use fake rate)
    //if (isMC && (theGoodMatchType == 1 || theGoodMatchType == 2 || theGoodMatchType == 4)) continue;

    //now the scale factors which are pt dependent
    float thisKFactorWeight = thisScale;    
    if (thisScale < 0) { 
      //if the photon is matched to the GEN object the GEN pt, otherwise use reco
      float photonTruePt = outtree.phoEt_;
      if (theGoodMatchType == 0) photonTruePt = theGoodMatchPt;      
      thisKFactorWeight = KFactorWeight(thisScale, photonTruePt);
    }
    //for photo efficiencies
    float thisScaleFactorWeight = 1.;
    if (isMC) thisScaleFactorWeight = ScaleFactorWeight(outtree.phoR9_, outtree.phoEt_);
    //if the photon comes from an electron use the fake rate scale factor
    //if (isMC && theGoodMatchType == 3) thisScaleFactorWeight *= 1.5;    
    outtree.kf_weight_  = thisKFactorWeight;
    outtree.hlt_weight_ = thisScaleFactorWeight;
      
    if (fInputSelection == "Lepton") {
      TLorentzVector tempBos;
      tempBos.SetPtEtaPhiE(intree.metCor_,0.,intree.metCorPhi_,intree.metCor_);
      TLorentzVector tempLep;
      tempLep.SetPtEtaPhiE(intree.lep1_.Pt(),0.,intree.lep1_.Phi(),intree.lep1_.Et());
      tempBos += tempLep;
      outtree.bosonPt_ = tempBos.Pt();
      outtree.bosonEta_ = tempBos.Eta();
      outtree.bosonPhi_ = tempBos.Phi();
      //Return transverse mass for the W
      float thisDphi = GetCorrDeltaPhi(outtree.metCorPhi_, intree.lep1_.Phi());
      outtree.bosonMass_ = sqrt(2*intree.metCor_*intree.lep1_.Pt()*(1-TMath::Cos(thisDphi)));
      //W case is easy metBosCor is W mom in transverse plane
      outtree.metBosCor_ = tempBos.Pt();
      outtree.metBosCorPhi_ = tempBos.Phi();
    }
    // Easy case muons
    if (fInputSelection == "DiLepton" && intree.lep1_.M() > 0.05) {
      TLorentzVector tempBos;
      TLorentzVector tempLep1;
      tempLep1.SetPtEtaPhiE(intree.lep1_.Pt(),intree.lep1_.Eta(),intree.lep1_.Phi(),intree.lep1_.E());
      TLorentzVector tempLep2;
      tempLep2.SetPtEtaPhiE(intree.lep2_.Pt(),intree.lep2_.Eta(),intree.lep2_.Phi(),intree.lep2_.E());
      tempBos = tempLep1 + tempLep2;
      outtree.bosonPt_ = tempBos.Pt();
      outtree.bosonEta_ = tempBos.Eta();
      outtree.bosonPhi_ = tempBos.Phi();
      outtree.bosonMass_ = tempBos.M();
      TLorentzVector tempPho;
      tempPho.SetPtEtaPhiM(outtree.phoEt_, outtree.phoEta_, outtree.phoPhi_, 0.);
      tempPho += tempBos;
      outtree.bosonPhoMass_ = tempPho.M();
      //Z case metBosCor is more tricky
      TLorentzVector tempMet;
      tempMet.SetPtEtaPhiE(intree.metCor_,0.,intree.metCorPhi_,intree.metCor_);
      tempMet += tempBos;
      outtree.metBosCor_ = tempBos.Pt();
      outtree.metBosCorPhi_ = tempBos.Phi();
    }
    // Some fun with electrons!
    if (fInputSelection == "DiLepton" && intree.lep1_.Pt() > 0. && intree.lep1_.M() < 0.05) {
      TLorentzVector tempBos;
      TLorentzVector tempLep1;
      tempLep1.SetPtEtaPhiE(intree.lep1_.Pt(),intree.lep1_.Eta(),intree.lep1_.Phi(),intree.lep1_.E());
      TLorentzVector tempLep2;
      // Lep 2 might be in fact an electron reconstructed as photon, otherwise second lepton
      if (intree.pho2_.Pt() > 0) {
        tempLep2.SetPtEtaPhiE(intree.pho2_.Pt(),intree.pho2_.Eta(),intree.pho2_.Phi(),intree.pho2_.E());
        if (intree.phoPassEleVeto_a2_ > 0.5) outtree.lepAsPho_ = 0;
        else outtree.lepAsPho_ = 1;
      }
      else 
        tempLep2.SetPtEtaPhiE(intree.lep2_.Pt(),intree.lep2_.Eta(),intree.lep2_.Phi(),intree.lep2_.E());
        
      tempBos = tempLep1 + tempLep2;
      outtree.bosonPt_ = tempBos.Pt();
      outtree.bosonEta_ = tempBos.Eta();
      outtree.bosonPhi_ = tempBos.Phi();
      outtree.bosonMass_ = tempBos.M();
      TLorentzVector tempPho;
      tempPho.SetPtEtaPhiM(outtree.phoEt_, outtree.phoEta_, outtree.phoPhi_, 0.);
      tempPho += tempBos;
      outtree.bosonPhoMass_ = tempPho.M();
      //Z case metBosCor is more tricky
      TLorentzVector tempMet;
      tempMet.SetPtEtaPhiE(intree.metCor_,0.,intree.metCorPhi_,intree.metCor_);
      tempMet += tempBos;
      outtree.metBosCor_ = tempBos.Pt();
      outtree.metBosCorPhi_ = tempBos.Phi();
    }
    
    outtree.nalljets_ = intree.njets_;
    //now we will loop on the jets and clean the collection wrt the sel photon: necessary for QCD tree
    float dphi,deta;
    for (unsigned int ijet = 0; ijet < intree.njets_; ijet++) {
      if (ijet == 0) {
        dphi = GetCorrDeltaPhi(outtree.phoPhi_, intree.jet1_.Phi());
        deta = outtree.phoEta_ - intree.jet1_.Eta();
        if (dphi*dphi+deta*deta < 0.3*0.3) {
          outtree.nalljets_ --;
          continue;
        }
        outtree.jet1Pt_ = intree.jet1_.Pt();
        outtree.jet1Eta_ = intree.jet1_.Eta();
        outtree.jet1Phi_ = intree.jet1_.Phi();
      }
      else if (ijet == 1) {
        dphi = GetCorrDeltaPhi(outtree.phoPhi_, intree.jet2_.Phi());
        deta = outtree.phoEta_ - intree.jet2_.Eta();
        if (dphi*dphi+deta*deta < 0.3*0.3) {
          outtree.nalljets_ --;
          continue;
        }
        outtree.jet2Pt_ = intree.jet2_.Pt();
        outtree.jet2Eta_ = intree.jet2_.Eta();
        outtree.jet2Phi_ = intree.jet2_.Phi();
      }
      else if (ijet == 2) {
        dphi = GetCorrDeltaPhi(outtree.phoPhi_, intree.jet3_.Phi());
        deta = outtree.phoEta_ - intree.jet3_.Eta();
        if (dphi*dphi+deta*deta < 0.3*0.3) {
          outtree.nalljets_ --;
          continue;
        }
        outtree.jet3Pt_ = intree.jet3_.Pt();
        outtree.jet3Eta_ = intree.jet3_.Eta();
        outtree.jet3Phi_ = intree.jet3_.Phi();
      }
      else if (ijet == 3) {
        dphi = GetCorrDeltaPhi(outtree.phoPhi_, intree.jet4_.Phi());
        deta = outtree.phoEta_ - intree.jet4_.Eta();
        if (dphi*dphi+deta*deta < 0.3*0.3) {
          outtree.nalljets_ --;
          continue;
        }
        outtree.jet4Pt_ = intree.jet4_.Pt();
        outtree.jet4Eta_ = intree.jet4_.Eta();
        outtree.jet4Phi_ = intree.jet4_.Phi();
      }
      else break;
    }

    //compute the min met
    float theMinMet, theMinMetProb;
    ComputeMinMet(intree, outtree, theMinMet, theMinMetProb);
    outtree.metMin_ = theMinMet;
    outtree.metMinProb_ = theMinMetProb;

    outtree.phoMetDeltaPhi_ = GetCorrDeltaPhi(outtree.phoPhi_, intree.metCorPhi_);
    if (outtree.nalljets_ > 0) outtree.jetMetDeltaPhi_ = GetCorrDeltaPhi(outtree.jet1Phi_, intree.metCorPhi_);
    if (outtree.nalljets_ > 0) outtree.phoJetDeltaPhi_ = GetCorrDeltaPhi(outtree.phoPhi_, outtree.jet1Phi_);
  
    //leptons
    outtree.nlep_ = intree.nlep_;
    for (unsigned int ilep = 0; ilep < intree.nlep_; ilep++) {
      if (ilep == 0) {
        outtree.lep1Pt_ = intree.lep1_.Pt();
        outtree.lep1Eta_ = intree.lep1_.Eta();
        outtree.lep1Phi_ = intree.lep1_.Phi();
        outtree.lep1Mass_ = intree.lep1_.M();
      }
      else if (ilep == 1) {
        outtree.lep2Pt_ = intree.lep2_.Pt();
        outtree.lep2Eta_ = intree.lep2_.Eta();
        outtree.lep2Phi_ = intree.lep2_.Phi();
        outtree.lep2Mass_ = intree.lep2_.M();
      }
      else break;
    }

    //cosmics
    outtree.ncosmics_ = intree.ncosmics_;
    if (intree.ncosmics_ > 0) {
      outtree.cosmic1Pt_ = intree.cosmic1_.Pt();
      outtree.cosmic1Eta_= intree.cosmic1_.Eta();
      outtree.cosmic1Phi_= intree.cosmic1_.Phi();
    }

    //fill the outtree content
    outtree.tree_->Fill();    
  }
    
  fOutFile   -> cd();
  outtree.tree_ -> Write();  
  return;
}

//--------------------------------------------------------------------------------------------------
void TreeReducer::InitResoFunctions()
{
  fMexReso = new TF1("Mex Sigma", "[0]*TMath::Sqrt(x) + [1]",0,5000);
  fMeyReso = new TF1("Mey Sigma", "[0]*TMath::Sqrt(x) + [1]",0,5000); 
  fPhotonReso = new TF1("Photon_Sigma","[0]*pow(x,[1])",0,10000);
  fMuonReso = new TF1("ptsigma_barrel", "pol3", 0, 10000);
  
  for (int iFunc = 0; iFunc < 12; iFunc++) {
    TString stringNum = "";
    stringNum += iFunc;
    fJetReso[iFunc] = new TF1 ("jetReso_"+stringNum,"TMath::Sqrt(((sign([0])*(([0]/x)^2))+(([1]^2)*(x^([2]-1.))))+([3]^2))",0,15000);
  }
  
}

//--------------------------------------------------------------------------------------------------
void TreeReducer::SetResoFunctions(bool isData)
{
  if(!isData)   fMexReso->SetParameters(0.53,3.08);
  if(!isData)   fMeyReso->SetParameters(0.53,3.4);
  if(isData)    fMexReso->SetParameters(0.61,0.37);
  if(isData)    fMeyReso->SetParameters(0.62,0.17);

  fPhotonReso->SetParameters(0.1138,-0.449);
  fMuonReso->SetParameters(0.006821,0.0001047,-1.213e-08,1.874e-12);  

  fJetReso[0]->SetParameters(2.866,0.3118,0.4075,0.01823);
  fJetReso[1]->SetParameters(2.91,0.2793,0.4629,0.001049);
  fJetReso[2]->SetParameters(2.768,0.3797,0.3144,0.02803);
  fJetReso[3]->SetParameters(2.934,0.3251,0.4401,0.0079); 
  fJetReso[4]->SetParameters(2.617,0.736,0.0899,-0.04179);
  fJetReso[5]->SetParameters(0.1406,1.477,-0.2062,-0.03656); 
  fJetReso[6]->SetParameters(1.959,1.099,-0.1357,-0.02382); 
  fJetReso[7]->SetParameters(4.113,0.4146,0.1918,0.02413); 
  fJetReso[8]->SetParameters(5.817,0.1547,0.5529,0.001136); 
  fJetReso[9]->SetParameters(4.894,0.3666,0.4251,-0.00215); 
  fJetReso[10]->SetParameters(3.624,0.2542,0.6046,0.02232); 
  fJetReso[11]->SetParameters(2.727,1.035,-0.1662,0);
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::PUWeight(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!sPUWeights)
    return 1.0;
  
  return sPUWeights->GetBinContent(sPUWeights->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::PUWeightUp(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!sPUWeightsUp)
    return 1.0;
  
  return sPUWeightsUp->GetBinContent(sPUWeightsUp->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::PUWeightDown(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!sPUWeightsDown)
    return 1.0;
  
  return sPUWeightsDown->GetBinContent(sPUWeightsDown->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::KFactorWeight(Float_t scale, Float_t phet)
{
  //k-factors description here
  //https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=283118
  //use linear interpolation
  //Zgamma
  if (scale == -1) {
    if (phet < 145.)
      return 1.0;
    else 
      return gZGammaKFactor -> Eval(phet);
  }
  //Wgamma
  else if (scale == -2) {
    if (phet < 145.)
      return 1.0;
    else 
      return gWGammaKFactor -> Eval(phet);
  }
  //QCD, taken from dedicated studies
  else if (scale == -3) {
    //use linear interpolation of fake rate
    return 
      gJetFakeRate -> Eval(phet);
  }
  //Ele, taken from dedicated studies
  else if (scale == -4) {
    //use linear interpolation of fake rate
    return 
      gEleFakeRate -> Eval(phet);
  }
  else
    return 1.;
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::ScaleFactorWeight(Float_t R9, Float_t phet)
{
  //photon ID efficiencies, taken from AN-12-160
  float preselEff = 1.;
  float idEff = 1.;
  //timing is taken from Zee studies
  float timingEff = 0.9631;
  
  if (R9 < 0.9) preselEff = 0.996*0.992;
  else preselEff = 0.998*0.999;
  if (R9 < 0.94) idEff = 0.992;
  else idEff = 1.002;
  
  float theEff = preselEff*idEff*timingEff;
  return theEff;
}


//--------------------------------------------------------------------------------------------------
bool TreeReducer::EventIsSelected(MitGPTree &tree, int treeType, int& theGoodPhoton, bool isEleFake)
{
  //only for ele fake rate studies
  int theEleVeto = !isEleFake;
  //standard selection
  if (treeType == 0) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == theEleVeto && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == theEleVeto && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == theEleVeto && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 1.5) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //lepton selection
  if (treeType == 1) {
    //photon number and lepton number
    if ( tree.nphotons_ == 0 || tree.nlep_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == theEleVeto && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == theEleVeto && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == theEleVeto && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 1.5) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //QCD selection:you want to use this also for the lepton selection, so relax met cut
  else if (treeType == 2) {
    //met and photon number
    if ( tree.nphotons_ == 0 ) return false;
    //Phase space & fake photons (no sigmaieta, sideband iso)
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == theEleVeto && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5 &&
      PhotonIsFake(tree.phoR9_a1_, tree.phoHadOverEm_a1_, tree.phoCoviEtaiEta_a1_, tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == theEleVeto && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5 &&
      PhotonIsFake(tree.phoR9_a2_, tree.phoHadOverEm_a2_, tree.phoCoviEtaiEta_a2_, tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == theEleVeto && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 1.5 &&
      PhotonIsFake(tree.phoR9_a3_, tree.phoHadOverEm_a3_, tree.phoCoviEtaiEta_a3_, tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //Beam halo selection: really do not care if you have leptons
  else if (treeType == 3) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == theEleVeto && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5 &&
      PhotonIsSelected(tree.phoR9_a1_, tree.phoHadOverEm_a1_, -1., tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == theEleVeto && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5 &&
      PhotonIsSelected(tree.phoR9_a2_, tree.phoHadOverEm_a2_, -1., tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == theEleVeto && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 1.5 &&
      PhotonIsSelected(tree.phoR9_a3_, tree.phoHadOverEm_a3_, -1., tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  else 
    return false;
}

//--------------------------------------------------------------------------------------------------
bool TreeReducer::EventIsSelectedEG(MitGPTree &tree, int treeType, int& theGoodPhoton, bool isEleFake)
{
  //only for ele fake rate studies
  int hasPixelSeed = isEleFake;
  //standard selection
  if (treeType == 0) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoHasPixelSeed_a1_ == hasPixelSeed && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -3 && tree.phoSeedTime_a1_ < 3) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoHasPixelSeed_a2_ == hasPixelSeed && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -3 && tree.phoSeedTime_a2_ < 3) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoHasPixelSeed_a3_ == hasPixelSeed && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -3 && tree.phoSeedTime_a3_ < 3) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //lepton selection
  if (treeType == 1) {
    //photon number and lepton number
    if ( tree.nphotons_ == 0 || tree.nlep_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoHasPixelSeed_a1_ == hasPixelSeed && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -3 && tree.phoSeedTime_a1_ < 3) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoHasPixelSeed_a2_ == hasPixelSeed && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -3 && tree.phoSeedTime_a2_ < 3) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoHasPixelSeed_a3_ == hasPixelSeed && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -3 && tree.phoSeedTime_a3_ < 3) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //QCD selection:you want to use this also for the lepton selection, so relax met cut
  else if (treeType == 2) {
    //met and photon number
    if ( tree.nphotons_ == 0 ) return false;
    //Phase space & fake photons (no sigmaieta, no pixel seed/ele veto, sideband iso)
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoHasPixelSeed_a1_ == hasPixelSeed && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -3 && tree.phoSeedTime_a1_ < 3 &&
      PhotonIsFakeEG(tree.pho1_.Et(), tree.phoHadOverEm_a1_, tree.phoCoviEtaiEta_a1_, tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoHasPixelSeed_a2_ == hasPixelSeed && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -3 && tree.phoSeedTime_a2_ < 3 &&
      PhotonIsFakeEG(tree.pho2_.Et(), tree.phoHadOverEm_a2_, tree.phoCoviEtaiEta_a2_, tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoHasPixelSeed_a3_ == hasPixelSeed && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -3 && tree.phoSeedTime_a3_ < 3 &&
      PhotonIsFakeEG(tree.pho3_.Et(), tree.phoHadOverEm_a3_, tree.phoCoviEtaiEta_a3_, tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //Beam halo selection: really do not care if you have leptons
  else if (treeType == 3 && !fIsBHStudy) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoHasPixelSeed_a1_ == hasPixelSeed && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -3 && tree.phoSeedTime_a1_ < 3 &&
      PhotonIsSelectedEG(tree.pho1_.Et(), tree.phoHadOverEm_a1_, -1., tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoHasPixelSeed_a2_ == hasPixelSeed && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -3 && tree.phoSeedTime_a2_ < 3 &&
      PhotonIsSelectedEG(tree.pho2_.Et(), tree.phoHadOverEm_a2_, -1., tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoHasPixelSeed_a3_ == hasPixelSeed && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -3 && tree.phoSeedTime_a3_ < 3 &&
      PhotonIsSelectedEG(tree.pho3_.Et(), tree.phoHadOverEm_a3_, -1., tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //Beam halo study: drop pixel seed, sigma ieta cut and beam halo sel
  else if (treeType == 3 && fIsBHStudy) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoSeedTime_a1_ > -3 && tree.phoSeedTime_a1_ < 3 &&
      PhotonIsSelectedEG(tree.pho1_.Et(), tree.phoHadOverEm_a1_, -1., tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoSeedTime_a2_ > -3 && tree.phoSeedTime_a2_ < 3 &&
      PhotonIsSelectedEG(tree.pho2_.Et(), tree.phoHadOverEm_a2_, -1., tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoSeedTime_a3_ > -3 && tree.phoSeedTime_a3_ < 3 &&
      PhotonIsSelectedEG(tree.pho3_.Et(), tree.phoHadOverEm_a3_, -1., tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  else 
    return false;
}


//--------------------------------------------------------------------------------------------------
bool TreeReducer::PhotonIsSelected(float R9, float HoverE, float CoviEtaiEta, float Iso1, float Iso2, float Iso3)
{
  int Cat ;
  Cat = 0;
  if (R9 < 0.94 ) Cat = 1;

  //Isolation cuts

  if (Cat == 0) {
    if (Iso1 > 6.0) return false ;
    if (Iso2 > 10.0) return false ;
    if (Iso3 > 3.8) return false ;
    if (CoviEtaiEta > 0.0108) return false ;
    if (HoverE > 0.124) return false ;
    if (R9 < 0.94 ) return false ; 
    return true;
  }
  else if (Cat == 1) {
    if (Iso1 > 4.7) return false ;
    if (Iso2 > 6.5) return false ;
    if (Iso3 > 2.5) return false ;
    if (CoviEtaiEta > 0.0102) return false ;
    if (HoverE > 0.092) return false ;
    if (R9 < 0.28 ) return false ; 
    return true;
  }
  else
    return false;
}

//--------------------------------------------------------------------------------------------------
bool TreeReducer::PhotonIsSelectedEG(float Pt, float HoverE, float CoviEtaiEta, float Iso1, float Iso2, float Iso3)
{

  //Isolation cuts
  if (Iso1 > 1.5) return false ;
  if (Iso2 > 1.0 + 0.04 * Pt) return false ;
  if (Iso3 > 0.7 + 0.005 * Pt) return false ;
  //Other cuts
  if (CoviEtaiEta > 0.011 ) return false ;
  if (HoverE > 0.05 ) return false ;
  return true;

}

//--------------------------------------------------------------------------------------------------
bool TreeReducer::PhotonIsFake(float R9, float HoverE, float CoviEtaiEta, float Iso1, float Iso2, float Iso3)
{
  int Cat ;
  Cat = 0;
  if (R9 < 0.94 ) Cat = 1;

  //Exclude sigma ieta far sideband
  if (CoviEtaiEta > 0.014) 
    return false;

  //Isolation cuts: first exclude far sideband
  if (Iso1 > 18.0 || Iso2 > 30.0 || Iso3 > 11.4)
    return false;

  if (Cat == 0) {
    // exclude signal region
    if (Iso1 <= 2.0 && Iso2 <= 3.0 && Iso3 <= 1.0) return false ;
    if (HoverE > 0.124) return false ;
    if (R9 < 0.94 ) return false ; 
    return true;
  }
  else if (Cat == 1) {
    // exclude signal region
    if (Iso1 <= 2.0 && Iso2 <= 3.0 && Iso3 <= 1.0) return false ;
    if (HoverE > 0.092) return false ;
    if (R9 < 0.28 ) return false ; 
    return true;
  }
  else
    return false;
}

//--------------------------------------------------------------------------------------------------
bool TreeReducer::PhotonIsFakeEG(float Pt, float HoverE, float CoviEtaiEta, float Iso1, float Iso2, float Iso3)
{
  
  //Exclude sigma ieta far sideband
  if (CoviEtaiEta > 0.013) 
    return false;

  //Exclude H/E far sideband
  if (HoverE > 0.05) 
    return false;

  //Isolation cuts: first exclude far sideband
  float Iso1Cut = std::min(5*2.6, 0.2*Pt);
  float Iso2Cut = std::min(5*(3.5+0.04*Pt), 0.2*Pt);
  float Iso3Cut = std::min(5*(1.3+0.005*Pt), 0.2*Pt);
  if (Iso1 > Iso1Cut || Iso2 > Iso2Cut || Iso3 > Iso3Cut)
    return false;

  // exclude signal region
  if (Iso1 <= 2.6 && Iso2 <= 3.5+0.04*Pt && Iso3 <= 1.3+0.005*Pt) 
    return false ;

  return true;

}

//--------------------------------------------------------------------------------------------------
float TreeReducer::GetCorrDeltaPhi(float phi1, float phi2)
{
  float corrDeltaPhi = TMath::Abs(phi1 - phi2);
  if (corrDeltaPhi > TMath::Pi())
    corrDeltaPhi = TMath::TwoPi() - corrDeltaPhi;     
  return corrDeltaPhi;
}

//--------------------------------------------------------------------------------------------------
void TreeReducer::ComputeMinMet(MitGPTree &intree, MitGPTreeReduced &outtree, float& theMinMet, float& theMinMetProb)
{
  float STPrime = intree.sumEt_;
  vector<float> reco_pt;
  vector<float> reco_phi;
  vector<float> sigma_ms;
  
  //first take care of the photon
  reco_pt.push_back(outtree.phoEt_);
  reco_phi.push_back(outtree.phoPhi_);
  sigma_ms.push_back(fPhotonReso->Eval(outtree.phoEt_)*outtree.phoEt_);
  STPrime -= outtree.phoEt_;

  //leptons
  for (unsigned int ilep = 0; ilep < intree.nlep_; ilep++) {
    if (ilep == 0) {
      reco_pt.push_back(intree.lep1_.Pt());
      reco_phi.push_back(intree.lep1_.Phi());
      STPrime -= intree.lep1_.Pt();
      if (intree.lid1_ == 13)   
        sigma_ms.push_back(fMuonReso->Eval(intree.lep1_.Pt())*intree.lep1_.Pt());
      else    
        sigma_ms.push_back(fPhotonReso->Eval(intree.lep1_.Pt())*intree.lep1_.Pt());
    }
    else if (ilep == 1) {
      reco_pt.push_back(intree.lep2_.Pt());
      reco_phi.push_back(intree.lep2_.Phi());
      STPrime -= intree.lep2_.Pt();
      if (intree.lid2_ == 13)   
        sigma_ms.push_back(fMuonReso->Eval(intree.lep2_.Pt())*intree.lep2_.Pt());
      else    
        sigma_ms.push_back(fPhotonReso->Eval(intree.lep2_.Pt())*intree.lep2_.Pt());
    }
    else if (ilep == 2) {
      reco_pt.push_back(intree.lep3_.Pt());
      reco_phi.push_back(intree.lep3_.Phi());
      STPrime -= intree.lep3_.Pt();
      if (intree.lid3_ == 13)   
        sigma_ms.push_back(fMuonReso->Eval(intree.lep3_.Pt())*intree.lep3_.Pt());
      else    
        sigma_ms.push_back(fPhotonReso->Eval(intree.lep3_.Pt())*intree.lep3_.Pt());
    }
    else break;
  }

  //jets
  for (unsigned int ijet = 0; ijet < outtree.nalljets_; ijet++) {
    if (ijet == 0) {
      reco_pt.push_back(outtree.jet1Pt_);
      reco_phi.push_back(outtree.jet1Phi_);
      STPrime -= outtree.jet1Pt_;
      sigma_ms.push_back(GetJetReso(outtree.jet1Pt_, outtree.jet1Eta_, outtree.isData_));
    }
    else if (ijet == 1) {
      reco_pt.push_back(outtree.jet2Pt_);
      reco_phi.push_back(outtree.jet2Phi_);
      STPrime -= outtree.jet2Pt_;
      sigma_ms.push_back(GetJetReso(outtree.jet2Pt_, outtree.jet2Eta_, outtree.isData_));
    }
    else if (ijet == 2) {
      reco_pt.push_back(outtree.jet3Pt_);
      reco_phi.push_back(outtree.jet3Phi_);
      STPrime -= outtree.jet3Pt_;
      sigma_ms.push_back(GetJetReso(outtree.jet3Pt_, outtree.jet3Eta_, outtree.isData_));
    }
    else if (ijet == 3) {
      reco_pt.push_back(outtree.jet4Pt_);
      reco_phi.push_back(outtree.jet4Phi_);
      STPrime -= outtree.jet4Pt_;
      sigma_ms.push_back(GetJetReso(outtree.jet4Pt_, outtree.jet4Eta_, outtree.isData_));
    }
    else break;
  }
  
  //This is the Missing Ex, Ey resolution
  float sigma_mex = fMexReso->Eval(TMath::Abs(STPrime));
  float sigma_mey = fMeyReso->Eval(TMath::Abs(STPrime));

  unsigned int Ndim = reco_pt.size();  
  
  //Set up the minimizer!
  ROOT::Math::Minimizer* min =  ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(100000);
  min->SetTolerance(0.001);

  //The function is called here (the function is at my .h file and i'll put the definition in this mail as well)
  SetMinMetParam(reco_pt, reco_phi, sigma_ms, sigma_mex, sigma_mey);
  ROOT::Math::Functor f(this, &TreeReducer::MinMetFunc, Ndim);
  min->SetFunction(f);

  // Set the free variables to be minimized! upper and lower boundaries.                                                                                                                                  
  for (int i=0; i<(int) Ndim; i++){
    TString name = "pt";
    name += i;
    double lower = TMath::Max((float) 0.,reco_pt[i] - 3*sigma_ms[i]);
    double upper = reco_pt[i] + (3*sigma_ms[i]);
    //std::cout << "lower: " << lower << " upper: " << upper << std::endl;                                                                                                      
    min->SetLimitedVariable(i,name.Data(),reco_pt[i],0.1,lower,upper);
  }

  min->Minimize();  
  theMinMetProb = TMath::Prob(MinMetFunc(min->X()),2);

  // Finding Parametrized MET:                                                                                                                                                  
  double ParPx = 0.;
  double ParPy = 0.;
  double RecoPx = 0.;
  double RecoPy = 0.;

  for(int j=0; j<(int) Ndim; j++){
    ParPx += min->X()[j]*cos(reco_phi[j]);
    ParPy += min->X()[j]*sin(reco_phi[j]);
    RecoPx += reco_pt[j]*cos(reco_phi[j]);
    RecoPy += reco_pt[j]*sin(reco_phi[j]);
  }

  theMinMet = sqrt((ParPx*ParPx)+(ParPy*ParPy));

  // Clean up memory allocation
  delete min;
  
  return;
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::GetJetReso(float pt, float eta, bool isData)
{
  float scale = 1.;
  int index = -1; //index for eta binning
  
  if (eta <= 0.3)              {index = 0; scale = 1.052;}
  if (eta > 0.3 && eta <= 0.5) {index = 1; scale = 1.052;}
  if (eta > 0.5 && eta <= 0.8) {index = 2; scale = 1.057;}
  if (eta > 0.8 && eta <= 1.1) {index = 3; scale = 1.057;}
  if (eta > 1.1 && eta <= 1.4) {index = 4; scale = 1.096;}
  if (eta > 1.4 && eta <= 1.7) {index = 5; scale = 1.096;}
  if (eta > 1.7 && eta <= 2.0) {index = 6; scale = 1.134;}
  if (eta > 2.0 && eta <= 2.3) {index = 7; scale = 1.134;}
  if (eta > 2.3 && eta <= 2.8) {index = 8; scale = 1.288;}
  if (eta > 2.8 && eta <= 3.2) {index = 9; scale = 1.288;}
  if (eta > 3.2 && eta <= 4.1) {index = 10; scale = 1.288;}
  if (eta > 4.1)               {index = 11; scale = 1.288;}
 
  float reso = -1;
  if(!isData) reso = pt * fJetReso[index]->Eval(pt);
  else        reso = pt * scale * fJetReso[index]->Eval(pt);
  return reso;
}

//--------------------------------------------------------------------------------------------------
void TreeReducer::SetMinMetParam(std::vector<float>& reco_pt, std::vector<float>& reco_phi, std::vector<float>& sigma_ms, float sigma_mex, float sigma_mey)
{
  MinMetRecoPtVec.clear();
  MinMetRecoPhiVec.clear();
  MinMetSigmaVec.clear();
  MinMetSigmaMex = sigma_mex;
  MinMetSigmaMey = sigma_mey;

  int Ndim = reco_pt.size();
  for(int i=0; i<Ndim; i++){
    MinMetRecoPtVec.push_back(reco_pt[i]);
    MinMetRecoPhiVec.push_back(reco_phi[i]);
    MinMetSigmaVec.push_back(sigma_ms[i]);
  }
}

//--------------------------------------------------------------------------------------------------
double TreeReducer::MinMetFunc(const double* par)
{
  int Ndim = MinMetRecoPtVec.size();
  double px = 0, py =0, arg = 0;
  for(int i=0; i<Ndim; i++){
    px += par[i]*cos(MinMetRecoPhiVec[i]);
    py += par[i]*sin(MinMetRecoPhiVec[i]);
    arg += pow((MinMetRecoPtVec[i]-par[i])/(MinMetSigmaVec[i]),2);
  }
  return arg + ((px*px)/(MinMetSigmaMex*MinMetSigmaMex) + (py*py)/(MinMetSigmaMey*MinMetSigmaMey));
}

//--------------------------------------------------------------------------------------------------
float TreeReducer::GetChHadEA(float eta)
{

  if (eta  < 1)                        return  0.012;
  else if (eta >= 1 && eta < 1.479)    return  0.010;
  else if (eta >= 1.479 && eta < 2.0)  return  0.014;
  else if (eta >= 2.0 && eta < 2.2)    return  0.012;
  else if (eta >= 2.2 && eta < 2.3)    return  0.016;
  else if (eta >= 2.3 && eta < 2.4)    return  0.020;
  else if (eta >= 2.4)                    return  0.012;
  else return 0;
  
}
