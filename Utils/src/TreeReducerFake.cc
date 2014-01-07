#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TPad.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TBox.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TTree.h>
#include "MitMonoPhoton/Utils/interface/TreeReducerFake.h"
#include "MitMonoPhoton/Core/MitGPTree.h"
#include "MitMonoPhoton/Core/MitGPTreeReducedFake.h"

#include "TLorentzVector.h"

ClassImp(mithep::TreeReducerFake)

using namespace std;
using namespace mithep;

const TH1D *TreeReducerFake::sPUWeights = 0;
const TH1D *TreeReducerFake::sPUWeightsUp = 0;
const TH1D *TreeReducerFake::sPUWeightsDown = 0;

//--------------------------------------------------------------------------------------------------
TreeReducerFake::TreeReducerFake(const Sample *mySample) :
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
TreeReducerFake::~TreeReducerFake()
{
  // Destructor
}

//--------------------------------------------------------------------------------------------------
void TreeReducerFake::MakeTree()
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
  bool isMC = false, isData = false;
  int treeType = 2; // Fake tree 
  float thisXsec   = *fSample->Xsec();
  float thisScale  = *fSample->Scale();
  if (thisXsec > 0)         isMC = true;

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
  MitGPTreeReducedFake outtree;
  TString outTreeName = (*fSample->Name());
  outtree.CreateTree(outTreeName.Data());

  // set up the input tree
  MitGPTree intree;
  intree.LoadTree(inFileName,treeType);
  intree.InitTree(0);

  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {

    outtree.InitVariables();
    intree.tree_-> GetEntry(iEntry);
    //apply the selection
    //define the good photon 
    int theGoodPhoton = -1;
    // photonIsFake flag: 
    // 0 not signal not fake
    // 1 is signal not fake
    // 2 is fake not signal
    // 3 is fake and signal 
    int photonIsFake = -1;    
    if ( !EventIsSelected(intree, treeType, theGoodPhoton, photonIsFake) ) continue;
    
    //get this event weights
    float thisPUWeight = 1.;
    float thisPUWeightUp = 1.;
    float thisPUWeightDown = 1.;
    if ( isMC ) thisPUWeight = PUWeight(intree.npu_);
    if ( isMC ) thisPUWeightUp = PUWeightUp(intree.npu_);
    if ( isMC ) thisPUWeightDown = PUWeightDown(intree.npu_);
    float thisKFactorWeight = thisScale;
    if (thisScale < 0) thisKFactorWeight = KFactorWeight(thisScale, intree.pho1_.Et());
    
    //get the relevant variables
    outtree.isData_ = isData;
    
    outtree.evt_weight_ = baseEvtWeight;
    outtree.hlt_weight_ = 1.;
    outtree.pu_weight_  = thisPUWeight;
    outtree.pu_weightup_  = thisPUWeightUp;
    outtree.pu_weightdo_  = thisPUWeightDown;
    outtree.kf_weight_  = thisKFactorWeight;
  
    outtree.nvtx_  = intree.nvtx_;
    outtree.metCor_= intree.metCor_;
    outtree.metCorPhi_= intree.metCorPhi_;
    outtree.metSig_= intree.metSig_;
  
    outtree.nphotons_ = intree.nphotons_;
    if ( theGoodPhoton == 0 ) {
      outtree.phoEt_ = intree.pho1_.Et();
      outtree.phoEta_ = intree.pho1_.Eta();
      outtree.phoPhi_ = intree.pho1_.Phi();
      outtree.phoCombIso1_ = intree.phoCombIso1_a1_;
      outtree.phoCombIso2_ = intree.phoCombIso2_a1_;
      outtree.phoCombIso3_ = intree.phoCombIso3_a1_;
      outtree.phoCoviEtaiEta_ = intree.phoCoviEtaiEta_a1_;
      outtree.phoR9_          = intree.phoR9_a1_;
      outtree.phoIsFake_      = photonIsFake;
    }
    else if ( theGoodPhoton == 1 ) {
      outtree.phoEt_ = intree.pho2_.Et();
      outtree.phoEta_ = intree.pho2_.Eta();
      outtree.phoPhi_ = intree.pho2_.Phi();
      outtree.phoCombIso1_ = intree.phoCombIso1_a2_;
      outtree.phoCombIso2_ = intree.phoCombIso2_a2_;
      outtree.phoCombIso3_ = intree.phoCombIso3_a2_;
      outtree.phoCoviEtaiEta_ = intree.phoCoviEtaiEta_a2_;
      outtree.phoR9_          = intree.phoR9_a2_;
      outtree.phoIsFake_      = photonIsFake;
    }
    else if ( theGoodPhoton == 2 ) {
      outtree.phoEt_ = intree.pho3_.Et();
      outtree.phoEta_ = intree.pho3_.Eta();
      outtree.phoPhi_ = intree.pho3_.Phi();
      outtree.phoCombIso1_ = intree.phoCombIso1_a3_;
      outtree.phoCombIso2_ = intree.phoCombIso2_a3_;
      outtree.phoCombIso3_ = intree.phoCombIso3_a3_;
      outtree.phoCoviEtaiEta_ = intree.phoCoviEtaiEta_a3_;
      outtree.phoR9_          = intree.phoR9_a3_;
      outtree.phoIsFake_      = photonIsFake;
    }
    else
      cout << "Error in the selection function! theGoodPhoton==" << theGoodPhoton  << endl; 
          
    outtree.nalljets_ = intree.njets_;
    //now we will loop on the jets and clean the collection wrt the sel photon
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

    outtree.phoMetDeltaPhi_ = GetCorrDeltaPhi(outtree.phoPhi_, intree.metCorPhi_);
    if (outtree.nalljets_ > 0) outtree.jetMetDeltaPhi_ = GetCorrDeltaPhi(outtree.jet1Phi_, intree.metCorPhi_);
    if (outtree.nalljets_ > 0) outtree.phoJetDeltaPhi_ = GetCorrDeltaPhi(outtree.phoPhi_, outtree.jet1Phi_);

    //leptons
    outtree.nlep_ = intree.nlep_;
  
    //fill the outtree content
    outtree.tree_->Fill();    
  }
    
  fOutFile   -> cd();
  outtree.tree_ -> Write();  
  return;
}

//--------------------------------------------------------------------------------------------------
float TreeReducerFake::PUWeight(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!sPUWeights)
    return 1.0;
  
  return sPUWeights->GetBinContent(sPUWeights->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
float TreeReducerFake::PUWeightUp(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!sPUWeightsUp)
    return 1.0;
  
  return sPUWeightsUp->GetBinContent(sPUWeightsUp->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
float TreeReducerFake::PUWeightDown(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!sPUWeightsDown)
    return 1.0;
  
  return sPUWeightsDown->GetBinContent(sPUWeightsDown->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
float TreeReducerFake::KFactorWeight(Float_t scale, Float_t phet)
{
  //Z>nunu, taken from https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=242281
  if (scale == -1) {
    if (phet < 145.)
      return 1.0;
    else 
      return (1.71 - 15.86/(phet - 122.93));
  }
  //QCD, taken from dedicated studies
  else if (scale == -2) {
    if (phet < 145.)
      return 1.0;
    else 
      return (0.0068 - 4.79e-06*phet);
  }
  else
    return 1.;
}

//--------------------------------------------------------------------------------------------------
bool TreeReducerFake::EventIsSelected(MitGPTree &tree, int treeType, int& theGoodPhoton, int& photonIsFake)
{
  //standard selection
  if (treeType == 0) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
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
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
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
  //QCD selection
  else if (treeType == 2) {
    if ( tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    photonIsFake  = 0;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5 &&
      tree.phoCombIso1_a1_ <= 18.0 && tree.phoCombIso2_a1_ <= 30.0 && tree.phoCombIso3_a1_ <= 11.4 &&
      tree.phoCoviEtaiEta_a1_ <= 0.014 && tree.phoHadOverEm_a1_ <= 0.124      
      ) { 
        if (PhotonIsSelected(tree.phoR9_a1_, tree.phoHadOverEm_a1_, tree.phoCoviEtaiEta_a1_, tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_))
          photonIsFake  = 1;
        if (PhotonIsFake(tree.phoR9_a1_, tree.phoHadOverEm_a1_, tree.phoCoviEtaiEta_a1_, tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_))
          photonIsFake += 2;
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5 &&
      tree.phoCombIso1_a2_ <= 18.0 && tree.phoCombIso2_a2_ <= 30.0 && tree.phoCombIso3_a2_ <= 11.4 &&
      tree.phoCoviEtaiEta_a2_ <= 0.014 && tree.phoHadOverEm_a2_ <= 0.124      
      ) { 
        if (PhotonIsSelected(tree.phoR9_a2_, tree.phoHadOverEm_a2_, tree.phoCoviEtaiEta_a2_, tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_))
          photonIsFake  = 1;
        if (PhotonIsFake(tree.phoR9_a2_, tree.phoHadOverEm_a2_, tree.phoCoviEtaiEta_a2_, tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_))
          photonIsFake += 2;
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 5. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 1.5 &&
      tree.phoCombIso1_a3_ <= 18.0 && tree.phoCombIso2_a3_ <= 30.0 && tree.phoCombIso3_a3_ <= 11.4 &&
      tree.phoCoviEtaiEta_a3_ <= 0.014 && tree.phoHadOverEm_a3_ <= 0.124      
      ) { 
        if (PhotonIsSelected(tree.phoR9_a3_, tree.phoHadOverEm_a3_, tree.phoCoviEtaiEta_a3_, tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_))
          photonIsFake  = 1;
        if (PhotonIsFake(tree.phoR9_a3_, tree.phoHadOverEm_a3_, tree.phoCoviEtaiEta_a3_, tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_))
          photonIsFake += 2;
        theGoodPhoton = 2;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //Beam halo selection
  else if (treeType == 3) {
    //met and photon number
    if ( tree.metCor_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 5. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 1.5 &&
      PhotonIsSelected(tree.phoR9_a1_, tree.phoHadOverEm_a1_, -1., tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 5. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 1.5 &&
      PhotonIsSelected(tree.phoR9_a2_, tree.phoHadOverEm_a2_, -1., tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
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
bool TreeReducerFake::PhotonIsSelected(float R9, float HoverE, float CoviEtaiEta, float Iso1, float Iso2, float Iso3)
{
  int Cat ;
  Cat = 0;
  if (R9 < 0.94 ) Cat = 1;

  //Isolation cuts

  if (Cat == 0) {
    if (Iso1 > 6.0) return false ;
    if (Iso2 > 10.0) return false ;
    if (Iso3 > 3.8) return false ;
    if (HoverE > 0.124) return false ;
    if (R9 < 0.94 ) return false ; 
    return true;
  }
  else if (Cat == 1) {
    if (Iso1 > 4.7) return false ;
    if (Iso2 > 6.5) return false ;
    if (Iso3 > 2.5) return false ;
    if (HoverE > 0.092) return false ;
    if (R9 < 0.28 ) return false ; 
    return true;
  }
  else
    return false;
}

//--------------------------------------------------------------------------------------------------
bool TreeReducerFake::PhotonIsFake(float R9, float HoverE, float CoviEtaiEta, float Iso1, float Iso2, float Iso3)
{
  int Cat ;
  Cat = 0;
  if (R9 < 0.94 ) Cat = 1;

  //Isolation cuts

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
float TreeReducerFake::GetCorrDeltaPhi(float phi1, float phi2)
{
  float corrDeltaPhi = TMath::Abs(phi1 - phi2);
  if (corrDeltaPhi > TMath::Pi())
    corrDeltaPhi = TMath::TwoPi() - corrDeltaPhi;     
  return corrDeltaPhi;
}
