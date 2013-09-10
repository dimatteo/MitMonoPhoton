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
#include "MitMonoPhoton/Utils/interface/TreeReducer.h"
#include "MitMonoPhoton/Core/MitGPTree.h"
#include "MitMonoPhoton/Core/MitGPTreeReduced.h"

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
  bool isMC = false, isData = false;
  int treeType = 0;
  float thisXsec   = *fSample->Xsec();
  float thisScale  = *fSample->Scale();
  if (thisXsec > 0)         isMC = true;
  else if (thisXsec == -1)  isData = true;
  else if (thisXsec == -2) {treeType = 2;} //QCD
  else                     {treeType = 3;} //BeamHalo

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
  if (treeType == 2) outTreeName = "QCD";
  if (treeType == 3) outTreeName = "BeamHalo";
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
    if ( !EventIsSelected(intree, treeType, theGoodPhoton) ) continue;

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
    outtree.met_   = intree.met_;
    outtree.metPhi_= intree.metPhi_;
    outtree.sumEt_ = intree.sumEt_;
    outtree.metSig_= intree.metSig_;
  
    outtree.nphotons_ = intree.nphotons_;
    if ( theGoodPhoton == 0 ) {
      outtree.phoEt_ = intree.pho1_.Et();
      outtree.phoEta_ = intree.pho1_.Eta();
      outtree.phoPhi_ = intree.pho1_.Phi();
    }
    else if ( theGoodPhoton == 1 ) {
      outtree.phoEt_ = intree.pho2_.Et();
      outtree.phoEta_ = intree.pho2_.Eta();
      outtree.phoPhi_ = intree.pho2_.Phi();
    }
    else if ( theGoodPhoton == 2 ) {
      outtree.phoEt_ = intree.pho3_.Et();
      outtree.phoEta_ = intree.pho3_.Eta();
      outtree.phoPhi_ = intree.pho3_.Phi();
    }
    else if ( theGoodPhoton == 3 ) {
      outtree.phoEt_ = intree.pho4_.Et();
      outtree.phoEta_ = intree.pho4_.Eta();
      outtree.phoPhi_ = intree.pho4_.Phi();
    }
    else
      cout << "Error in the selection function! theGoodPhoton==" << theGoodPhoton  << endl; 
  
    outtree.phoMetDeltaPhi_ = GetCorrDeltaPhi(outtree.phoPhi_, intree.metPhi_);
    
    std::vector<float> jet_Pt;
    std::vector<float> jet_Eta;
    std::vector<float> jet_Phi;
    outtree.nalljets_ = intree.njets_;
    for (unsigned int ijet = 0; ijet < intree.njets_; ijet++) {
      if (ijet == 0) {
        outtree.jet1Pt_ = intree.jet1_.Pt();
        outtree.jet1Eta_ = intree.jet1_.Eta();
        if (fabs(intree.jet1_.Eta()) > 2.4) continue;
        jet_Pt.push_back(intree.jet1_.Pt());
        jet_Eta.push_back(intree.jet1_.Eta());
        jet_Phi.push_back(intree.jet1_.Phi());
      }
      else if (ijet == 1) {
        outtree.jet2Pt_ = intree.jet2_.Pt();
        outtree.jet2Eta_ = intree.jet2_.Eta();
        if (fabs(intree.jet2_.Eta()) > 2.4) continue;
        jet_Pt.push_back(intree.jet2_.Pt());
        jet_Eta.push_back(intree.jet2_.Eta());
        jet_Phi.push_back(intree.jet2_.Phi());
      }
      else if (ijet == 2) {
        outtree.jet3Pt_ = intree.jet3_.Pt();
        outtree.jet3Eta_ = intree.jet3_.Eta();
        if (fabs(intree.jet3_.Eta()) > 2.4) continue;
        jet_Pt.push_back(intree.jet3_.Pt());
        jet_Eta.push_back(intree.jet3_.Eta());
        jet_Phi.push_back(intree.jet3_.Phi());
      }
      else if (ijet == 3) {
        outtree.jet4Pt_ = intree.jet4_.Pt();
        outtree.jet4Eta_ = intree.jet4_.Eta();
        if (fabs(intree.jet4_.Eta()) > 2.4) continue;
        jet_Pt.push_back(intree.jet4_.Pt());
        jet_Eta.push_back(intree.jet4_.Eta());
        jet_Phi.push_back(intree.jet4_.Phi());
      }
      else break;
    }
    outtree.njets_ = jet_Pt.size();
    if (jet_Pt.size() > 0) {
      outtree.leadJetPt_ = jet_Pt[0];
      outtree.leadJetEta_= jet_Eta[0];
      outtree.leadJetPhi_= jet_Phi[0];
    }
    if (jet_Pt.size() > 1) {
      outtree.trailJetPt_ = jet_Pt[1];
      outtree.trailJetEta_= jet_Eta[1];
      outtree.trailJetPhi_= jet_Phi[1];
    }
  
    //leptons
    outtree.nlep_ = intree.nlep_;
    for (unsigned int ilep = 0; ilep < intree.nlep_; ilep++) {
      if (ilep == 0) {
        outtree.lep1Pt_ = intree.lep1_.Pt();
        outtree.lep1Eta_ = intree.lep1_.Eta();
        outtree.lep1Id_ = intree.lid1_;
      }
      else if (ilep == 1) {
        outtree.lep2Pt_ = intree.lep2_.Pt();
        outtree.lep2Eta_ = intree.lep2_.Eta();
        outtree.lep2Id_ = intree.lid2_;
      }
      else if (ilep == 3) {
        outtree.lep3Pt_ = intree.lep3_.Pt();
        outtree.lep3Eta_ = intree.lep3_.Eta();
        outtree.lep3Id_ = intree.lid3_;
      }
      else break;
    }

    //cosmics
    outtree.ncosmics_ = intree.ncosmics_;

    //fill the outtree content
    outtree.tree_->Fill();    
  }
    
  fOutFile   -> cd();
  fOutFile   -> Write();  
  return;
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
bool TreeReducer::EventIsSelected(MitGPTree &tree, int treeType, int& theGoodPhoton)
{
  //standard selection
  if (treeType == 0) {
    //met and photon number
    if ( tree.met_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 8. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 3.0) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 8. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 3.0) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 8. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 3.0) { 
        theGoodPhoton = 2;
        break;
      }
      else if ( ipho == 3 && abs(tree.pho4_.Eta()) < 1.479 && tree.pho4_.Et() > 140 && 
      tree.phoPassEleVeto_a4_ == 1 && tree.phoIsTrigger_a4_ == 1 && 
      abs(tree.phoLeadTimeSpan_a4_) < 8. && tree.phoCoviEtaiEta_a4_ > 0.001 && tree.phoCoviPhiiPhi_a4_ > 0.001 && 
      tree.phoMipIsHalo_a4_ == 0 && tree.phoSeedTime_a4_ > -1.5 && tree.phoSeedTime_a4_ < 3.0) { 
        theGoodPhoton = 3;
        break;
      }
      else break;
    }
    if ( theGoodPhoton < 0 ) return false;
    return true;
  }
  //QCD selection
  else if (treeType == 2) {
    //met and photon number
    if ( tree.met_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 8. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a1_, tree.phoHadOverEm_a1_, tree.phoCoviEtaiEta_a1_, tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 8. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a2_, tree.phoHadOverEm_a2_, tree.phoCoviEtaiEta_a2_, tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 8. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a3_, tree.phoHadOverEm_a3_, tree.phoCoviEtaiEta_a3_, tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else if ( ipho == 3 && abs(tree.pho4_.Eta()) < 1.479 && tree.pho4_.Et() > 140 && 
      tree.phoPassEleVeto_a4_ == 1 && tree.phoIsTrigger_a4_ == 1 && 
      abs(tree.phoLeadTimeSpan_a4_) < 8. && tree.phoCoviEtaiEta_a4_ > 0.001 && tree.phoCoviPhiiPhi_a4_ > 0.001 && 
      tree.phoMipIsHalo_a4_ == 0 && tree.phoSeedTime_a4_ > -1.5 && tree.phoSeedTime_a4_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a4_, tree.phoHadOverEm_a4_, tree.phoCoviEtaiEta_a4_, tree.phoCombIso1_a4_, tree.phoCombIso2_a4_, tree.phoCombIso3_a4_)) { 
        theGoodPhoton = 3;
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
    if ( tree.met_ < 100 || tree.nphotons_ == 0 ) return false;
    //Phase space & photons
    theGoodPhoton = -1;
    for ( unsigned int ipho = 0; ipho < tree.nphotons_; ipho++ ) {
      if ( ipho == 0 && abs(tree.pho1_.Eta()) < 1.479 && tree.pho1_.Et() > 140 && 
      tree.phoPassEleVeto_a1_ == 1 && tree.phoIsTrigger_a1_ == 1 && 
      abs(tree.phoLeadTimeSpan_a1_) < 8. && tree.phoCoviEtaiEta_a1_ > 0.001 && tree.phoCoviPhiiPhi_a1_ > 0.001 && 
      tree.phoMipIsHalo_a1_ == 0 && tree.phoSeedTime_a1_ > -1.5 && tree.phoSeedTime_a1_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a1_, tree.phoHadOverEm_a1_, -1., tree.phoCombIso1_a1_, tree.phoCombIso2_a1_, tree.phoCombIso3_a1_)) { 
        theGoodPhoton = 0;
        break;
      }
      else if ( ipho == 1 && abs(tree.pho2_.Eta()) < 1.479 && tree.pho2_.Et() > 140 && 
      tree.phoPassEleVeto_a2_ == 1 && tree.phoIsTrigger_a2_ == 1 && 
      abs(tree.phoLeadTimeSpan_a2_) < 8. && tree.phoCoviEtaiEta_a2_ > 0.001 && tree.phoCoviPhiiPhi_a2_ > 0.001 && 
      tree.phoMipIsHalo_a2_ == 0 && tree.phoSeedTime_a2_ > -1.5 && tree.phoSeedTime_a2_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a2_, tree.phoHadOverEm_a2_, -1., tree.phoCombIso1_a2_, tree.phoCombIso2_a2_, tree.phoCombIso3_a2_)) { 
        theGoodPhoton = 1;
        break;
      }
      else if ( ipho == 2 && abs(tree.pho3_.Eta()) < 1.479 && tree.pho3_.Et() > 140 && 
      tree.phoPassEleVeto_a3_ == 1 && tree.phoIsTrigger_a3_ == 1 && 
      abs(tree.phoLeadTimeSpan_a3_) < 8. && tree.phoCoviEtaiEta_a3_ > 0.001 && tree.phoCoviPhiiPhi_a3_ > 0.001 && 
      tree.phoMipIsHalo_a3_ == 0 && tree.phoSeedTime_a3_ > -1.5 && tree.phoSeedTime_a3_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a3_, tree.phoHadOverEm_a3_, -1., tree.phoCombIso1_a3_, tree.phoCombIso2_a3_, tree.phoCombIso3_a3_)) { 
        theGoodPhoton = 2;
        break;
      }
      else if ( ipho == 3 && abs(tree.pho4_.Eta()) < 1.479 && tree.pho4_.Et() > 140 && 
      tree.phoPassEleVeto_a4_ == 1 && tree.phoIsTrigger_a4_ == 1 && 
      abs(tree.phoLeadTimeSpan_a4_) < 8. && tree.phoCoviEtaiEta_a4_ > 0.001 && tree.phoCoviPhiiPhi_a4_ > 0.001 && 
      tree.phoMipIsHalo_a4_ == 0 && tree.phoSeedTime_a4_ > -1.5 && tree.phoSeedTime_a4_ < 3.0 &&
      PhotonIsSelected(tree.phoR9_a4_, tree.phoHadOverEm_a4_, -1., tree.phoCombIso1_a4_, tree.phoCombIso2_a4_, tree.phoCombIso3_a4_)) { 
        theGoodPhoton = 3;
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
float TreeReducer::GetCorrDeltaPhi(float phi1, float phi2)
{
  float corrDeltaPhi = TMath::Abs(phi1 - phi2);
  if (corrDeltaPhi > TMath::Pi())
    corrDeltaPhi = TMath::TwoPi() - corrDeltaPhi;     
  return corrDeltaPhi;
}
