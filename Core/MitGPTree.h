#ifndef MitGPTree_H
#define MitGPTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitGPTree {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

  /// DON'T CHANGE ORDER
  enum DataType {
    data,
    dmg,
    dmj,
    top,
    dyll,
    wjets,
    vv,
    wgamma,
    vvv,
    gj,
    qcd,
    other
  };

  enum Selection {
    GoodPhoton      = 1UL<<0,    // good photon
    Fake            = 1UL<<1,    // fake
    BeamHalo        = 1UL<<2,    // beam halo
    DiLepton        = 1UL<<3     // dilepton
  };

  /// variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   nvtx_;
  unsigned int   cuts_;
  float          scale1fb_;
  float          met_;
  float          metPhi_;
  float          sumEt_;
  float          metSig_;
  DataType       dstype_;
  float          metCor_;
  float          metCorPhi_;

  unsigned int   nlep_;
  LorentzVector  lep1_;
  int            lid1_;
  float          etaOfSTA_l1_;
  float          phiOfSTA_l1_;
  LorentzVector  lep2_;
  int            lid2_;
  float          etaOfSTA_l2_;
  float          phiOfSTA_l2_;
  LorentzVector  lep3_;
  int            lid3_;
  float          etaOfSTA_l3_;
  float          phiOfSTA_l3_;

  unsigned int   nphotons_;
  LorentzVector  pho1_;
  float phoCombIso1_a1_;
  float phoCombIso2_a1_;
  float phoCombIso3_a1_;
  float phoWorstIso_a1_;
  bool  phoPassEleVeto_a1_;
  bool  phoHasPixelSeed_a1_;
  float phoCoviEtaiEta_a1_;
  float phoCoviPhiiPhi_a1_;
  float phoR9_a1_;
  float phoSeedTime_a1_;
  float phoHadOverEm_a1_;
  float phoLeadTimeSpan_a1_;
  float phoSubLeadTimeSpan_a1_;
  float phoMipChi2_a1_;
  float phoMipTotEnergy_a1_;
  float phoMipSlope_a1_;
  float phoMipIntercept_a1_;
  int   phoMipNhitCone_a1_;
  bool  phoMipIsHalo_a1_;
  int   phoMatchType_a1_;
  float phoMatchPt_a1_;
  bool  phoIsTrigger_a1_;
  LorentzVector  pho2_;
  float phoCombIso1_a2_;
  float phoCombIso2_a2_;
  float phoCombIso3_a2_;
  float phoWorstIso_a2_;
  bool  phoPassEleVeto_a2_;
  bool  phoHasPixelSeed_a2_;
  float phoCoviEtaiEta_a2_;
  float phoCoviPhiiPhi_a2_;
  float phoR9_a2_;
  float phoSeedTime_a2_;
  float phoHadOverEm_a2_;
  float phoLeadTimeSpan_a2_;
  float phoSubLeadTimeSpan_a2_;
  float phoMipChi2_a2_;
  float phoMipTotEnergy_a2_;
  float phoMipSlope_a2_;
  float phoMipIntercept_a2_;
  int   phoMipNhitCone_a2_;
  bool  phoMipIsHalo_a2_;
  int   phoMatchType_a2_;
  float phoMatchPt_a2_;
  bool  phoIsTrigger_a2_;
  LorentzVector  pho3_;
  float phoCombIso1_a3_;
  float phoCombIso2_a3_;
  float phoCombIso3_a3_;
  float phoWorstIso_a3_;
  bool  phoPassEleVeto_a3_;
  bool  phoHasPixelSeed_a3_;
  float phoCoviEtaiEta_a3_;
  float phoCoviPhiiPhi_a3_;
  float phoR9_a3_;
  float phoSeedTime_a3_;
  float phoHadOverEm_a3_;
  float phoLeadTimeSpan_a3_;
  float phoSubLeadTimeSpan_a3_;
  float phoMipChi2_a3_;
  float phoMipTotEnergy_a3_;
  float phoMipSlope_a3_;
  float phoMipIntercept_a3_;
  int   phoMipNhitCone_a3_;
  bool  phoMipIsHalo_a3_;
  int   phoMatchType_a3_;
  float phoMatchPt_a3_;
  bool  phoIsTrigger_a3_;

  unsigned int   njets_;
  LorentzVector  jet1_;
  float          jet1Btag_;
  LorentzVector  jet2_;
  float          jet2Btag_;
  LorentzVector  jet3_;
  float          jet3Btag_;
  LorentzVector  jet4_;
  float          jet4Btag_;

  unsigned int   ntracks_;
  LorentzVector  track1_;
  LorentzVector  track2_;
  LorentzVector  track3_;

  unsigned int   ncosmics_;
  LorentzVector  cosmic1_;
  float          etaOfSTA_c1_;
  float          phiOfSTA_c1_;
  LorentzVector  cosmic2_;
  float          etaOfSTA_c2_;
  float          phiOfSTA_c2_;
  LorentzVector  cosmic3_;
  float          etaOfSTA_c3_;
  float          phiOfSTA_c3_;

  float          Q_;
  float          id1_;
  float          x1_;
  float          pdf1_;
  float          id2_;
  float          x2_;
  float          pdf2_;
  int            processId_;
  float          npu_;
  float          npuPlusOne_;
  float          npuMinusOne_;
  float          rho_;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitGPTree():
    lepPtr1_(&lep1_),lepPtr2_(&lep2_),lepPtr3_(&lep3_),
    phoPtr1_(&pho1_),phoPtr2_(&pho2_),phoPtr3_(&pho3_),
    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),
    trackPtr1_(&track1_),trackPtr2_(&track2_),trackPtr3_(&track3_),
    cosmicPtr1_(&cosmic1_),cosmicPtr2_(&cosmic2_),cosmicPtr3_(&cosmic3_){}
  /// default destructor
  ~MitGPTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitGPTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job DMTree0/1
    // type == 0/1/2/3 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    f_ = TFile::Open(file);
    assert(f_);
    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("MPhoTree"));
    else if(type == 1) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("MPhoTreeLepton"));
    else if(type == 2) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("MPhoTreePhFake"));
    else if(type == 3) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("MPhoTreeBeamHalo"));
    else               tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("tree"));
    assert(tree_);
  }

  /// create a MitGPTree
  void CreateTree(int type = -1){
    assert(type==type); // just to suppress warnings
    // to create three different ntuples in the same job DMTree0/1
    // type == 0/1/2/3 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    if     (type == 0) tree_ = new TTree("MPhoTree","Smurf ntuples");
    else if(type == 1) tree_ = new TTree("MPhoTreeLepton","Smurf ntuples");
    else if(type == 2) tree_ = new TTree("MPhoTreePhFake","Smurf ntuples");
    else if(type == 3) tree_ = new TTree("MPhoTreeBeamHalo","Smurf ntuples");
    else               tree_ = new TTree("tree","Smurf ntuples");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("event"        , &event_        ,   "event/i");
    tree_->Branch("run"          , &run_          ,   "run/i");
    tree_->Branch("lumi"         , &lumi_         ,   "lumi/i");
    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("scale1fb"     , &scale1fb_     ,   "scale1fb/F");
    tree_->Branch("met"          , &met_          ,   "met/F");
    tree_->Branch("metPhi"       , &metPhi_       ,   "metPhi/F");
    tree_->Branch("sumEt"        , &sumEt_        ,   "sumEt/F");
    tree_->Branch("metSig"       , &metSig_       ,   "metSig/F");
    tree_->Branch("dstype"       , &dstype_       ,   "dstype/I");
    tree_->Branch("metCor"          , &metCor_          ,   "metCor/F");
    tree_->Branch("metCorPhi"       , &metCorPhi_       ,   "metCorPhi/F");

    tree_->Branch("nlep"         , &nlep_         ,   "nlep/i");
    tree_->Branch("lep1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lid1"         , &lid1_         ,   "lid1/I");
    tree_->Branch("etaOfSTA_l1"  , &etaOfSTA_l1_  ,   "etaOfSTA_l1/F");
    tree_->Branch("phiOfSTA_l1"  , &phiOfSTA_l1_  ,   "phiOfSTA_l1/F");
    tree_->Branch("lep2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lid2"         , &lid2_         ,   "lid2/I");
    tree_->Branch("etaOfSTA_l2"  , &etaOfSTA_l2_  ,   "etaOfSTA_l2/F");
    tree_->Branch("phiOfSTA_l2"  , &phiOfSTA_l2_  ,   "phiOfSTA_l2/F");
    tree_->Branch("lep3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
    tree_->Branch("lid3"         , &lid3_         ,   "lid3/I");
    tree_->Branch("etaOfSTA_l3"  , &etaOfSTA_l3_  ,   "etaOfSTA_l3/F");
    tree_->Branch("phiOfSTA_l3"  , &phiOfSTA_l3_  ,   "phiOfSTA_l3/F");

    tree_->Branch("nphotons"     , &nphotons_     ,   "nphotons/i");
    tree_->Branch("pho1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr1_);
    tree_->Branch("phoCombIso1_a1"            , &phoCombIso1_a1_            , "phoCombIso1_a1/F");
    tree_->Branch("phoCombIso2_a1"            , &phoCombIso2_a1_            , "phoCombIso2_a1/F");
    tree_->Branch("phoCombIso3_a1"            , &phoCombIso3_a1_            , "phoCombIso3_a1/F");
    tree_->Branch("phoWorstIso_a1"            , &phoWorstIso_a1_            , "phoWorstIso_a1/F");
    tree_->Branch("phoPassEleVeto_a1"         , &phoPassEleVeto_a1_         , "phoPassEleVeto_a1/B");
    tree_->Branch("phoHasPixelSeed_a1"        , &phoHasPixelSeed_a1_        , "phoHasPixelSeed_a1/B");
    tree_->Branch("phoCoviEtaiEta_a1"         , &phoCoviEtaiEta_a1_         , "phoCoviEtaiEta_a1/F");
    tree_->Branch("phoCoviPhiiPhi_a1"         , &phoCoviPhiiPhi_a1_         , "phoCoviPhiiPhi_a1/F");
    tree_->Branch("phoR9_a1"                  , &phoR9_a1_                  , "phoR9_a1/F");
    tree_->Branch("phoSeedTime_a1"            , &phoSeedTime_a1_            , "phoSeedTime_a1/F");
    tree_->Branch("phoHadOverEm_a1"           , &phoHadOverEm_a1_           , "phoHadOverEm_a1/F");
    tree_->Branch("phoLeadTimeSpan_a1"        , &phoLeadTimeSpan_a1_        , "phoLeadTimeSpan_a1/F");
    tree_->Branch("phoSubLeadTimeSpan_a1"     , &phoSubLeadTimeSpan_a1_     , "phoSubLeadTimeSpan_a1/F");
    tree_->Branch("phoMipChi2_a1"             , &phoMipChi2_a1_             , "phoMipChi2_a1/F");
    tree_->Branch("phoMipTotEnergy_a1"        , &phoMipTotEnergy_a1_        , "phoMipTotEnergy_a1/F");
    tree_->Branch("phoMipSlope_a1"            , &phoMipSlope_a1_            , "phoMipSlope_a1/F");
    tree_->Branch("phoMipIntercept_a1"        , &phoMipIntercept_a1_        , "phoMipIntercept_a1/F");
    tree_->Branch("phoMipNhitCone_a1"         , &phoMipNhitCone_a1_         , "phoMipNhitCone_a1/I");
    tree_->Branch("phoMipIsHalo_a1"           , &phoMipIsHalo_a1_           , "phoMipIsHalo_a1/B");
    tree_->Branch("phoMatchType_a1"           , &phoMatchType_a1_           , "phoMatchType_a1/i");
    tree_->Branch("phoMatchPt_a1"             , &phoMatchPt_a1_             , "phoMatchPt_a1/F");
    tree_->Branch("phoIsTrigger_a1"           , &phoIsTrigger_a1_           , "phoIsTrigger_a1_/B");
    tree_->Branch("pho2"                      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr2_);
    tree_->Branch("phoCombIso1_a2"            , &phoCombIso1_a2_            , "phoCombIso1_a2/F");
    tree_->Branch("phoCombIso2_a2"            , &phoCombIso2_a2_            , "phoCombIso2_a2/F");
    tree_->Branch("phoCombIso3_a2"            , &phoCombIso3_a2_            , "phoCombIso3_a2/F");
    tree_->Branch("phoWorstIso_a2"            , &phoWorstIso_a2_            , "phoWorstIso_a2/F");
    tree_->Branch("phoPassEleVeto_a2"         , &phoPassEleVeto_a2_         , "phoPassEleVeto_a2/B");
    tree_->Branch("phoHasPixelSeed_a2"        , &phoHasPixelSeed_a2_        , "phoHasPixelSeed_a2/B");
    tree_->Branch("phoCoviEtaiEta_a2"         , &phoCoviEtaiEta_a2_         , "phoCoviEtaiEta_a2/F");
    tree_->Branch("phoCoviPhiiPhi_a2"         , &phoCoviPhiiPhi_a2_         , "phoCoviPhiiPhi_a2/F");
    tree_->Branch("phoR9_a2"                  , &phoR9_a2_                  , "phoR9_a2/F");
    tree_->Branch("phoSeedTime_a2"            , &phoSeedTime_a2_            , "phoSeedTime_a2/F");
    tree_->Branch("phoHadOverEm_a2"           , &phoHadOverEm_a2_           , "phoHadOverEm_a2/F");
    tree_->Branch("phoLeadTimeSpan_a2"        , &phoLeadTimeSpan_a2_        , "phoLeadTimeSpan_a2/F");
    tree_->Branch("phoSubLeadTimeSpan_a2"     , &phoSubLeadTimeSpan_a2_     , "phoSubLeadTimeSpan_a2/F");
    tree_->Branch("phoMipChi2_a2"             , &phoMipChi2_a2_             , "phoMipChi2_a2/F");
    tree_->Branch("phoMipTotEnergy_a2"        , &phoMipTotEnergy_a2_        , "phoMipTotEnergy_a2/F");
    tree_->Branch("phoMipSlope_a2"            , &phoMipSlope_a2_            , "phoMipSlope_a2/F");
    tree_->Branch("phoMipIntercept_a2"        , &phoMipIntercept_a2_        , "phoMipIntercept_a2/F");
    tree_->Branch("phoMipNhitCone_a2"         , &phoMipNhitCone_a2_         , "phoMipNhitCone_a2/I");
    tree_->Branch("phoMipIsHalo_a2"           , &phoMipIsHalo_a2_           , "phoMipIsHalo_a2/B");
    tree_->Branch("phoMatchType_a2"           , &phoMatchType_a2_           , "phoMatchType_a2/i");
    tree_->Branch("phoMatchPt_a2"             , &phoMatchPt_a2_             , "phoMatchPt_a2/F");
    tree_->Branch("phoIsTrigger_a2"           , &phoIsTrigger_a2_           , "phoIsTrigger_a2_/B");
    tree_->Branch("pho3"                      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr3_);
    tree_->Branch("phoCombIso1_a3"            , &phoCombIso1_a3_            , "phoCombIso1_a3/F");
    tree_->Branch("phoCombIso2_a3"            , &phoCombIso2_a3_            , "phoCombIso2_a3/F");
    tree_->Branch("phoCombIso3_a3"            , &phoCombIso3_a3_            , "phoCombIso3_a3/F");
    tree_->Branch("phoWorstIso_a3"            , &phoWorstIso_a3_            , "phoWorstIso_a3/F");
    tree_->Branch("phoPassEleVeto_a3"         , &phoPassEleVeto_a3_         , "phoPassEleVeto_a3/B");
    tree_->Branch("phoHasPixelSeed_a3"        , &phoHasPixelSeed_a3_        , "phoHasPixelSeed_a3/B");
    tree_->Branch("phoCoviEtaiEta_a3"         , &phoCoviEtaiEta_a3_         , "phoCoviEtaiEta_a3/F");
    tree_->Branch("phoCoviPhiiPhi_a3"         , &phoCoviPhiiPhi_a3_         , "phoCoviPhiiPhi_a3/F");
    tree_->Branch("phoR9_a3"                  , &phoR9_a3_                  , "phoR9_a3/F");
    tree_->Branch("phoSeedTime_a3"            , &phoSeedTime_a3_            , "phoSeedTime_a3/F");
    tree_->Branch("phoHadOverEm_a3"           , &phoHadOverEm_a3_           , "phoHadOverEm_a3/F");
    tree_->Branch("phoLeadTimeSpan_a3"        , &phoLeadTimeSpan_a3_        , "phoLeadTimeSpan_a3/F");
    tree_->Branch("phoSubLeadTimeSpan_a3"     , &phoSubLeadTimeSpan_a3_     , "phoSubLeadTimeSpan_a3/F");
    tree_->Branch("phoMipChi2_a3"             , &phoMipChi2_a3_             , "phoMipChi2_a3/F");
    tree_->Branch("phoMipTotEnergy_a3"        , &phoMipTotEnergy_a3_        , "phoMipTotEnergy_a3/F");
    tree_->Branch("phoMipSlope_a3"            , &phoMipSlope_a3_            , "phoMipSlope_a3/F");
    tree_->Branch("phoMipIntercept_a3"        , &phoMipIntercept_a3_        , "phoMipIntercept_a3/F");
    tree_->Branch("phoMipNhitCone_a3"         , &phoMipNhitCone_a3_         , "phoMipNhitCone_a3/I");
    tree_->Branch("phoMipIsHalo_a3"           , &phoMipIsHalo_a3_           , "phoMipIsHalo_a3/B");
    tree_->Branch("phoMatchType_a3"           , &phoMatchType_a3_          , "phoMatchType_a3/i");
    tree_->Branch("phoMatchPt_a3"             , &phoMatchPt_a3_          , "phoMatchPt_a3/F");
    tree_->Branch("phoIsTrigger_a3"           , &phoIsTrigger_a3_           , "phoIsTrigger_a3_/B");

    tree_->Branch("njets"        , &njets_        ,   "njets/i");
    tree_->Branch("jet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1Btag"     , &jet1Btag_     ,   "jet1Btag/F");
    tree_->Branch("jet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet2Btag"     , &jet2Btag_     ,   "jet2Btag/F");
    tree_->Branch("jet3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
    tree_->Branch("jet3Btag"     , &jet3Btag_     ,   "jet3Btag/F");
    tree_->Branch("jet4"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr4_);
    tree_->Branch("jet4Btag"     , &jet4Btag_     ,   "jet4Btag/F");

    tree_->Branch("ntracks"        , &ntracks_        ,   "ntracks/i");
    tree_->Branch("track1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &trackPtr1_);
    tree_->Branch("track2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &trackPtr2_);
    tree_->Branch("track3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &trackPtr3_);

    tree_->Branch("ncosmics"        , &ncosmics_        ,   "ncosmics/i");
    tree_->Branch("cosmic1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &cosmicPtr1_);
    tree_->Branch("etaOfSTA_c1"     , &etaOfSTA_c1_     ,   "etaOfSTA_c1/F");
    tree_->Branch("phiOfSTA_c1"     , &phiOfSTA_c1_     ,   "phiOfSTA_c1/F");
    tree_->Branch("cosmic2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &cosmicPtr2_);
    tree_->Branch("etaOfSTA_c2"     , &etaOfSTA_c2_     ,   "etaOfSTA_c2/F");
    tree_->Branch("phiOfSTA_c2"     , &phiOfSTA_c2_     ,   "phiOfSTA_c2/F");
    tree_->Branch("cosmic3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &cosmicPtr3_);
    tree_->Branch("etaOfSTA_c3"     , &etaOfSTA_c3_     ,   "etaOfSTA_c3/F");
    tree_->Branch("phiOfSTA_c3"     , &phiOfSTA_c3_     ,   "phiOfSTA_c3/F");

    tree_->Branch("Q",             &Q_	  ,     "Q/F");
    tree_->Branch("id1",           &id1_  ,     "id1/F");
    tree_->Branch("x1",            &x1_	  ,     "x1/F");
    tree_->Branch("pdf1",          &pdf1_ ,     "pdf1/F");
    tree_->Branch("id2",           &id2_  ,     "id2/F");
    tree_->Branch("x2",            &x2_	  ,     "x2/F");
    tree_->Branch("pdf2",          &pdf2_ ,     "pdf2/F");
    tree_->Branch("processId",     &processId_  , "processId/I");
    tree_->Branch("npu",           &npu_        , "npu/F");
    tree_->Branch("npuPlusOne",    &npuPlusOne_ , "npuPlusOne/F");
    tree_->Branch("npuMinusOne",   &npuMinusOne_, "npuMinusOne/F");
    tree_->Branch("rho",           &rho_,         "rho/F");

  }

  // initialze a MitGPTree
  void InitTree(int type = -1){
    assert(type==type); // just to suppress warnings
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("scale1fb",      &scale1fb_);
    tree_->SetBranchAddress("met",           &met_);
    tree_->SetBranchAddress("metPhi",        &metPhi_);
    tree_->SetBranchAddress("sumEt",         &sumEt_);
    tree_->SetBranchAddress("metSig",        &metSig_);
    tree_->SetBranchAddress("dstype",        &dstype_);
    tree_->SetBranchAddress("metCor"          , &metCor_   );
    tree_->SetBranchAddress("metCorPhi"       , &metCorPhi_);

    tree_->SetBranchAddress("nlep",          &nlep_);
    tree_->SetBranchAddress("lep1",          &lepPtr1_);
    tree_->SetBranchAddress("lid1",          &lid1_);
    tree_->SetBranchAddress("etaOfSTA_l1",   &etaOfSTA_l1_);
    tree_->SetBranchAddress("phiOfSTA_l1",   &phiOfSTA_l1_);
    tree_->SetBranchAddress("lep2",          &lepPtr2_);
    tree_->SetBranchAddress("lid2",          &lid2_);
    tree_->SetBranchAddress("etaOfSTA_l2",   &etaOfSTA_l2_);
    tree_->SetBranchAddress("phiOfSTA_l2",   &phiOfSTA_l2_);
    tree_->SetBranchAddress("lep3",          &lepPtr3_);
    tree_->SetBranchAddress("lid3",          &lid3_);
    tree_->SetBranchAddress("etaOfSTA_l3",   &etaOfSTA_l3_);
    tree_->SetBranchAddress("phiOfSTA_l3",   &phiOfSTA_l3_);

    tree_->SetBranchAddress("nphotons"                  , &nphotons_);
    tree_->SetBranchAddress("pho1"                      , &phoPtr1_);
    tree_->SetBranchAddress("phoCombIso1_a1"            , &phoCombIso1_a1_);
    tree_->SetBranchAddress("phoCombIso2_a1"            , &phoCombIso2_a1_);
    tree_->SetBranchAddress("phoCombIso3_a1"            , &phoCombIso3_a1_);
    tree_->SetBranchAddress("phoWorstIso_a1"            , &phoWorstIso_a1_);
    tree_->SetBranchAddress("phoPassEleVeto_a1"         , &phoPassEleVeto_a1_);
    tree_->SetBranchAddress("phoHasPixelSeed_a1"        , &phoHasPixelSeed_a1_);
    tree_->SetBranchAddress("phoCoviEtaiEta_a1"         , &phoCoviEtaiEta_a1_);
    tree_->SetBranchAddress("phoCoviPhiiPhi_a1"         , &phoCoviPhiiPhi_a1_);
    tree_->SetBranchAddress("phoR9_a1"                  , &phoR9_a1_);
    tree_->SetBranchAddress("phoSeedTime_a1"            , &phoSeedTime_a1_);
    tree_->SetBranchAddress("phoHadOverEm_a1"           , &phoHadOverEm_a1_);
    tree_->SetBranchAddress("phoLeadTimeSpan_a1"        , &phoLeadTimeSpan_a1_);
    tree_->SetBranchAddress("phoSubLeadTimeSpan_a1"     , &phoSubLeadTimeSpan_a1_);
    tree_->SetBranchAddress("phoMipChi2_a1"             , &phoMipChi2_a1_);
    tree_->SetBranchAddress("phoMipTotEnergy_a1"        , &phoMipTotEnergy_a1_);
    tree_->SetBranchAddress("phoMipSlope_a1"            , &phoMipSlope_a1_);
    tree_->SetBranchAddress("phoMipIntercept_a1"        , &phoMipIntercept_a1_);
    tree_->SetBranchAddress("phoMipNhitCone_a1"         , &phoMipNhitCone_a1_);
    tree_->SetBranchAddress("phoMipIsHalo_a1"           , &phoMipIsHalo_a1_);
    tree_->SetBranchAddress("phoMatchType_a1"           , &phoMatchType_a1_);
    tree_->SetBranchAddress("phoMatchPt_a1"             , &phoMatchPt_a1_);
    tree_->SetBranchAddress("phoIsTrigger_a1"           , &phoIsTrigger_a1_);
    tree_->SetBranchAddress("pho2"                      , &phoPtr2_);
    tree_->SetBranchAddress("phoCombIso1_a2"            , &phoCombIso1_a2_);
    tree_->SetBranchAddress("phoCombIso2_a2"            , &phoCombIso2_a2_);
    tree_->SetBranchAddress("phoCombIso3_a2"            , &phoCombIso3_a2_ );
    tree_->SetBranchAddress("phoWorstIso_a2"            , &phoWorstIso_a2_);
    tree_->SetBranchAddress("phoPassEleVeto_a2"         , &phoPassEleVeto_a2_);
    tree_->SetBranchAddress("phoHasPixelSeed_a2"        , &phoHasPixelSeed_a2_);
    tree_->SetBranchAddress("phoCoviEtaiEta_a2"         , &phoCoviEtaiEta_a2_);
    tree_->SetBranchAddress("phoCoviPhiiPhi_a2"         , &phoCoviPhiiPhi_a2_);
    tree_->SetBranchAddress("phoR9_a2"                  , &phoR9_a2_);
    tree_->SetBranchAddress("phoSeedTime_a2"            , &phoSeedTime_a2_);
    tree_->SetBranchAddress("phoHadOverEm_a2"           , &phoHadOverEm_a2_);
    tree_->SetBranchAddress("phoLeadTimeSpan_a2"        , &phoLeadTimeSpan_a2_);
    tree_->SetBranchAddress("phoSubLeadTimeSpan_a2"     , &phoSubLeadTimeSpan_a2_);
    tree_->SetBranchAddress("phoMipChi2_a2"             , &phoMipChi2_a2_);
    tree_->SetBranchAddress("phoMipTotEnergy_a2"        , &phoMipTotEnergy_a2_);
    tree_->SetBranchAddress("phoMipSlope_a2"            , &phoMipSlope_a2_);
    tree_->SetBranchAddress("phoMipIntercept_a2"        , &phoMipIntercept_a2_);
    tree_->SetBranchAddress("phoMipNhitCone_a2"         , &phoMipNhitCone_a2_);
    tree_->SetBranchAddress("phoMipIsHalo_a2"           , &phoMipIsHalo_a2_);
    tree_->SetBranchAddress("phoMatchType_a2"           , &phoMatchType_a2_);
    tree_->SetBranchAddress("phoMatchPt_a2"             , &phoMatchPt_a2_);
    tree_->SetBranchAddress("phoIsTrigger_a2"           , &phoIsTrigger_a2_);
    tree_->SetBranchAddress("pho3"                      , &phoPtr3_);
    tree_->SetBranchAddress("phoCombIso1_a3"            , &phoCombIso1_a3_);
    tree_->SetBranchAddress("phoCombIso2_a3"            , &phoCombIso2_a3_);
    tree_->SetBranchAddress("phoCombIso3_a3"            , &phoCombIso3_a3_ );
    tree_->SetBranchAddress("phoWorstIso_a3"            , &phoWorstIso_a3_);
    tree_->SetBranchAddress("phoPassEleVeto_a3"         , &phoPassEleVeto_a3_);
    tree_->SetBranchAddress("phoHasPixelSeed_a3"        , &phoHasPixelSeed_a3_);
    tree_->SetBranchAddress("phoCoviEtaiEta_a3"         , &phoCoviEtaiEta_a3_);
    tree_->SetBranchAddress("phoCoviPhiiPhi_a3"         , &phoCoviPhiiPhi_a3_);
    tree_->SetBranchAddress("phoR9_a3"                  , &phoR9_a3_);
    tree_->SetBranchAddress("phoSeedTime_a3"            , &phoSeedTime_a3_);
    tree_->SetBranchAddress("phoHadOverEm_a3"           , &phoHadOverEm_a3_);
    tree_->SetBranchAddress("phoLeadTimeSpan_a3"        , &phoLeadTimeSpan_a3_);
    tree_->SetBranchAddress("phoSubLeadTimeSpan_a3"     , &phoSubLeadTimeSpan_a3_);
    tree_->SetBranchAddress("phoMipChi2_a3"             , &phoMipChi2_a3_);
    tree_->SetBranchAddress("phoMipTotEnergy_a3"        , &phoMipTotEnergy_a3_);
    tree_->SetBranchAddress("phoMipSlope_a3"            , &phoMipSlope_a3_);
    tree_->SetBranchAddress("phoMipIntercept_a3"        , &phoMipIntercept_a3_);
    tree_->SetBranchAddress("phoMipNhitCone_a3"         , &phoMipNhitCone_a3_);
    tree_->SetBranchAddress("phoMipIsHalo_a3"           , &phoMipIsHalo_a3_);
    tree_->SetBranchAddress("phoMatchType_a3"           , &phoMatchType_a3_);
    tree_->SetBranchAddress("phoMatchPt_a3"             , &phoMatchPt_a3_);
    tree_->SetBranchAddress("phoIsTrigger_a3"           , &phoIsTrigger_a3_);

    tree_->SetBranchAddress("njets",         &njets_);
    tree_->SetBranchAddress("jet1",          &jetPtr1_);
    tree_->SetBranchAddress("jet1Btag",      &jet1Btag_);
    tree_->SetBranchAddress("jet2",          &jetPtr2_);
    tree_->SetBranchAddress("jet2Btag",      &jet2Btag_);
    tree_->SetBranchAddress("jet3",          &jetPtr3_);
    tree_->SetBranchAddress("jet3Btag",      &jet3Btag_);
    tree_->SetBranchAddress("jet4",          &jetPtr4_);
    tree_->SetBranchAddress("jet4Btag",      &jet4Btag_);

    tree_->SetBranchAddress("ntracks",       &ntracks_);
    tree_->SetBranchAddress("track1",        &trackPtr1_);
    tree_->SetBranchAddress("track2",        &trackPtr2_);
    tree_->SetBranchAddress("track3",        &trackPtr3_);

    tree_->SetBranchAddress("ncosmics",       &ncosmics_);
    tree_->SetBranchAddress("cosmic1",        &cosmicPtr1_);
    tree_->SetBranchAddress("etaOfSTA_c1",    &etaOfSTA_c1_);
    tree_->SetBranchAddress("phiOfSTA_c1",    &phiOfSTA_c1_);
    tree_->SetBranchAddress("cosmic2",        &cosmicPtr2_);
    tree_->SetBranchAddress("etaOfSTA_c2",    &etaOfSTA_c2_);
    tree_->SetBranchAddress("phiOfSTA_c2",    &phiOfSTA_c2_);
    tree_->SetBranchAddress("cosmic3",        &cosmicPtr3_);
    tree_->SetBranchAddress("etaOfSTA_c3",    &etaOfSTA_c3_);
    tree_->SetBranchAddress("phiOfSTA_c3",    &phiOfSTA_c3_);

    tree_->SetBranchAddress("Q",	          &Q_);
    tree_->SetBranchAddress("id1",	        &id1_);
    tree_->SetBranchAddress("x1",	          &x1_);
    tree_->SetBranchAddress("pdf1",	        &pdf1_);
    tree_->SetBranchAddress("id2",	        &id2_);
    tree_->SetBranchAddress("x2",	          &x2_);
    tree_->SetBranchAddress("pdf2",	        &pdf2_);
    tree_->SetBranchAddress("processId",    &processId_);
    tree_->SetBranchAddress("npu",	        &npu_);
    tree_->SetBranchAddress("npuPlusOne",   &npuPlusOne_);
    tree_->SetBranchAddress("npuMinusOne",  &npuMinusOne_);
    tree_->SetBranchAddress("rho",          &rho_);

    gErrorIgnoreLevel = currentState;
  }

  private:

  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* lepPtr3_;
  LorentzVector* phoPtr1_;
  LorentzVector* phoPtr2_;
  LorentzVector* phoPtr3_;
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
  LorentzVector* trackPtr1_;
  LorentzVector* trackPtr2_;
  LorentzVector* trackPtr3_;
  LorentzVector* cosmicPtr1_;
  LorentzVector* cosmicPtr2_;
  LorentzVector* cosmicPtr3_;

}; 

inline void 
MitGPTree::InitVariables(){
  // inizialize variables
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  nvtx_          = 0;
  scale1fb_      = 0;
  met_           = -999.;
  metPhi_        = -999.;
  sumEt_         = -999.;
  metSig_        = -999.;
  dstype_        = data;
  metCor_        = -999.;      
  metCorPhi_     = -999.;   

  nlep_          = 0;
  lep1_       	 = LorentzVector();
  lid1_          = 0;
  etaOfSTA_l1_   = -999.;
  phiOfSTA_l1_   = -999.;
  lep2_       	 = LorentzVector();
  lid2_          = 0;
  etaOfSTA_l2_   = -999.;
  phiOfSTA_l2_   = -999.;
  lep3_       	 = LorentzVector();
  lid3_          = 0;
  etaOfSTA_l3_   = -999.;
  phiOfSTA_l3_   = -999.;

  nphotons_      = 0;
  pho1_       	 = LorentzVector();
  phoCombIso1_a1_ = -1.0;
  phoCombIso2_a1_ = -1.0;
  phoCombIso3_a1_ = -1.0;
  phoWorstIso_a1_ = -1.0;
  phoPassEleVeto_a1_ = false;
  phoHasPixelSeed_a1_ = false;
  phoCoviEtaiEta_a1_ = -1.0;
  phoCoviPhiiPhi_a1_ = -1.0;
  phoR9_a1_ = -1.0;
  phoSeedTime_a1_ = -999.0;
  phoHadOverEm_a1_ = -1.0;  
  phoLeadTimeSpan_a1_ = -999.0;
  phoSubLeadTimeSpan_a1_ = -999.0;
  phoMipChi2_a1_ = -1.0;
  phoMipTotEnergy_a1_ = -1.0;
  phoMipSlope_a1_ = -999.0;
  phoMipIntercept_a1_ = -999.0;
  phoMipNhitCone_a1_ = 0;
  phoMipIsHalo_a1_ = false;
  phoMatchType_a1_ = -1.0;
  phoMatchPt_a1_ = -1.0;
  phoIsTrigger_a1_ = false;
  pho2_       	 = LorentzVector();
  phoCombIso1_a2_ = -1.0;
  phoCombIso2_a2_ = -1.0;
  phoCombIso3_a2_ = -1.0;
  phoWorstIso_a2_ = -1.0;
  phoPassEleVeto_a2_ = false;
  phoHasPixelSeed_a2_ = false;
  phoCoviEtaiEta_a2_ = -1.0;
  phoCoviPhiiPhi_a2_ = -1.0;
  phoR9_a2_ = -1.0;
  phoSeedTime_a2_ = -999.0;
  phoHadOverEm_a2_ = -1.0;  
  phoLeadTimeSpan_a2_ = -999.0;
  phoSubLeadTimeSpan_a2_ = -999.0;
  phoMipChi2_a2_ = -1.0;
  phoMipTotEnergy_a2_ = -1.0;
  phoMipSlope_a2_ = -999.0;
  phoMipIntercept_a2_ = -999.0;
  phoMipNhitCone_a2_ = 0;
  phoMipIsHalo_a2_ = false;
  phoMatchType_a2_ = -1.0;
  phoMatchPt_a2_ = -1.0;
  phoIsTrigger_a2_ = false;
  pho3_       	 = LorentzVector();
  phoCombIso1_a3_ = -1.0;
  phoCombIso2_a3_ = -1.0;
  phoCombIso3_a3_ = -1.0;
  phoWorstIso_a3_ = -1.0;
  phoPassEleVeto_a3_ = false;
  phoHasPixelSeed_a3_ = false;
  phoCoviEtaiEta_a3_ = -1.0;
  phoCoviPhiiPhi_a3_ = -1.0;
  phoR9_a3_ = -1.0;
  phoSeedTime_a3_ = -999.0;
  phoHadOverEm_a3_ = -1.0;  
  phoLeadTimeSpan_a3_ = -999.0;
  phoSubLeadTimeSpan_a3_ = -999.0;
  phoMipChi2_a3_ = -1.0;
  phoMipTotEnergy_a3_ = -1.0;
  phoMipSlope_a3_ = -999.0;
  phoMipIntercept_a3_ = -999.0;
  phoMipNhitCone_a3_ = 0;
  phoMipIsHalo_a3_ = false;
  phoMatchType_a3_ = -1.0;
  phoMatchPt_a3_ = -1.0;
  phoIsTrigger_a3_ = false;

  njets_ = 0;
  jet1_     = LorentzVector();
  jet1Btag_ = -999.;
  jet2_     = LorentzVector();
  jet2Btag_ = -999.;
  jet3_     = LorentzVector();
  jet3Btag_ = -999.;
  jet4_     = LorentzVector();
  jet4Btag_ = -999.;

  ntracks_ = 0;
  track1_ = LorentzVector();
  track2_ = LorentzVector();
  track3_ = LorentzVector();

  ncosmics_ = 0;
  cosmic1_ = LorentzVector();
  etaOfSTA_c1_   = -999.;
  phiOfSTA_c1_   = -999.;
  cosmic2_ = LorentzVector();
  etaOfSTA_c2_   = -999.;
  phiOfSTA_c2_   = -999.;
  cosmic3_ = LorentzVector();
  etaOfSTA_c3_   = -999.;
  phiOfSTA_c3_   = -999.;

  Q_		 = -999.;
  id1_  	 = -999.;
  x1_		 = -999.;
  pdf1_ 	 = -999.;  
  id2_  	 = -999.;  
  x2_		 = -999.;
  pdf2_ 	 = -999.;  
  processId_	 = 0;
  npu_           = -999.;
  npuPlusOne_    = -999.;
  npuMinusOne_   = -999.;
  rho_           = -1.0;
}

#endif
