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

  /// variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   nvtx_;
  float          scale1fb_;
  float          met_;
  float          metPhi_;
  float          sumEt_;
  float          metSig_;
  DataType       dstype_;

  unsigned int   nlep_;
  LorentzVector  lep1_;
  int            lid1_;
  LorentzVector  lep2_;
  int            lid2_;
  LorentzVector  lep3_;
  int            lid3_;

  unsigned int   nphotons_;
  LorentzVector  pho1_;
  float phoHCALisoDR03_a1_;
  float phoECALisoDR03_a1_;
  float phoHollowConeTKisoDR03_a1_;
  float phoHCALisoDR04_a1_;
  float phoECALisoDR04_a1_;
  float phoHollowConeTKisoDR04_a1_;
  float phoCoviEtaiEta_a1_;
  float phoR9_a1_;
  float phoSeedTime_a1_;
  float phoHadOverEm_a1_;
  LorentzVector  pho2_;
  float phoHCALisoDR03_a2_;
  float phoECALisoDR03_a2_;
  float phoHollowConeTKisoDR03_a2_;
  float phoHCALisoDR04_a2_;
  float phoECALisoDR04_a2_;
  float phoHollowConeTKisoDR04_a2_;
  float phoCoviEtaiEta_a2_;
  float phoR9_a2_;
  float phoSeedTime_a2_;
  float phoHadOverEm_a2_;
  LorentzVector  pho3_;
  float phoHCALisoDR03_a3_;
  float phoECALisoDR03_a3_;
  float phoHollowConeTKisoDR03_a3_;
  float phoHCALisoDR04_a3_;
  float phoECALisoDR04_a3_;
  float phoHollowConeTKisoDR04_a3_;
  float phoCoviEtaiEta_a3_;
  float phoR9_a3_;
  float phoSeedTime_a3_;
  float phoHadOverEm_a3_;
  LorentzVector  pho4_;
  float phoHCALisoDR03_a4_;
  float phoECALisoDR03_a4_;
  float phoHollowConeTKisoDR03_a4_;
  float phoHCALisoDR04_a4_;
  float phoECALisoDR04_a4_;
  float phoHollowConeTKisoDR04_a4_;
  float phoCoviEtaiEta_a4_;
  float phoR9_a4_;
  float phoSeedTime_a4_;
  float phoHadOverEm_a4_;

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

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitGPTree():
    lepPtr1_(&lep1_),lepPtr2_(&lep2_),lepPtr3_(&lep3_),
    phoPtr1_(&pho1_),phoPtr2_(&pho2_),phoPtr3_(&pho3_),phoPtr4_(&pho4_),
    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),
    trackPtr1_(&track1_),trackPtr2_(&track2_),trackPtr3_(&track3_){}
  /// default destructor
  ~MitGPTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitGPTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job DMTree0/1
    // type == 0/1 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    f_ = TFile::Open(file);
    assert(f_);
    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->Get("DMTree0"));
    else if(type == 1) tree_ = dynamic_cast<TTree*>(f_->Get("DMTree1"));
    else               tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
    assert(tree_);
  }

  /// create a MitGPTree
  void CreateTree(int type = -1){
    assert(type==type); // just to suppress warnings
    // to create three different ntuples in the same job DMTree0/1
    // type == 0/1 add all variables
    // type = -1 (default) add a minimum set of variables with tree as name
    if     (type == 0) tree_ = new TTree("DMTree0","Smurf ntuples");
    else if(type == 1) tree_ = new TTree("DMTree1","Smurf ntuples");
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

    tree_->Branch("nlep"         , &nlep_         ,   "nlep/i");
    tree_->Branch("lep1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lid1"         , &lid1_         ,   "lid1/I");
    tree_->Branch("lep2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lid2"         , &lid2_         ,   "lid2/I");
    tree_->Branch("lep3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
    tree_->Branch("lid3"         , &lid3_         ,   "lid3/I");

    tree_->Branch("nphotons"     , &nphotons_     ,   "nphotons/i");
    tree_->Branch("pho1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr1_);
    tree_->Branch("phoHCALisoDR03_a1"         , &phoHCALisoDR03_a1_         , "phoHCALisoDR03_a1/F");
    tree_->Branch("phoECALisoDR03_a1"         , &phoECALisoDR03_a1_         , "phoECALisoDR03_a1/F");
    tree_->Branch("phoHollowConeTKisoDR03_a1" , &phoHollowConeTKisoDR03_a1_ , "phoHollowConeTKisoDR03_a1/F");
    tree_->Branch("phoHCALisoDR04_a1"         , &phoHCALisoDR04_a1_         , "phoHCALisoDR04_a1/F");
    tree_->Branch("phoECALisoDR04_a1"         , &phoECALisoDR04_a1_         , "phoECALisoDR04_a1/F");
    tree_->Branch("phoHollowConeTKisoDR04_a1" , &phoHollowConeTKisoDR04_a1_ , "phoHollowConeTKisoDR04_a1/F");
    tree_->Branch("phoCoviEtaiEta_a1"         , &phoCoviEtaiEta_a1_         , "phoCoviEtaiEta_a1/F");
    tree_->Branch("phoR9_a1"                  , &phoR9_a1_                  , "phoR9_a1/F");
    tree_->Branch("phoSeedTime_a1"            , &phoSeedTime_a1_            , "phoSeedTime_a1/F");
    tree_->Branch("phoHadOverEm_a1"           , &phoHadOverEm_a1_           , "phoHadOverEm_a1/F");
    tree_->Branch("pho2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr2_);
    tree_->Branch("phoHCALisoDR03_a2"         , &phoHCALisoDR03_a2_         , "phoHCALisoDR03_a2/F");
    tree_->Branch("phoECALisoDR03_a2"         , &phoECALisoDR03_a2_         , "phoECALisoDR03_a2/F");
    tree_->Branch("phoHollowConeTKisoDR03_a2" , &phoHollowConeTKisoDR03_a2_ , "phoHollowConeTKisoDR03_a2/F");
    tree_->Branch("phoHCALisoDR04_a2"         , &phoHCALisoDR04_a2_         , "phoHCALisoDR04_a2/F");
    tree_->Branch("phoECALisoDR04_a2"         , &phoECALisoDR04_a2_         , "phoECALisoDR04_a2/F");
    tree_->Branch("phoHollowConeTKisoDR04_a2" , &phoHollowConeTKisoDR04_a2_ , "phoHollowConeTKisoDR04_a2/F");
    tree_->Branch("phoCoviEtaiEta_a2"         , &phoCoviEtaiEta_a2_         , "phoCoviEtaiEta_a2/F");
    tree_->Branch("phoR9_a2"                  , &phoR9_a2_                  , "phoR9_a2/F");
    tree_->Branch("phoSeedTime_a2"            , &phoSeedTime_a2_            , "phoSeedTime_a2/F");
    tree_->Branch("phoHadOverEm_a2"           , &phoHadOverEm_a2_           , "phoHadOverEm_a2/F");
    tree_->Branch("pho3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr3_);
    tree_->Branch("phoHCALisoDR03_a3"         , &phoHCALisoDR03_a3_         , "phoHCALisoDR03_a3/F");
    tree_->Branch("phoECALisoDR03_a3"         , &phoECALisoDR03_a3_         , "phoECALisoDR03_a3/F");
    tree_->Branch("phoHollowConeTKisoDR03_a3" , &phoHollowConeTKisoDR03_a3_ , "phoHollowConeTKisoDR03_a3/F");
    tree_->Branch("phoHCALisoDR04_a3"         , &phoHCALisoDR04_a3_         , "phoHCALisoDR04_a3/F");
    tree_->Branch("phoECALisoDR04_a3"         , &phoECALisoDR04_a3_         , "phoECALisoDR04_a3/F");
    tree_->Branch("phoHollowConeTKisoDR04_a3" , &phoHollowConeTKisoDR04_a3_ , "phoHollowConeTKisoDR04_a3/F");
    tree_->Branch("phoCoviEtaiEta_a3"         , &phoCoviEtaiEta_a3_         , "phoCoviEtaiEta_a3/F");
    tree_->Branch("phoR9_a3"                  , &phoR9_a3_                  , "phoR9_a3/F");
    tree_->Branch("phoSeedTime_a3"            , &phoSeedTime_a3_            , "phoSeedTime_a3/F");
    tree_->Branch("phoHadOverEm_a3"           , &phoHadOverEm_a3_           , "phoHadOverEm_a3/F");
    tree_->Branch("pho4"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr4_);
    tree_->Branch("phoHCALisoDR03_a4"         , &phoHCALisoDR03_a4_         , "phoHCALisoDR03_a4/F");
    tree_->Branch("phoECALisoDR03_a4"         , &phoECALisoDR03_a4_         , "phoECALisoDR03_a4/F");
    tree_->Branch("phoHollowConeTKisoDR03_a4" , &phoHollowConeTKisoDR03_a4_ , "phoHollowConeTKisoDR03_a4/F");
    tree_->Branch("phoHCALisoDR04_a4"         , &phoHCALisoDR04_a4_         , "phoHCALisoDR04_a4/F");
    tree_->Branch("phoECALisoDR04_a4"         , &phoECALisoDR04_a4_         , "phoECALisoDR04_a4/F");
    tree_->Branch("phoHollowConeTKisoDR04_a4" , &phoHollowConeTKisoDR04_a4_ , "phoHollowConeTKisoDR04_a4/F");
    tree_->Branch("phoCoviEtaiEta_a4"         , &phoCoviEtaiEta_a4_         , "phoCoviEtaiEta_a4/F");
    tree_->Branch("phoR9_a4"                  , &phoR9_a4_                  , "phoR9_a4/F");
    tree_->Branch("phoSeedTime_a4"            , &phoSeedTime_a4_            , "phoSeedTime_a4/F");
    tree_->Branch("phoHadOverEm_a4"           , &phoHadOverEm_a4_           , "phoHadOverEm_a4/F");

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

    tree_->SetBranchAddress("nlep",          &nlep_);
    tree_->SetBranchAddress("lep1",          &lepPtr1_);
    tree_->SetBranchAddress("lid1",          &lid1_);
    tree_->SetBranchAddress("lep2",          &lepPtr2_);
    tree_->SetBranchAddress("lid2",          &lid2_);
    tree_->SetBranchAddress("lep3",          &lepPtr3_);
    tree_->SetBranchAddress("lid3",          &lid3_);

    tree_->SetBranchAddress("nphotons"                  , &nphotons_);
    tree_->SetBranchAddress("pho1"                      , &phoPtr1_);
    tree_->SetBranchAddress("phoHCALisoDR03_a1"         , &phoHCALisoDR03_a1_);
    tree_->SetBranchAddress("phoECALisoDR03_a1"         , &phoECALisoDR03_a1_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR03_a1" , &phoHollowConeTKisoDR03_a1_ );
    tree_->SetBranchAddress("phoHCALisoDR04_a1"         , &phoHCALisoDR04_a1_);
    tree_->SetBranchAddress("phoECALisoDR04_a1"         , &phoECALisoDR04_a1_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR04_a1" , &phoHollowConeTKisoDR04_a1_ );
    tree_->SetBranchAddress("phoCoviEtaiEta_a1"         , &phoCoviEtaiEta_a1_);
    tree_->SetBranchAddress("phoR9_a1"                  , &phoR9_a1_);
    tree_->SetBranchAddress("phoSeedTime_a1"            , &phoSeedTime_a1_);
    tree_->SetBranchAddress("phoHadOverEm_a1"           , &phoHadOverEm_a1_);
    tree_->SetBranchAddress("pho1"                      , &phoPtr1_);
    tree_->SetBranchAddress("phoHCALisoDR03_a2"         , &phoHCALisoDR03_a2_);
    tree_->SetBranchAddress("phoECALisoDR03_a2"         , &phoECALisoDR03_a2_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR03_a2" , &phoHollowConeTKisoDR03_a2_ );
    tree_->SetBranchAddress("phoHCALisoDR04_a2"         , &phoHCALisoDR04_a2_);
    tree_->SetBranchAddress("phoECALisoDR04_a2"         , &phoECALisoDR04_a2_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR04_a2" , &phoHollowConeTKisoDR04_a2_ );
    tree_->SetBranchAddress("phoCoviEtaiEta_a2"         , &phoCoviEtaiEta_a2_);
    tree_->SetBranchAddress("phoR9_a2"                  , &phoR9_a2_);
    tree_->SetBranchAddress("phoSeedTime_a2"            , &phoSeedTime_a2_);
    tree_->SetBranchAddress("phoHadOverEm_a2"           , &phoHadOverEm_a2_);
    tree_->SetBranchAddress("pho3"                      , &phoPtr3_);
    tree_->SetBranchAddress("phoHCALisoDR03_a3"         , &phoHCALisoDR03_a3_);
    tree_->SetBranchAddress("phoECALisoDR03_a3"         , &phoECALisoDR03_a3_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR03_a3" , &phoHollowConeTKisoDR03_a3_ );
    tree_->SetBranchAddress("phoHCALisoDR04_a3"         , &phoHCALisoDR04_a3_);
    tree_->SetBranchAddress("phoECALisoDR04_a3"         , &phoECALisoDR04_a3_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR04_a3" , &phoHollowConeTKisoDR04_a3_ );
    tree_->SetBranchAddress("phoCoviEtaiEta_a3"         , &phoCoviEtaiEta_a3_);
    tree_->SetBranchAddress("phoR9_a3"                  , &phoR9_a3_);
    tree_->SetBranchAddress("phoSeedTime_a3"            , &phoSeedTime_a3_);
    tree_->SetBranchAddress("phoHadOverEm_a3"           , &phoHadOverEm_a3_);
    tree_->SetBranchAddress("pho4"                      , &phoPtr4_);
    tree_->SetBranchAddress("phoHCALisoDR03_a4"         , &phoHCALisoDR03_a4_);
    tree_->SetBranchAddress("phoECALisoDR03_a4"         , &phoECALisoDR03_a4_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR03_a4" , &phoHollowConeTKisoDR03_a4_ );
    tree_->SetBranchAddress("phoHCALisoDR04_a4"         , &phoHCALisoDR04_a4_);
    tree_->SetBranchAddress("phoECALisoDR04_a4"         , &phoECALisoDR04_a4_);
    tree_->SetBranchAddress("phoHollowConeTKisoDR04_a4" , &phoHollowConeTKisoDR04_a4_ );
    tree_->SetBranchAddress("phoCoviEtaiEta_a4"         , &phoCoviEtaiEta_a4_);
    tree_->SetBranchAddress("phoR9_a4"                  , &phoR9_a4_);
    tree_->SetBranchAddress("phoSeedTime_a4"            , &phoSeedTime_a4_);
    tree_->SetBranchAddress("phoHadOverEm_a4"           , &phoHadOverEm_a4_);

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

    tree_->SetBranchAddress("Q",	     &Q_);
    tree_->SetBranchAddress("id1",	     &id1_);
    tree_->SetBranchAddress("x1",	     &x1_);
    tree_->SetBranchAddress("pdf1",	     &pdf1_);
    tree_->SetBranchAddress("id2",	     &id2_);
    tree_->SetBranchAddress("x2",	     &x2_);
    tree_->SetBranchAddress("pdf2",	     &pdf2_);
    tree_->SetBranchAddress("processId",     &processId_);
    tree_->SetBranchAddress("npu",	     &npu_);
    tree_->SetBranchAddress("npuPlusOne",    &npuPlusOne_);
    tree_->SetBranchAddress("npuMinusOne",   &npuMinusOne_);

    gErrorIgnoreLevel = currentState;
  }

  private:

  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* lepPtr3_;
  LorentzVector* phoPtr1_;
  LorentzVector* phoPtr2_;
  LorentzVector* phoPtr3_;
  LorentzVector* phoPtr4_;
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
  LorentzVector* trackPtr1_;
  LorentzVector* trackPtr2_;
  LorentzVector* trackPtr3_;
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

  nlep_          = 0;
  lep1_       	 = LorentzVector();
  lid1_          = 0;
  lep2_       	 = LorentzVector();
  lid2_          = 0;
  lep3_       	 = LorentzVector();
  lid3_          = 0;

  nphotons_      = 0;
  pho1_       	 = LorentzVector();
  phoHCALisoDR03_a1_ = 0.0;
  phoECALisoDR03_a1_ = 0.0;
  phoHollowConeTKisoDR03_a1_ = 0.0;
  phoHCALisoDR04_a1_ = 0.0;
  phoECALisoDR04_a1_ = 0.0;
  phoHollowConeTKisoDR04_a1_ = 0.0;
  phoCoviEtaiEta_a1_ = 0.0;
  phoR9_a1_ = 0.0;
  phoSeedTime_a1_ = 0.0;
  phoHadOverEm_a1_ = 0.0;  
  pho2_       	 = LorentzVector();
  phoHCALisoDR03_a2_ = 0.0;
  phoECALisoDR03_a2_ = 0.0;
  phoHollowConeTKisoDR03_a2_ = 0.0;
  phoHCALisoDR04_a2_ = 0.0;
  phoECALisoDR04_a2_ = 0.0;
  phoHollowConeTKisoDR04_a2_ = 0.0;
  phoCoviEtaiEta_a2_ = 0.0;
  phoR9_a2_ = 0.0;
  phoSeedTime_a2_ = 0.0;
  phoHadOverEm_a2_ = 0.0;  
  pho3_       	 = LorentzVector();
  phoHCALisoDR03_a3_ = 0.0;
  phoECALisoDR03_a3_ = 0.0;
  phoHollowConeTKisoDR03_a3_ = 0.0;
  phoHCALisoDR04_a3_ = 0.0;
  phoECALisoDR04_a3_ = 0.0;
  phoHollowConeTKisoDR04_a3_ = 0.0;
  phoCoviEtaiEta_a3_ = 0.0;
  phoR9_a3_ = 0.0;
  phoSeedTime_a3_ = 0.0;
  phoHadOverEm_a3_ = 0.0;  
  pho4_       	 = LorentzVector();
  phoHCALisoDR03_a4_ = 0.0;
  phoECALisoDR03_a4_ = 0.0;
  phoHollowConeTKisoDR03_a4_ = 0.0;
  phoHCALisoDR04_a4_ = 0.0;
  phoECALisoDR04_a4_ = 0.0;
  phoHollowConeTKisoDR04_a4_ = 0.0;
  phoCoviEtaiEta_a4_ = 0.0;
  phoR9_a4_ = 0.0;
  phoSeedTime_a4_ = 0.0;
  phoHadOverEm_a4_ = 0.0;  

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
}

#endif
