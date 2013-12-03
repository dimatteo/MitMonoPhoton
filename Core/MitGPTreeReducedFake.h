#ifndef MitGPTreeReducedFake_H
#define MitGPTreeReducedFake_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitGPTreeReducedFake {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

  /// variables
  bool           isData_;

  float          evt_weight_;
  float          hlt_weight_;
  float          pu_weight_;
  float          pu_weightup_;
  float          pu_weightdo_;
  float          kf_weight_;

  unsigned int   nvtx_;
  float          metCor_;
  float          metCorPhi_;
  float          metSig_;

  int            nphotons_;
  float          phoEt_;
  float          phoEta_;
  float          phoPhi_;
  float          phoCombIso1_; 
  float          phoCombIso2_; 
  float          phoCombIso3_; 
  float          phoCoviEtaiEta_; 
  float          phoR9_; 
  int            phoIsFake_; 

  float          phoMetDeltaPhi_;
  float          jetMetDeltaPhi_;
  float          phoJetDeltaPhi_;
  
  //for MET syst
  unsigned int   nalljets_;
  float  jet1Pt_;
  float  jet1Eta_;
  float  jet1Phi_;
  float  jet2Pt_;
  float  jet2Eta_;
  float  jet2Phi_;
  float  jet3Pt_;
  float  jet3Eta_;
  float  jet3Phi_;
  float  jet4Pt_;
  float  jet4Eta_;
  float  jet4Phi_;

  //for Lep syst
  unsigned int   nlep_;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitGPTreeReducedFake(){}
  /// default destructor
  ~MitGPTreeReducedFake(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitGPTreeReducedFake
  void LoadTree(const char* file, const char* tree_name){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(tree_name));
    assert(tree_);
  }

  /// create a MitGPTreeReducedFake
  void CreateTree(const char* tree_name){
    tree_ = new TTree(tree_name,"Ana ntuple");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("isData"            , &isData_            ,   "isData/B");

    tree_->Branch("evt_weight"        , &evt_weight_        ,   "evt_weight/F");
    tree_->Branch("hlt_weight"        , &hlt_weight_        ,   "hlt_weight/F");
    tree_->Branch("pu_weight"         , &pu_weight_         ,   "pu_weight/F");
    tree_->Branch("pu_weightup"       , &pu_weightup_       ,   "pu_weightup/F");
    tree_->Branch("pu_weightdo"       , &pu_weightdo_       ,   "pu_weightdo/F");
    tree_->Branch("kf_weight"         , &kf_weight_         ,   "kf_weight/F");
  
    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("metCor"       , &metCor_       ,   "metCor/F");
    tree_->Branch("metCorPhi"    , &metCorPhi_    ,   "metCorPhi/F");
    tree_->Branch("metSig"       , &metSig_       ,   "metSig/F");
  
    tree_->Branch("nphotons"        , &nphotons_        ,   "nphotons/i");
    tree_->Branch("phoEt"           , &phoEt_           ,   "phoEt/F");
    tree_->Branch("phoEta"          , &phoEta_          ,   "phoEta/F");
    tree_->Branch("phoPhi"          , &phoPhi_          ,   "phoPhi/F");
    tree_->Branch("phoCombIso1"     , &phoCombIso1_     ,   "phoCombIso1/F"); 
    tree_->Branch("phoCombIso2"     , &phoCombIso2_     ,   "phoCombIso2/F"); 
    tree_->Branch("phoCombIso3"     , &phoCombIso3_     ,   "phoCombIso3/F"); 
    tree_->Branch("phoCoviEtaiEta"  , &phoCoviEtaiEta_  ,   "phoCoviEtaiEta/F"); 
    tree_->Branch("phoR9"           , &phoR9_           ,   "phoR9/F"); 
    tree_->Branch("phoIsFake"       , &phoIsFake_       ,   "phoIsFake/i"); 

    tree_->Branch("phoMetDeltaPhi", &phoMetDeltaPhi_,   "phoMetDeltaPhi/F");
    tree_->Branch("jetMetDeltaPhi", &jetMetDeltaPhi_,   "jetMetDeltaPhi/F");
    tree_->Branch("phoJetDeltaPhi", &phoJetDeltaPhi_,   "phoJetDeltaPhi/F");
  
    tree_->Branch("nalljets"          , &nalljets_          ,   "nalljets/i");
    tree_->Branch("jet1Pt"            , &jet1Pt_            ,   "jet1Pt/F");
    tree_->Branch("jet1Eta"           , &jet1Eta_           ,   "jet1Eta/F");
    tree_->Branch("jet1Phi"           , &jet1Phi_           ,   "jet1Phi/F");
    tree_->Branch("jet2Pt"            , &jet2Pt_            ,   "jet2Pt/F");
    tree_->Branch("jet2Eta"           , &jet2Eta_           ,   "jet2Eta/F");
    tree_->Branch("jet2Phi"           , &jet2Phi_           ,   "jet2Phi/F");
    tree_->Branch("jet3Pt"            , &jet3Pt_            ,   "jet3Pt/F");
    tree_->Branch("jet3Eta"           , &jet3Eta_           ,   "jet3Eta/F");
    tree_->Branch("jet3Phi"           , &jet3Phi_           ,   "jet3Phi/F");
    tree_->Branch("jet4Pt"            , &jet4Pt_            ,   "jet4Pt/F");
    tree_->Branch("jet4Eta"           , &jet4Eta_           ,   "jet4Eta/F");
    tree_->Branch("jet4Phi"           , &jet4Phi_           ,   "jet4Phi/F");

    tree_->Branch("nlep"              , &nlep_              ,   "nlep/i");

  }

  // initialze a MitGPTreeReducedFake
  void InitTree(){
    assert(tree_);
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("isData"            , &isData_);

    tree_->SetBranchAddress("evt_weight"        , &evt_weight_ );
    tree_->SetBranchAddress("hlt_weight"        , &hlt_weight_ );
    tree_->SetBranchAddress("pu_weight"         , &pu_weight_  );
    tree_->SetBranchAddress("pu_weightup"       , &pu_weightup_  );
    tree_->SetBranchAddress("pu_weightdo"       , &pu_weightdo_  );
    tree_->SetBranchAddress("kf_weight"         , &kf_weight_ );
  
    tree_->SetBranchAddress("nvtx"         , &nvtx_      );
    tree_->SetBranchAddress("metCor"       , &metCor_    );
    tree_->SetBranchAddress("metCorPhi"    , &metCorPhi_ );
    tree_->SetBranchAddress("metSig"       , &metSig_    );
  
    tree_->SetBranchAddress("nphotons"       , &nphotons_);
    tree_->SetBranchAddress("phoEt"          , &phoEt_     );
    tree_->SetBranchAddress("phoEta"         , &phoEta_    );
    tree_->SetBranchAddress("phoPhi"         , &phoPhi_    );
    tree_->SetBranchAddress("phoCombIso1"    , &phoCombIso1_); 
    tree_->SetBranchAddress("phoCombIso2"    , &phoCombIso2_); 
    tree_->SetBranchAddress("phoCombIso3"    , &phoCombIso3_); 
    tree_->SetBranchAddress("phoCoviEtaiEta" , &phoCoviEtaiEta_); 
    tree_->SetBranchAddress("phoR9"          , &phoR9_          ); 
    tree_->SetBranchAddress("phoIsFake"      , &phoIsFake_      ); 
  
    tree_->SetBranchAddress("phoMetDeltaPhi", &phoMetDeltaPhi_ );
    tree_->SetBranchAddress("jetMetDeltaPhi", &jetMetDeltaPhi_ );
    tree_->SetBranchAddress("phoJetDeltaPhi", &phoJetDeltaPhi_ );
  
    tree_->SetBranchAddress("nalljets"          , &nalljets_          );
    tree_->SetBranchAddress("jet1Pt"            , &jet1Pt_            );
    tree_->SetBranchAddress("jet1Eta"           , &jet1Eta_           );
    tree_->SetBranchAddress("jet1Phi"           , &jet1Phi_           );
    tree_->SetBranchAddress("jet2Pt"            , &jet2Pt_            );
    tree_->SetBranchAddress("jet2Eta"           , &jet2Eta_           );
    tree_->SetBranchAddress("jet2Phi"           , &jet2Phi_           );
    tree_->SetBranchAddress("jet3Pt"            , &jet3Pt_            );
    tree_->SetBranchAddress("jet3Eta"           , &jet3Eta_           );
    tree_->SetBranchAddress("jet3Phi"           , &jet3Phi_           );
    tree_->SetBranchAddress("jet4Pt"            , &jet4Pt_            );
    tree_->SetBranchAddress("jet4Eta"           , &jet4Eta_           );
    tree_->SetBranchAddress("jet4Phi"           , &jet4Phi_           );

    tree_->SetBranchAddress("nlep"              , &nlep_            );

    gErrorIgnoreLevel = currentState;
  }

  private:

}; 

inline void 
MitGPTreeReducedFake::InitVariables(){
  // inizialize variables
  isData_ = false;
  
  evt_weight_ = -1.;
  hlt_weight_ = -1.;
  pu_weight_  = -1.;
  pu_weightup_ = -1.;
  pu_weightdo_ = -1.;
  kf_weight_  = -1.;

  nvtx_  = 0;
  metCor_   = -1.;
  metCorPhi_= -100.;
  metSig_= -1.;

  nphotons_ = 0;
  phoEt_ = -1.;
  phoEta_= -100.;
  phoPhi_= -100.;
  phoCombIso1_ = -1.;
  phoCombIso2_ = -1.; 
  phoCombIso3_ = -1.;
  phoCoviEtaiEta_ = -1.;
  phoR9_          = -1.;
  phoIsFake_      = -1.; 

  phoMetDeltaPhi_ = -100.;
  jetMetDeltaPhi_ = -100.;
  phoJetDeltaPhi_ = -100.;
  
  nalljets_ = 0;
  jet1Pt_ = -1.;
  jet1Eta_ = -100.;
  jet1Phi_ = -100.;
  jet2Pt_ = -1.;
  jet2Eta_ = -100.;
  jet2Phi_ = -100.;
  jet3Pt_ = -1.;
  jet3Eta_ = -100.;
  jet3Phi_ = -100.;
  jet4Pt_ = -1.;
  jet4Eta_ = -100.;
  jet4Phi_ = -100.;

  nlep_ = 0;

}

#endif
