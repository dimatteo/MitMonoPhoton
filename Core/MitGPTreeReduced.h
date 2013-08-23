#ifndef MitGPTreeReduced_H
#define MitGPTreeReduced_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitGPTreeReduced {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

  /// variables
  bool           isData_;

  float          evt_weight_;
  float          hlt_weight_;
  float          pu_weight_;
  float          kf_weight_;

  unsigned int   nvtx_;
  float          met_;
  float          metPhi_;
  float          sumEt_;
  float          metSig_;

  float          phoEt_;
  float          phoEta_;
  float          phoPhi_;

  float          phoMetDeltaPhi_;

  float          jetEt_;
  float          jetEta_;
  float          jetPhi_;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitGPTreeReduced(){}
  /// default destructor
  ~MitGPTreeReduced(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitGPTreeReduced
  void LoadTree(const char* file, const char* tree_name){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(tree_name));
    assert(tree_);
  }

  /// create a MitGPTreeReduced
  void CreateTree(const char* tree_name){
    tree_ = new TTree(tree_name,"Ana ntuple");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("isData"            , &isData_            ,   "isData/B");

    tree_->Branch("evt_weight"        , &evt_weight_        ,   "evt_weight/F");
    tree_->Branch("hlt_weight"        , &hlt_weight_        ,   "hlt_weight/F");
    tree_->Branch("pu_weight"         , &pu_weight_         ,   "pu_weight/F");
    tree_->Branch("kf_weight"         , &kf_weight_         ,   "kf_weight/F");
  
    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("met"          , &met_          ,   "met/F");
    tree_->Branch("metPhi"       , &metPhi_       ,   "metPhi/F");
    tree_->Branch("sumEt"        , &sumEt_        ,   "sumEt/F");
    tree_->Branch("metSig"       , &metSig_       ,   "metSig/F");
  
    tree_->Branch("phoEt"         , &phoEt_         ,   "phoEt/F");
    tree_->Branch("phoEta"        , &phoEta_        ,   "phoEta/F");
    tree_->Branch("phoPhi"        , &phoPhi_        ,   "phoPhi/F");
  
    tree_->Branch("phoMetDeltaPhi", &phoMetDeltaPhi_,   "phoMetDeltaPhi/F");
  
    tree_->Branch("jetEt"         , &jetEt_         ,   "jetEt/F");
    tree_->Branch("jetEta"        , &jetEta_        ,   "jetEta/F");
    tree_->Branch("jetPhi"        , &jetPhi_        ,   "jetPhi/F");

  }

  // initialze a MitGPTreeReduced
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
    tree_->SetBranchAddress("kf_weight"         , &kf_weight_ );
  
    tree_->SetBranchAddress("nvtx"         , &nvtx_      );
    tree_->SetBranchAddress("met"          , &met_       );
    tree_->SetBranchAddress("metPhi"       , &metPhi_    );
    tree_->SetBranchAddress("sumEt"        , &sumEt_     );
    tree_->SetBranchAddress("metSig"       , &metSig_    );
  
    tree_->SetBranchAddress("phoEt"         , &phoEt_     );
    tree_->SetBranchAddress("phoEta"        , &phoEta_    );
    tree_->SetBranchAddress("phoPhi"        , &phoPhi_    );
  
    tree_->SetBranchAddress("phoMetDeltaPhi", &phoMetDeltaPhi_ );
  
    tree_->SetBranchAddress("jetEt"         , &jetEt_    );
    tree_->SetBranchAddress("jetEta"        , &jetEta_   );
    tree_->SetBranchAddress("jetPhi"        , &jetPhi_   );

    gErrorIgnoreLevel = currentState;
  }

  private:

}; 

inline void 
MitGPTreeReduced::InitVariables(){
  // inizialize variables
  isData_ = false;
  
  evt_weight_ = -1.;
  hlt_weight_ = -1.;
  pu_weight_  = -1.;
  kf_weight_  = -1.;

  nvtx_  = 0;
  met_   = -1.;
  metPhi_= -100.;
  sumEt_ = -1.;
  metSig_= -1.;

  phoEt_ = -1.;
  phoEta_= -100.;
  phoPhi_= -100.;

  phoMetDeltaPhi_ = -100.;

  jetEt_ = -1.;
  jetEta_= -100.;
  jetPhi_= -100.;
}

#endif
