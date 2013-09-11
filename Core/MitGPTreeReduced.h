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
  float          pu_weightup_;
  float          pu_weightdo_;
  float          kf_weight_;

  unsigned int   nvtx_;
  float          met_;
  float          metPhi_;
  float          sumEt_;
  float          metSig_;

  int            nphotons_;
  float          phoEt_;
  float          phoEta_;
  float          phoPhi_;

  float          phoMetDeltaPhi_;

  unsigned int   njets_;
  float          leadJetPt_;
  float          leadJetEta_;
  float          leadJetPhi_;
  float          trailJetPt_;
  float          trailJetEta_;
  float          trailJetPhi_;
  
  //for MET syst
  unsigned int   nalljets_;
  float  jet1Pt_;
  float  jet1Eta_;
  float  jet2Pt_;
  float  jet2Eta_;
  float  jet3Pt_;
  float  jet3Eta_;
  float  jet4Pt_;
  float  jet4Eta_;

  //for Lep syst
  unsigned int   nlep_;
  float  lep1Pt_;
  float  lep1Eta_;
  int    lep1Id_;
  float  lep2Pt_;
  float  lep2Eta_;
  int    lep2Id_;
  float  lep3Pt_;
  float  lep3Eta_;
  int    lep3Id_;

  unsigned int   ncosmics_;

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
    tree_->Branch("pu_weightup"       , &pu_weightup_       ,   "pu_weightup/F");
    tree_->Branch("pu_weightdo"       , &pu_weightdo_       ,   "pu_weightdo/F");
    tree_->Branch("kf_weight"         , &kf_weight_         ,   "kf_weight/F");
  
    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("met"          , &met_          ,   "met/F");
    tree_->Branch("metPhi"       , &metPhi_       ,   "metPhi/F");
    tree_->Branch("sumEt"        , &sumEt_        ,   "sumEt/F");
    tree_->Branch("metSig"       , &metSig_       ,   "metSig/F");
  
    tree_->Branch("nphotons"      , &nphotons_      ,   "nphotons/i");
    tree_->Branch("phoEt"         , &phoEt_         ,   "phoEt/F");
    tree_->Branch("phoEta"        , &phoEta_        ,   "phoEta/F");
    tree_->Branch("phoPhi"        , &phoPhi_        ,   "phoPhi/F");
  
    tree_->Branch("phoMetDeltaPhi", &phoMetDeltaPhi_,   "phoMetDeltaPhi/F");
  
    tree_->Branch("njets"         , &njets_         ,   "njets/i");
    tree_->Branch("leadJetPt"         , &leadJetPt_         ,   "leadJetPt/F");
    tree_->Branch("leadJetEta"        , &leadJetEta_        ,   "leadJetEta/F");
    tree_->Branch("leadJetPhi"        , &leadJetPhi_        ,   "leadJetPhi/F");
    tree_->Branch("trailJetPt"        , &trailJetPt_        ,   "trailJetPt/F");
    tree_->Branch("trailJetEta"       , &trailJetEta_       ,   "trailJetEta/F");
    tree_->Branch("trailJetPhi"       , &trailJetPhi_       ,   "trailJetPhi/F");

    tree_->Branch("nalljets"          , &nalljets_          ,   "nalljets/i");
    tree_->Branch("jet1Pt"            , &jet1Pt_            ,   "jet1Pt/F");
    tree_->Branch("jet1Eta"           , &jet1Eta_           ,   "jet1Eta/F");
    tree_->Branch("jet2Pt"            , &jet2Pt_            ,   "jet2Pt/F");
    tree_->Branch("jet2Eta"           , &jet2Eta_           ,   "jet2Eta/F");
    tree_->Branch("jet3Pt"            , &jet3Pt_            ,   "jet3Pt/F");
    tree_->Branch("jet3Eta"           , &jet3Eta_           ,   "jet3Eta/F");
    tree_->Branch("jet4Pt"            , &jet4Pt_            ,   "jet4Pt/F");
    tree_->Branch("jet4Eta"           , &jet4Eta_           ,   "jet4Eta/F");

    tree_->Branch("nlep"              , &nlep_              ,   "nlep/i");
    tree_->Branch("lep1Pt"            , &lep1Pt_            ,   "lep1Pt/F");
    tree_->Branch("lep1Eta"           , &lep1Eta_           ,   "lep1Eta/F");
    tree_->Branch("lep1Id"            , &lep1Id_            ,   "lep1Id_/i");
    tree_->Branch("lep2Pt"            , &lep2Pt_            ,   "lep2Pt/F");
    tree_->Branch("lep2Eta"           , &lep2Eta_           ,   "lep2Eta/F");
    tree_->Branch("lep2Id"            , &lep2Id_            ,   "lep2Id_/i");
    tree_->Branch("lep3Pt"            , &lep3Pt_            ,   "lep3Pt/F");
    tree_->Branch("lep3Eta"           , &lep3Eta_           ,   "lep3Eta/F");
    tree_->Branch("lep3Id"            , &lep3Id_            ,   "lep3Id_/i");

    tree_->Branch("ncosmics"        , &ncosmics_        ,   "ncosmics/i");

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
    tree_->SetBranchAddress("pu_weightup"       , &pu_weightup_  );
    tree_->SetBranchAddress("pu_weightdo"       , &pu_weightdo_  );
    tree_->SetBranchAddress("kf_weight"         , &kf_weight_ );
  
    tree_->SetBranchAddress("nvtx"         , &nvtx_      );
    tree_->SetBranchAddress("met"          , &met_       );
    tree_->SetBranchAddress("metPhi"       , &metPhi_    );
    tree_->SetBranchAddress("sumEt"        , &sumEt_     );
    tree_->SetBranchAddress("metSig"       , &metSig_    );
  
    tree_->SetBranchAddress("nphotons"      , &nphotons_);
    tree_->SetBranchAddress("phoEt"         , &phoEt_     );
    tree_->SetBranchAddress("phoEta"        , &phoEta_    );
    tree_->SetBranchAddress("phoPhi"        , &phoPhi_    );
  
    tree_->SetBranchAddress("phoMetDeltaPhi", &phoMetDeltaPhi_ );
  
    tree_->SetBranchAddress("njets"         , &njets_    );
    tree_->SetBranchAddress("leadJetPt"         , &leadJetPt_    );
    tree_->SetBranchAddress("leadJetEta"        , &leadJetEta_   );
    tree_->SetBranchAddress("leadJetPhi"        , &leadJetPhi_   );
    tree_->SetBranchAddress("trailJetPt"        , &trailJetPt_    );
    tree_->SetBranchAddress("trailJetEta"       , &trailJetEta_   );
    tree_->SetBranchAddress("trailJetPhi"       , &trailJetPhi_   );

    tree_->SetBranchAddress("nalljets"          , &nalljets_          );
    tree_->SetBranchAddress("jet1Pt"            , &jet1Pt_            );
    tree_->SetBranchAddress("jet1Eta"           , &jet1Eta_           );
    tree_->SetBranchAddress("jet2Pt"            , &jet2Pt_            );
    tree_->SetBranchAddress("jet2Eta"           , &jet2Eta_           );
    tree_->SetBranchAddress("jet3Pt"            , &jet3Pt_            );
    tree_->SetBranchAddress("jet3Eta"           , &jet3Eta_           );
    tree_->SetBranchAddress("jet4Pt"            , &jet4Pt_            );
    tree_->SetBranchAddress("jet4Eta"           , &jet4Eta_           );

    tree_->SetBranchAddress("nlep"              , &nlep_            );
    tree_->SetBranchAddress("lep1Pt"            , &lep1Pt_          );
    tree_->SetBranchAddress("lep1Eta"           , &lep1Eta_         );
    tree_->SetBranchAddress("lep1Id"            , &lep1Id_          );
    tree_->SetBranchAddress("lep2Pt"            , &lep2Pt_          );
    tree_->SetBranchAddress("lep2Eta"           , &lep2Eta_         );
    tree_->SetBranchAddress("lep2Id"            , &lep2Id_          );
    tree_->SetBranchAddress("lep3Pt"            , &lep3Pt_          );
    tree_->SetBranchAddress("lep3Eta"           , &lep3Eta_         );
    tree_->SetBranchAddress("lep3Id"            , &lep3Id_          );

    tree_->SetBranchAddress("ncosmics",       &ncosmics_);

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
  pu_weightup_ = -1.;
  pu_weightdo_ = -1.;
  kf_weight_  = -1.;

  nvtx_  = 0;
  met_   = -1.;
  metPhi_= -100.;
  sumEt_ = -1.;
  metSig_= -1.;

  nphotons_ = 0;
  phoEt_ = -1.;
  phoEta_= -100.;
  phoPhi_= -100.;

  phoMetDeltaPhi_ = -100.;

  njets_ = 0;
  leadJetPt_ = -1.;
  leadJetEta_= -100.;
  leadJetPhi_= -100.;
  trailJetPt_ = -1.;
  trailJetEta_= -100.;
  trailJetPhi_= -100.;
  
  njets_ = 0;
  jet1Pt_ = -1.;
  jet1Eta_ = -100.;
  jet2Pt_ = -1.;
  jet2Eta_ = -100.;
  jet3Pt_ = -1.;
  jet3Eta_ = -100.;
  jet4Pt_ = -1.;
  jet4Eta_ = -100.;

  nlep_ = 0;
  lep1Pt_ = -1.;
  lep1Eta_ = -100.;
  lep1Id_ = 0;
  lep2Pt_ = -1.;
  lep2Eta_ = -100.;
  lep2Id_ = 0;
  lep3Pt_ = -1.;
  lep3Eta_ = -100.;
  lep3Id_ = 0;

  ncosmics_ = 0;

}

#endif
