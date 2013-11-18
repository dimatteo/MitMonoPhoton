//--------------------------------------------------------------------------------------------------
// $Id: TreeReducer.h,v 1.5 2013/10/27 06:24:53 dimatteo Exp $
//
// TreeReducer
//
// This class implements TreeReducer which allows to reduce the trees with various additions.
//
// Authors: L. Di Matteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOPHOTON_UTILS_TREEREDUCER_H
#define MITMONOPHOTON_UTILS_TREEREDUCER_H

#include <TH1D.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitMonoPhoton/Core/MitGPTree.h"

namespace mithep 
{
  class TreeReducer
  {
  public:
    TreeReducer(const Sample *mySample);
    ~TreeReducer();

    // Core method
    void                 MakeTree();

    // Various parameters
    void                 SetPUTarget(const TH1D *h) { fPUTarget = h; }
    void                 SetPUTargetUp(const TH1D *h) { fPUTargetUp = h; }
    void                 SetPUTargetDown(const TH1D *h) { fPUTargetDown = h; }
    void                 SetPhoFakeRate(const TGraphErrors *g) { gPhoFakeRate = g; }
    void                 SetOutput(TFile *f)  { fOutFile = f; }
    void                 SetVerbose(bool b) { fVerbose = b; }
    void                 SetInputBaseDir(const TString s)  { fInputBaseDir = s; }
    void                 SetSelectionMode(const TString s)  { fInputSelection = s; }
    void                 SetLumi(float l)  { fLumi = l; }
            
  private:
    // Members
    const Sample       * fSample;   // sample to be reduced
    const TH1D         * fPUTarget; // target PU histo
    const TH1D         * fPUTargetUp; // target PU histo up
    const TH1D         * fPUTargetDown; // target PU histo down
    const TGraphErrors * gPhoFakeRate; // photon fake rate function
    TFile              * fOutFile;  // target out file

    bool           fVerbose;        // print out what the reducer is doing
    TString        fInputBaseDir;   // set the input samples base dir
    TString        fInputSelection; // set the selection mode
    float          fLumi;           // set the target lumi

    // Auxiliary functions
    static const TH1D   *sPUWeights;
    static const TH1D   *sPUWeightsUp;
    static const TH1D   *sPUWeightsDown;
    float                PUWeight(Float_t npu);  // PU reweighting function
    float                PUWeightUp(Float_t npu);  // PU reweighting function
    float                PUWeightDown(Float_t npu);  // PU reweighting function
    float                KFactorWeight(Float_t scale, Float_t phet);  // KFactor reweighting function
    float                ScaleFactorWeight(Float_t R9, Float_t phet);  // Scale reweighting function
    bool                 EventIsSelected(MitGPTree &tree, int treeType, int& theGoodPhoton);    // Selection function
    bool                 PhotonIsSelected(float R9, float HoverE, float CovIetaIeta, float Iso1, float Iso2, float Iso3);    // Photon id function
    bool                 PhotonIsFake(float R9, float HoverE, float CovIetaIeta, float Iso1, float Iso2, float Iso3);    // Photon fake id function
    float                GetCorrDeltaPhi(float phi1, float phi2);    // Corr delta phi calculator

    ClassDef(TreeReducer, 0)  // TreeReducer with various options
  };
}
#endif
