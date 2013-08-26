//--------------------------------------------------------------------------------------------------
// $Id: TreeReducer.h,v 1.2 2013/08/23 22:47:34 dimatteo Exp $
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
    void                 SetOutput(TFile *f)  { fOutFile = f; }
    void                 SetVerbose(bool b) { fVerbose = b; }
    void                 SetInputBaseDir(const TString s)  { fInputBaseDir = s; }
    void                 SetLumi(float l)  { fLumi = l; }
            
  private:
    // Members
    const Sample * fSample;   // sample to be reduced
    const TH1D   * fPUTarget; // target PU histo
    TFile        * fOutFile;  // target out file

    bool           fVerbose;      // print out what the reducer is doing
    TString        fInputBaseDir; // set the input samples base dir
    float          fLumi;         // set the target lumi

    // Auxiliary functions
    static const TH1D   *sPUWeights;
    float                PUWeight(Float_t npu);  // PU reweighting function
    float                KFactorWeight(Float_t scale, Float_t phet);  // KFactor reweighting function
    bool                 EventIsSelected(MitGPTree &tree, int treeType);    // Selection function
    bool                 PhotonIsSelected(float R9, float HoverE, float CovIetaIeta, float Iso1, float Iso2, float Iso3);    // Photon id function
    float                GetCorrDeltaPhi(float phi1, float phi2);    // Corr delta phi calculator

    ClassDef(TreeReducer, 0)  // TreeReducer with various options
  };
}
#endif
