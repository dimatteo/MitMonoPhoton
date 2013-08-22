//--------------------------------------------------------------------------------------------------
// $Id: TreeReducer.h,v 1.5 2012/03/28 12:10:21 paus Exp $
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
    void                 SetPuTarget(const TH1D *h) { fPuTarget = h; }
    void                 SetOutput(TFile *f)  { fOutFile = f; }
    void                 SetVerbose(bool b) { fVerbose = b; }
    void                 SetInputBaseDir(const TString s)  { fInputBaseDir = s; }
    void                 SetLumi(float l)  { fLumi = l; }
            
  private:
    // Members
    const Sample * fSample;   // sample to be reduced
    const TH1D   * fPuTarget; // target PU histo
    TFile        * fOutFile;  // target out file

    bool           fVerbose;      // print out what the reducer is doing
    TString        fInputBaseDir; // set the input samples base dir
    float          fLumi;         // set the target lumi

    // Auxiliary functions
    static const TH1D   *sPuWeights;
    float                PuWeight(Float_t npu);    // PU reweighting function
    bool                 EventIsSelected(MitGPTree &tree, int treeType);    // Selection function
    bool                 PhotonIsSelected(float R9, float HoverE, float Iso1, float Iso2, float Iso3);    // Photon id function
    float                GetCorrDeltaPhi(float phi1, float phi2);    // Corr delta phi calculator

    ClassDef(TreeReducer, 0)  // TreeReducer with various options
  };
}
#endif
