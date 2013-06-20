//--------------------------------------------------------------------------------------------------
// Perform a plot task using a specified set of samples. Nice and clean.
//
// Authors: C.Paus                                                                        (Aug 2010)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Plot/interface/PlotTask.h"
#endif

using namespace std;
using namespace mithep;

void plot();

void plot(const char *name, const char* title, int logy,
          double xmin, double xmax, double ymin, double ymax,
          int nRebin, double lumi, TString draw="", TString cut="", int nbins=100);

//==================================================================================================
void plotsMonoPhoton(double lumi = 4500.0)
{
  // setup graphics stuff before starting
  MitStyle::Init();
  gROOT->LoadMacro("$CMSSW_BASE/src/MitMonoPhoton/macros/plot.C+");

  // plot from TTree named hHggNtuple
  TString nTuple = "hMonoPhotonTree";
  TString cut    = "pfmet > 0 && nPhotons > 0";

  printf("\n Plotting from tuple:  > %s <  with cuts: > %s \n\n",nTuple.Data(),cut.Data());

  plot(nTuple.Data(),"pfmet","MET [GeV]",0, 80.,220., 0., -1.,1,lumi,"pfmet",cut,140);

  // plot from TH1D
  plot("hPhotonEt","hPhotonEt","photon E_T [GeV]",  0, 0.,  0., 0., -1.,1,lumi);

  return;
}
