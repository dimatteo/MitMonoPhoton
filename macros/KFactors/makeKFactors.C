//This macro creates the TGraphs necessary to account for
//the photon pt dependent k-factors for Zgamma and Wgamma
//source of calculation is
//AN 12/439


void makeKFactors(){
  
  const int nBins = 7;

  //Stuff for all
  float theXVal[nBins] = 
  {152.5, 175, 220, 325, 550, 850};
  float theXErr[nBins] = 
  {7.5, 15, 30, 75, 150, 150};

  //Stuff for ZGamma
  float ZGammaVal[nBins] = 
  {1.320, 1.430, 1.499, 1.528, 1.642, 1.880};
  float ZGammaStat[nBins] = 
  {0.013, 0.008, 0.007, 0.006, 0.010, 0.045};
  float ZGammaSyst[nBins] = 
  {0.108, 0.137, 0.178, 0.194, 0.319, 0.911};
  TGraphErrors* ZGamma = new TGraphErrors(nBins);
  ZGamma -> SetName("ZGamma");
  
  //Stuff for WGamma
  float WGammaVal[nBins] = 
  {1.527, 1.546, 1.642, 1.675, 1.830, 1.982};
  float WGammaStat[nBins] = 
  {0.059, 0.031, 0.036, 0.006, 0.011, 0.057};
  float WGammaSyst[nBins] = 
  {0.251, 0.289, 0.305, 0.383, 0.746, 1.643};
  TGraphErrors* WGamma = new TGraphErrors(nBins);
  WGamma -> SetName("WGamma");

  for (int iPoint = 0; iPoint < nBins; iPoint++) {
    
    float thisZerror = sqrt(ZGammaStat[iPoint]*ZGammaStat[iPoint]+ZGammaSyst[iPoint]*ZGammaSyst[iPoint]);
    ZGamma -> SetPoint(iPoint, theXVal[iPoint], ZGammaVal[iPoint]);
    ZGamma -> SetPointError(iPoint, theXErr[iPoint], thisZerror);

    float thisWerror = sqrt(WGammaStat[iPoint]*WGammaStat[iPoint]+WGammaSyst[iPoint]*WGammaSyst[iPoint]);
    WGamma -> SetPoint(iPoint, theXVal[iPoint], WGammaVal[iPoint]);
    WGamma -> SetPointError(iPoint, theXErr[iPoint], thisWerror);

  }

  TFile* outFile = new TFile("KFactors.root","RECREATE");
  ZGamma -> Write();
  WGamma -> Write();
  
}
