//This macro creates the TGraphs necessary to account for
//the photon pt dependent k-factors for Zgamma and Wgamma
//source of calculation is
//https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=283118

void makeKFactors(){
  
  const int nBins = 7;

  //Stuff for all
  float theXVal[nBins] = 
  {152.5, 175, 220, 325, 550, 850};
  float theXErr[nBins] = 
  {7.5, 15, 30, 75, 150, 150};

  //Stuff for ZGamma
  float ZGammaVal[nBins] = 
  {1.721, 1.615, 1.565, 1.548, 1.605, 1.839};
  float ZGammaStat[nBins] = 
  {0.031, 0.023, 0.009, 0.014, 0.014, 0.048};
  float ZGammaSyst[nBins] = 
  {0.116, 0.122, 0.156, 0.168, 0.230, 0.578};
  TGraphErrors* ZGamma = new TGraphErrors(nBins);
  ZGamma -> SetName("ZGamma");
  
  //Stuff for WGamma
  float WGammaVal[nBins] = 
  {1.310, 1.267, 1.281, 1.449, 1.534, 1.953};
  float WGammaStat[nBins] = 
  {0.152, 0.069, 0.132, 0.051, 0.031, 0.080};
  float WGammaSyst[nBins] = 
  {0.240, 0.184, 0.196, 0.319, 0.360, 1.104};
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
