//This macro creates the TGraphs necessary to account for
//the photon pt ele->photon fake rate

void makeEleFake(){
  
  const int nBins = 1;

  //Stuff for all
  float theXVal[nBins] = 
  {400};
  float theXErr[nBins] = 
  {1000};

  //Stuff for EleFake
  float theYVal[nBins] = 
  {0.827};
  float theYErr[nBins] = 
  {0.024};
  TGraphErrors* EleFake = new TGraphErrors(nBins);
  EleFake -> SetName("EleFake");
  
  for (int iPoint = 0; iPoint < nBins; iPoint++) {
    
    float thisVal = (1-theYVal[iPoint])/theYVal[iPoint];
    float thisError = sqrt(2)*theYErr[iPoint]/theYVal[iPoint]*thisVal;
    EleFake -> SetPoint(iPoint, theXVal[iPoint],thisVal);
    EleFake -> SetPointError(iPoint, theXErr[iPoint], thisError);

  }

  TFile* outFile = new TFile("EleFake.root","RECREATE");
  EleFake -> Write();
  
}
