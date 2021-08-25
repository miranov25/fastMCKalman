/*
this is test development macor - not used for anything yet

  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  .L $fastMCKalman/fastMCKalman/gasProp/GetVar.cxx+
  .L $fastMCKalman/fastMCKalman/gasProp/testInterpolationGas.C+
  TFile f("all_merged.root");
  tree = (TTree*) f.Get("gastree");
  makeGasFunction(2000);

 */

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"
#include "AliExternalTrackParam.h"
#include "AliNDLocalRegression.h"
#include "GetVar.h"



TTree * tree = 0;
map<string,AliNDLocalRegression*> mapRegressionDtStr;
map<int,TString> mapId;
map<int,AliNDLocalRegression*> mapRegressionDtInt;

double getDt(int id,  float fquencher, float ep){
  AliNDLocalRegression * regresion=mapRegressionDtInt[id];
  if (regresion==0) return 0;
  Double_t inData[2]={fquencher,ep};
  return regresion->Eval(inData);
}
char * getIdName(int id){return (char*) mapId[id].Data();}
void initTree(){
  TFile *f = TFile::Open("all_merged.root");
  tree = (TTree*) f->Get("gastree");
}
void testInterpolation(){
  TFile f("all_merged.root");
  tree = (TTree*) f.Get("gastree");
  tree->Draw("sqrt(p)*dt:(fquencher):E/p","sim==1&&nH==4&&nC==1&&Agas==40&&Zgas==18&&B==0&&E<5000&&E/p>0.2","colz");
  //::Info("makeElectronFitRegion","BEGIN");
  //timerFunction.Start();
  TString range0="(20,0,100,20,0,1)";
  TString kernel0="5:0.01+(E/p)*0.05";
  TString valueString=TString::Format("normdT:0.1");
  TString fitString=TString::Format("fquencher:E/p");
  const char *fitSelection="sim==1&&nH==4&&nC==1&&Agas==40&&Zgas==18&&B==0";
  //
  tree->SetAlias("normdT","sqrt(p)*dt");
  tree->SetAlias("fitSelection",fitSelection);
  AliNDLocalRegression *fit = AliNDLocalRegression::MakeRegression(tree,"normdTFit",range0, valueString.Data(), fitString,fitSelection, kernel0,0.001);
  //
}



void makeGasFunction(int minEntries){
  int nH[]={1,4,6,8,10};
  int nC[]={1,2,3,4};
  int Agas[]={4,40};
  int Zgas[]={2,18};
  float B[]={0,0.5};
  //
  tree->SetAlias("normdT","sqrt(p)*dt");
  Int_t id=0;
  for (int iH=0; iH<5; iH++)
    for (int iC=0; iC<4; iC++)
      for (int iAgas=0; iAgas<2;iAgas++)
        for (int iZgas=0; iZgas<2;iZgas++)
          for (int iB=0; iB<2; iB++){
            TString selection=TString::Format("sim==1&&nH==%d&&nC==%d&&Agas==%d&&Zgas==%d&&B==%0.2f",nH[iH],nC[iC], Agas[iAgas],Zgas[iZgas],B[iB]);
            Int_t entries=tree->Draw("fquencher:E/p",selection,"goff");
            if (entries<minEntries) continue;
            printf("%s\t%d\n",selection.Data(),entries);
            TString range0="(20,0,100,40,0.005,1)";
            TString kernel0="5:0.01+(E/p)*0.05";
            TString valueString=TString::Format("normdT:0.1");
            TString fitString=TString::Format("fquencher:E/p");
            //
            AliNDLocalRegression *fit = AliNDLocalRegression::MakeRegression(tree,Form("normdTFit_%d",id),range0, valueString.Data(), fitString,selection+"&&E/p<1&&E/p>0.005", kernel0,0.001);
            mapId[id]=selection.Data();
            mapRegressionDtStr[selection.Data()]=fit;
            mapRegressionDtInt[id]=fit;
            id++;
          }
}

void rangeSimul(int nPoints=10000, Int_t nSteps=100){
  // Todo:
  //  - add landau smearing
  //  - find the landau width
  //  - compare with atima
  const float  kLandau=1;                         // landau factor
  const Float_t E0=2.e-8;                         // E0 in GeV - 20 eV
  const float bgStart=0.015;
  float  xTimesRho=7.8350968e-05;                /// for neon https://www.sciencedirect.com/science/article/pii/S0168900210008910 0.6 % of X/X0, Argon, CO2 obtained from the Geant
  TVectorF pArray(nSteps);
  TVectorF sigmapArray(nSteps);
  TVectorF dEdxArray(nSteps);
  TTreeSRedirector * pcstream = new TTreeSRedirector("rangeSimul.root","recreate");
  //
  for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
    float mass = gRandom->Rndm();    /// generate random mass
    float Z    = gRandom->Rndm()*4;    /// generate random Z
    float bg=bgStart;
    float p=bg*mass;
    float sigmaP=0;
    for (int iStep=0; iStep<nSteps; iStep++){
      //Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction);
      Double_t dEdx=Z*Z*AliExternalTrackParam::BetheBlochGeant(bg);
      Double_t dE=(dEdx*xTimesRho);
      Double_t sigmadERel=kLandau/TMath::Sqrt(dE/E0);
      pArray[iStep]=p;
      sigmapArray[iStep]=sigmaP;
      dEdxArray[iStep]=dEdx;
      Double_t dPdx=dEdx*xTimesRho*TMath::Sqrt(1.+1./(bg*bg));
      sigmaP=TMath::Sqrt(sigmaP*sigmaP +dPdx*dPdx*sigmadERel*sigmadERel);    // not sure about error propagation
      p+=dPdx;
      bg=p/mass;
    }
    (*pcstream)<<"range"<<
      "mass="<<mass<<
      "Z="<<Z<<
      "pArray.="<<&pArray<<
      "sigmapArray.="<<&sigmapArray<<
      "dEdxArray.="<<&dEdxArray<<
      "\n";
  }
  delete pcstream;
}
