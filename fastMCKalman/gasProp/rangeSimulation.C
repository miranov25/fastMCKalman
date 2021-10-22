/*
  this is range simulation macro - for th moment does not work - to check with the ATIMA

  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  .L $fastMCKalman/fastMCKalman/gasProp/rangeSimulation.C

  rangeSimul(10000,100);

 */

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"
#include "AliExternalTrackParam.h"
#include "AliNDLocalRegression.h"
TTree * tree =0;
//
//const Float_t  xx0=7.8350968e-05;    //  radiation length Ar/CO2 per cm
//const Float_t xrho=0.0016265266;     //  density  (q/cm&2) per cm

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
    float mass = gRandom->Rndm()*2;    /// generate random mass
    float Z    = gRandom->Rndm()*4;    /// generate random Z
    int pidCode= (gRandom->Rndm()<0.2) ? -1: (gRandom->Rndm()*9.);
    if (pidCode>=0){
      mass=AliPID::ParticleMass(pidCode);
      Z=AliPID::ParticleCharge(pidCode);
    }
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
      "pidCode="<<pidCode<<
      "xx0="<<xTimesRho<<
      "mass="<<mass<<
      "Z="<<Z<<
      "pArray.="<<&pArray<<
      "sigmapArray.="<<&sigmapArray<<
      "dEdxArray.="<<&dEdxArray<<
      "\n";
  }
  delete pcstream;
}

void initTree(){
  TFile * f= TFile::Open("rangeSimul.root");
  tree  = (TTree*)f->Get("range");

}