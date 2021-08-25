/*
  //this is test development macor - not used for anything yet

  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  .L $fastMCKalman/fastMCKalman/gasProp/GetVar.cxx+g
  GetVar::Init();
   GetVar::initTree()
  GetVar::makeGasFunction(2000);
   .L $fastMCKalman/fastMCKalman/gasProp/GetGasVariableTest.C+
makeGasTree(10000);

 */
#include "GetVar.h"
#include "TF2.h"

void makeTestTree(){
  GetVar::initTree();
  GetVar::makeGasFunction(2000);
  GetVar::tree->Draw("normdT/normdTFit_0:E/p:fquencher",GetVar::mapId[0]+"&&E/p<1&&E/p>0.005", "colz");
}
void makeTestF2(){
  TF2 f2("f2","GetVar::getDtNorm(0,x,y)",0.0,100,0.005,1);
  f2.Draw("colz");
}

void makeGasTree(Int_t nPoints){
  TTreeSRedirector * pcstream = new TTreeSRedirector("gasTree.root","recreate");
  const float pMin=1000;
  const float pMax=10000;
  const float eMin=100;
  const float eMax=1000;
  const float length=1; // length in cm
  for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
    float p=pMin+gRandom->Rndm()*(pMax-pMin);
    float E=eMin+gRandom->Rndm()*(eMax-eMin);
    int   id=gRandom->Rndm()*GetVar::mapParam.size();
    float quencher=gRandom->Rndm()*100;
    // GetNumberOfX0(int fquencher=10, int p=1000 , int nC=1, int nH=4, int Zgas=18, int Agas=40, float length=0);
    float x0=GetVar::GetNumberOfX0(quencher, p, GetVar::mapParam[id][0], GetVar::mapParam[id][1], GetVar::mapParam[id][2],GetVar::mapParam[id][3],length);
    float rho=GetVar::GetDensity(quencher, p, GetVar::mapParam[id][0], GetVar::mapParam[id][1], GetVar::mapParam[id][2],GetVar::mapParam[id][3]);
    float sigmaDNorm=GetVar::getDtNorm(id,quencher,E/p);
    float sigmaD=GetVar::getDtNorm(id,quencher,E/p)/TMath::Sqrt(p);
    (*pcstream)<<"gasTree"<<
      "p="<<p<<
      "E="<<E<<
      "id="<<id<<
      "nC="<<GetVar::mapParam[id][0]<<
      "nH="<<GetVar::mapParam[id][1]<<
      "Z="<<GetVar::mapParam[id][2]<<
      "A="<<GetVar::mapParam[id][3]<<
      "quencher="<<quencher<<
      "x0="<<x0<<
      "rho="<<rho<<
      "sigmaDNorm="<<sigmaDNorm<<
      "sigmaD="<<sigmaD<<              //micron/sqrt(cm)
      "\n";
  }
  delete pcstream;
}