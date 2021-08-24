/*
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  .L $fastMCKalman/fastMCKalman/gasProp/testInterpolationGas.C+
  TFile f("all_merged.root");
  tree = (TTree*) f.Get("gastree");
  makeGasFunction(2000);
 */

#include "TFile.h"
#include "TTree.h"
#include "AliNDLocalRegression.h"


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
