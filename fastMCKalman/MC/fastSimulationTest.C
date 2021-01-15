/*
   .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    testPCStream(5000,kTRUE);
    //initTreeFast()
    .> a.log
     testPCStream(5000,kTRUE);
     .>

 */
#include "fastSimulation.h"
#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"

TChain * treeFast = 0;
TChain * treeTurn=0;
void testDummy(){
  double r[]={0,0,0};
  double p[]={1,0,0};
  //particle.simulateParticle(geom, r,p,211,200,90);
}

void testPCStream(Int_t nParticles, bool dumpStream){
    const Float_t smearR=200;
    const Float_t smearZ=200;
  TStopwatch timer;
  timer.Start();
  fastGeometry geom= fastGeometry(201.);
  geom.fBz=5;
  fastParticle particle=fastParticle(200);

  float resol[2]={0.001,0.001};
  geom.setLayerRadiusPower(0,80, 4,80, 1.2,0.006,0.16,resol);
  resol[0]=0.1;
  resol[1]=0.1;
  geom.setLayerRadiusPower(80,160,80,250,1.0,0.0001,0.08,resol);


  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticle.root","recreate");
  particle.fgStreamer=pcstream;
  TTree * tree = 0;
  for (Int_t i=0; i<nParticles; i++){
    double r[]     = {0,0,0};
    Bool_t  isSecondary=gRandom->Rndm()<0.5;
    //isSecondary=kTRUE;
    if (isSecondary){
        r[0]=2*(gRandom->Rndm()-0.5)*smearR;
        r[1]=2*(gRandom->Rndm()-0.5)*smearR;
        r[2]=2*(gRandom->Rndm()-0.5)*smearZ;
    }
    double pt      = 0.02/(0.002+gRandom->Rndm());
    double phi     = gRandom->Rndm()*TMath::TwoPi();
    double theta = (gRandom->Rndm()-0.5)*3;
    double p[]={pt*sin(phi),pt*cos(phi),pt*theta};
    int    pidCode=int(gRandom->Rndm()*5);
    float  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    int    pdgCode = AliPID::ParticleCode(pidCode)*charge;
    particle.simulateParticle(geom, r,p,pdgCode, 250,161);
    particle.reconstructParticle(geom,pdgCode,160);
    //particle.simulateParticle(geom, r,p,211, 250,161);
    //particle.reconstructParticle(geom,211,160);
    if (dumpStream==kFALSE) continue;
    if (tree) tree->Fill();
    else {
      (*pcstream) << "fastPart" <<
                  "i=" << i <<
                  "isSecondary="<<isSecondary<<
                  "pidCode="<<pidCode<<
                  "charge="<<charge<<
                  "part.=" << &particle <<
                  "\n";
      tree=  ((*pcstream) << "fastPart").GetTree();
    }
  }
  delete pcstream;
  timer.Print();
}

void setAliases(TTree & tree){
  tree.SetAlias("gxIn","cos(part.fParamIn[].fAlpha)*part.fParamIn[].fX");
  tree.SetAlias("gyIn","sin(part.fParamIn[].fAlpha)*part.fParamIn[].fX");
  tree.SetAlias("gxMC","cos(part.fParamMC[].fAlpha)*part.fParamMC[].fX");
  tree.SetAlias("gyMC","sin(part.fParamMC[].fAlpha)*part.fParamMC[].fX");
  tree.SetAlias("gzMC","part.fParamMC[].fP[1]");
  tree.SetAlias("rMC","part.fParamMC[].fX");

  tree.SetAlias("ptMC","part.fParamMC[0].fData.Pt()");
  tree.SetAlias("ptRec","part.fParamIn[0].fData.Pt()");
}

void initTreeFast(){
  treeFast  = AliXRDPROOFtoolkit::MakeChainRandom("fastParticle.list","fastPart",0,10000);
   treeTurn  = AliXRDPROOFtoolkit::MakeChainRandom("fastParticle.list","turn",0,10000);
  treeFast->SetMarkerStyle(21);
  treeFast->SetMarkerSize(0.5);
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  setAliases(*treeFast);
}

void drawDisplay(){
  TCanvas * canvas1 = new TCanvas("c1","c1",600,500);
  gSystem->mkdir("fig/");
  treeFast->Draw("gyMC:gxMC:part.fParamMC[0].fP[4]","","colz",100);
  gPad->SaveAs("fig/gygzMC.png");
  treeFast->Draw("gyIn:gxIn:part.fParamMC[0].fP[4]","","colz",100);
  gPad->SaveAs("fig/gygzIn.png");
}


/*

void testDataFrame(){
  RDataFrame d(100); // a RDF that will generate 100 entries (currently empty)
  int x = -1;
  auto d_with_columns = d.Define("x", [&x] { return ++x; })
                       .Define("xx", [&x] { return x*x; });
  d_with_columns.Snapshot("myNewTree", "newfile.root");
  //
  RDataFrame df("particles", "particle.root");
  //auto zMean = d.Define("z", "sqrt(x*x + y*y)").Mean("z");
}

*/

/*

void checkMaterialAliRoot(){
  Double_t x0[3]={0,0,0};
  Double_t x1[3]={6,0,0};
  Double_t param[10];
  AliTrackerBase::MeanMaterialBudget(x0,x1,param);
  Double_t xrho=param[0]*param[4], xx0=param[1];
  // root [34] xrho
  // (Double_t)4.66295523219595331e-01
  // root [35] xx0
  // (Double_t)1.32765905087294800e-02

}*/
