/*
   .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    testPCStream(10000,kTRUE);  //setup for the looper development
    testAlice(10000,kTRUE);   // ALICE setup
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
TChain * treeUnit0=0;
void testDummy(){
  double r[]={0,0,0};
  double p[]={1,0,0};
  //particle.simulateParticle(geom, r,p,211,200,90);
}

/// test for looper development with continous tracking
/// \param nParticles
/// \param dumpStream
void testPCStream(Int_t nParticles, bool dumpStream){
    const Int_t   nLayerTPC=250;
      const Float_t kMinPt=0.02;
  const Float_t kMax1Pt=1./100.;
  const Float_t smearR=200;
  const Float_t smearZ=200;
    const Float_t  xx0=7.8350968e-05;
  const Float_t xrho=0.0016265266;
  TStopwatch timer;
  timer.Start();
  fastGeometry geom(nLayerTPC+1);
  geom.fBz=5;
  fastParticle particle(nLayerTPC+1);

  float resol[2]={0.001,0.001};
  resol[0]=0.1;
  resol[1]=0.1;
  geom.setLayerRadiusPower(0,nLayerTPC,1,nLayerTPC,1.0,xx0,xrho,resol);

  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticle.root","recreate");
  particle.fgStreamer=pcstream;
  TTree * tree = 0;
  for (Int_t i=0; i<nParticles; i++){
    double r[]     = {0,0,0};
    Bool_t  isSecondary=gRandom->Rndm()<0.5;
    isSecondary=kFALSE;
    if (isSecondary){
        r[0]=2*(gRandom->Rndm()-0.5)*smearR;
        r[1]=2*(gRandom->Rndm()-0.5)*smearR;
        r[2]=2*(gRandom->Rndm()-0.5)*smearZ;
    }
    double pt      = kMinPt/(kMax1Pt*kMinPt+gRandom->Rndm());
    double phi     = gRandom->Rndm()*TMath::TwoPi();
    double theta = (gRandom->Rndm()-0.5)*3;
    double p[]={pt*sin(phi),pt*cos(phi),pt*theta};
    int    pidCode=int(gRandom->Rndm()*5);
    float  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    int    pdgCode = AliPID::ParticleCode(pidCode)*charge;
    particle.simulateParticle(geom, r,p,pdgCode, 250,nLayerTPC);
    particle.reconstructParticle(geom,pdgCode,nLayerTPC);
    particle.reconstructParticleRotate0(geom,pdgCode,nLayerTPC);
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

/// testAlice configuration
/// \param nParticles
/// \param dumpStream
void testAlice(Int_t nParticles, bool dumpStream){
  const Float_t smearR=10;
  const Float_t smearZ=10;
  const Float_t resolY=0.1;
  const Float_t resolZ=0.1;
  const Float_t bz=5;
  const Int_t   nLayerTPC=160;
  const Float_t kMinPt=0.02;
  const Float_t kMax1Pt=1./100.;
  //const Float_t kMax1Pt=1./0.1;
  const Float_t kThetaMax=2.;
  const Float_t kMaxLength=300;
  // values as obtaied from the gromManager see -  https://alice.its.cern.ch/jira/browse/PWGPP-613?focusedCommentId=263181&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-263181
  const Float_t  xx0=7.8350968e-05;
  const Float_t xrho=0.0016265266;
  //
  TStopwatch timer;
  timer.Start();
  fastGeometry geom= fastGeometry(nLayerTPC);
  geom.fBz=bz;
  fastParticle particle=fastParticle(nLayerTPC);
  // ITS emulation  0 6 layers wit  O(0.001 cm) resolution
  //     X0 is radiation length per layer and rho is the integrated density per layer - obtained from Geant simulation
  float resol[2]={0.001,0.001};
  geom.setLayerRadiusPower(0,5, 4,37, 1.0,0.006,0.16, resol);
  // TPC emulation
  //    X0 is radiation length per layer and rho is the integrated density per layer - obtained from Geant simulation
  resol[0]=0.1;
  resol[1]=0.1;
  //geom.setLayerRadiusPower(6,nLayerTPC,89,260,1.0,0.000025,0.0009,resol);
  geom.setLayerRadiusPower(6,nLayerTPC,89,260,1.0,xx0,xrho,resol);
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticleALICE.root","recreate");
  particle.fgStreamer=pcstream;
  TTree * tree = 0;
  for (Int_t i=0; i<nParticles; i++){
    double r[]     = {0,0,0};
    Bool_t  isSecondary=gRandom->Rndm()<0.5;
    isSecondary=kFALSE;
    if (isSecondary){
        r[0]=2*(gRandom->Rndm()-0.5)*smearR;
        r[1]=2*(gRandom->Rndm()-0.5)*smearR;
        r[2]=2*(gRandom->Rndm()-0.5)*smearZ;
    }
    double pt      = kMinPt/(kMax1Pt*kMinPt+gRandom->Rndm());
    double phi     = gRandom->Rndm()*TMath::TwoPi();
    double theta = 2.*(gRandom->Rndm()-0.5)*kThetaMax;
    double p[]={pt*sin(phi),pt*cos(phi),pt*theta};
    int    pidCode=int(gRandom->Rndm()*5);                   // PID code of particles - 0-electron ,1-muon, 2-pion, 3-Kaon
    //pidCode=2;
    float  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    int    pdgCode = AliPID::ParticleCode(pidCode)*charge;  // PID code covnerted to the PdgCode
    particle.simulateParticle(geom, r,p,pdgCode, kMaxLength,nLayerTPC);
    particle.reconstructParticle(geom,pdgCode,nLayerTPC);
    particle.reconstructParticleRotate0(geom,pdgCode,nLayerTPC);
    Float_t mass =AliPID::ParticleMass(pidCode);

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


/// testAlice configuration
/// \param nParticles
/// \param dumpStream
void testAlice3(Int_t nParticles, bool dumpStream){
  const Float_t smearR=10;
  const Float_t smearZ=10;
  const Float_t resolY=0.1;
  const Float_t resolZ=0.1;
  const Float_t bz=5;
  const Int_t   nLayerTPC=160;
  const Float_t kMinPt=0.02;
  const Float_t kMax1Pt=1./100.;
  //const Float_t kMax1Pt=1./0.2;
  const Float_t kThetaMax=2.;
  const Float_t kMaxLength=300;
  // values as obtaied from the gromManager see -  https://alice.its.cern.ch/jira/browse/PWGPP-613?focusedCommentId=263181&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-263181
  const Float_t  xx0=7.8350968e-05;
  const Float_t xrho=0.0016265266;
  //
  TStopwatch timer;
  timer.Start();
  fastGeometry geom= fastGeometry(nLayerTPC);
  geom.fBz=bz;
  fastParticle particle=fastParticle(nLayerTPC);
  // ITS emulation  0 6 layers wit  O(0.001 cm) resolution
  //     X0 is radiation length per layer and rho is the integrated density per layer - obtained from Geant simulation
  float resol[2]={0.001,0.001};
  geom.setLayerRadiusPower(0,5, 4,37, 1.0,0.006,0.16, resol);
  // TPC emulation
  //    X0 is radiation length per layer and rho is the integrated density per layer - obtained from Geant simulation
  resol[0]=0.1;
  resol[1]=0.1;
  //geom.setLayerRadiusPower(6,nLayerTPC,89,260,1.0,0.000025,0.0009,resol);
  geom.setLayerRadiusPower(6,nLayerTPC,89,260,1.0,xx0,xrho,resol);
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticleALICE.root","recreate");
  particle.fgStreamer=pcstream;
  TTree * tree = 0;
  for (Int_t i=0; i<nParticles; i++){
    double r[]     = {0,0,0};
    Bool_t  isSecondary=gRandom->Rndm()<0.5;
    isSecondary=kFALSE;
    if (isSecondary){
        r[0]=2*(gRandom->Rndm()-0.5)*smearR;
        r[1]=2*(gRandom->Rndm()-0.5)*smearR;
        r[2]=2*(gRandom->Rndm()-0.5)*smearZ;
    }
    double pt      = kMinPt/(kMax1Pt*kMinPt+gRandom->Rndm());
    double phi     = gRandom->Rndm()*TMath::TwoPi();
    double theta = 2.*(gRandom->Rndm()-0.5)*kThetaMax;
    double p[]={pt*sin(phi),pt*cos(phi),pt*theta};
    int    pidCode=int(gRandom->Rndm()*5);                   // PID code of particles - 0-electron ,1-muon, 2-pion, 3-Kaon
    float  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    int    pdgCode = AliPID::ParticleCode(pidCode)*charge;  // PID code covnerted to the PdgCode
    particle.simulateParticle(geom, r,p,pdgCode, kMaxLength,nLayerTPC);
    particle.reconstructParticle(geom,pdgCode,nLayerTPC);
    particle.reconstructParticleRotate0(geom,pdgCode,nLayerTPC);
    Float_t mass =AliPID::ParticleMass(pidCode);

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
  fastParticle::setAliases(tree);
//  tree.SetAlias("gxIn","cos(part.fParamIn[].fAlpha)*part.fParamIn[].fX");
//  tree.SetAlias("gyIn","sin(part.fParamIn[].fAlpha)*part.fParamIn[].fX");
//  tree.SetAlias("gxMC","cos(part.fParamMC[].fAlpha)*part.fParamMC[].fX");
//  tree.SetAlias("gyMC","sin(part.fParamMC[].fAlpha)*part.fParamMC[].fX");
//  tree.SetAlias("gzMC","part.fParamMC[].fP[1]");
//  tree.SetAlias("rMC","part.fParamMC[].fX");
//  tree.SetAlias("ptMC","part.fParamMC[0].fData.Pt()");
//  tree.SetAlias("pMC","part.fParamMC[0].fData.P()");
//  tree.SetAlias("ptRec","part.fParamIn[0].fData.Pt()");
//  //
//  tree.SetAlias("layer","Iteration$");
//  tree.SetAlias("c0MC","sqrt(part.fParamMC[].fC[0])");
//  tree.SetAlias("c2MC","sqrt(part.fParamMC[].fC[2])");
//  tree.SetAlias("c14MC","sqrt(part.fParamMC[].fC[14])");
//  tree.SetAlias("c0In","sqrt(part.fParamIn[].fC[0])");
//  tree.SetAlias("c2In","sqrt(part.fParamIn[].fC[2])");
//  tree.SetAlias("c0InRot","sqrt(part.fParamInRot[].fC[0])");
//  tree.SetAlias("c2InRot","sqrt(part.fParamInRot[].fC[2])");
//  tree.SetAlias("dEdxExp","AliExternalTrackParam::BetheBlochAleph(pMC/AliPID::ParticleMass(pidCode))");
//  tree.SetAlias("dEdxExpSolid","AliExternalTrackParam::BetheBlochSolid(pMC/AliPID::ParticleMass(pidCode))");
//  tree.SetAlias("dEdxExpSolidL","AliExternalTrackParam::BetheBlochSolid(part.fParamMC[].fData.P()/AliPID::ParticleMass(pidCode))");
//  tree.SetAlias("dEdxExpSolidL1","AliExternalTrackParam::BetheBlochSolid(part.fParamMC[Iteration$-1].fData.P()/AliPID::ParticleMass(pidCode))");
//  tree.SetAlias("elossTPCIn","(part.fParamIn[159].fData.GetP()-part.fParamIn[7].fData.GetP())/part.fParamMC[1].fData.GetP()");
//  tree.SetAlias("elossTPCMC","(part.fParamMC[159].fData.GetP()-part.fParamMC[7].fData.GetP())/part.fParamMC[1].fData.GetP()");
}

void initTreeFast(){
  treeFast  = AliXRDPROOFtoolkit::MakeChainRandom("fastParticle.list","fastPart",0,10000);
  treeTurn  = AliXRDPROOFtoolkit::MakeChainRandom("fastParticle.list","turn",0,10000);
  treeUnit0  = AliXRDPROOFtoolkit::MakeChainRandom("fastParticle.list","UnitTestDumpCorrectForMaterial",0,10000);
  treeFast->SetMarkerStyle(21);
  treeFast->SetMarkerSize(0.5);
   treeUnit0->SetMarkerStyle(21);
  treeUnit0->SetMarkerSize(0.5);

  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  setAliases(*treeFast);
  //
  treeUnit0->SetAlias("dEdxOutIn","AliExternalTrackParam::BetheBlochSolid(0+paramStepRK.P()/mass)/AliExternalTrackParam::BetheBlochSolid(paramIn.P()/mass)");
  treeUnit0->SetAlias("dEdxIn","AliExternalTrackParam::BetheBlochSolid(paramIn.P()/mass+0)");
  treeUnit0->SetAlias("dEdxOut","AliExternalTrackParam::BetheBlochSolid(paramRK.P()/mass+0)");

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
