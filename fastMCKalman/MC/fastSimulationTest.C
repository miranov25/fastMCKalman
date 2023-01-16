/*
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
      gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
   .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    AliPDG::AddParticlesToPdgDataBase();
    testTPC(10000,kTRUE);            //setup for the looper development
    testAlice(10000,kTRUE);          // ALICE setup
    testAlice3Werner(50000,kTRUE)    // ALICE3 setup
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
//#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"
#include "AliPID.h"
const Float_t kDecayFraction=0.5;
const Float_t kRandomPDGFraction=0.5;

TChain * treeFast = 0;
TChain * treeTurn=0;
TChain * treeUnit0=0;
TChain * treeSeed=0;
void testDummy(){
  double r[]={0,0,0};
  double p[]={1,0,0};
  //particle.simulateParticle(geom, r,p,211,200,90);
}

/// test for looper development with continous tracking - ALICE TPC gas cylinder without ITS - emulation of the gas detectors
/// \param nParticles
/// \param dumpStream
void testTPC(Int_t nParticles, bool dumpStream=1){

  const Int_t   nLayerTPC=250;
  const Int_t   nPoints=nLayerTPC*4;     ///maximum number of points a track can have, different from nLayerTPC for Loopers/Secondaries
  const Float_t kMinPt=0.02;
  const Float_t kMax1Pt=1./100.;
  const Float_t kFlatPtMax=50;
  const Float_t kFlatPtFraction=0.3;
  const Float_t smearR=200;
  const Float_t smearZ=200;
  const Float_t  xx0=7.8350968e-05;
  const Float_t  xrho=0.0016265266;
  const Float_t kNominalFraction=0.3;     // fraction of nominal properties
  const Float_t kMaterialScaling=10;      // material random scaling to
  const Float_t kMaxResol=0.2;            // point resolution scan max resolution
  const Float_t kMinResol=0.01;           //  point resolution scan min resolution
  const Float_t kDefResol=0.1;           //  point resolution scan min resolution

  TStopwatch timer;
  timer.Start();
  fastGeometry geom(nLayerTPC+1);
  geom.fBz=5;


  float resol[2]={0.001,0.001};
  resol[0]=0.1;
  resol[1]=0.1;
  geom.setLayerRadiusPower(0,nLayerTPC,1,nLayerTPC,1.0,xx0,xrho,resol);

  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticle.root","recreate");
  TTree * tree = 0;
  for (Int_t i=0; i<nParticles; i++){
    fastParticle particle(nLayerTPC+1);
    particle.fAddMSsmearing=true;
    particle.fAddPadsmearing=true;
    particle.fUseMCInfo=true;
    particle.fgStreamer=pcstream;
    particle.gid=i;
    // generate scan detector properties
    Float_t matScaling  =(gRandom->Rndm()<kNominalFraction) ? 1.:  (gRandom->Rndm()*kMaterialScaling)+1;
    Float_t resolScan=(gRandom->Rndm()<kNominalFraction) ? kDefResol: gRandom->Rndm()*kMaxResol;
    for (size_t iLayer=0; iLayer<geom.fLayerX0.size();iLayer++) {
      geom.fLayerX0[iLayer] = xx0 * matScaling;
      geom.fLayerRho[iLayer] = xrho * matScaling;
      geom.fLayerResolRPhi[iLayer] = resolScan;
      geom.fLayerResolZ[iLayer] = resolScan;
    }
    double r[]     = {0,0,0};
    Bool_t  isSecondary=gRandom->Rndm()<0.5;
    // isSecondary=kFALSE;
    if (isSecondary){
        r[0]=2*(gRandom->Rndm()-0.5)*smearR;
        r[1]=2*(gRandom->Rndm()-0.5)*smearR;
        r[2]=2*(gRandom->Rndm()-0.5)*smearZ;
    }
    double pt      = kMinPt/(kMax1Pt*kMinPt+gRandom->Rndm());
    if (gRandom->Rndm()<kFlatPtFraction) pt= gRandom->Rndm()*kFlatPtMax;
    double phi     = gRandom->Rndm()*TMath::TwoPi();
    double theta = (gRandom->Rndm()-0.5)*3;
    double p[]={pt*sin(phi),pt*cos(phi),pt*theta};
    int    pidCode=int(gRandom->Rndm()*5);          //avoid unrecognized pdg codes
    short  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    int64_t   pdgCode = AliPID::ParticleCode(pidCode)*charge;
    if (gRandom->Rndm()<kRandomPDGFraction) pdgCode=0;
    Bool_t  hasDecay=(gRandom->Rndm()<kDecayFraction);
    Float_t decayLength= hasDecay ?gRandom->Rndm()*geom.fLayerRadius[geom.fLayerRadius.size()-1]:0;
    particle.fDecayLength=decayLength;
    particle.simulateParticle(geom, r,p,pdgCode,nPoints,nPoints);
    particle.reconstructParticle(geom,pdgCode,nPoints);
    fastParticle particle0 = particle;
    particle.reconstructParticleFull(geom,pdgCode,nPoints);
    particle.reconstructParticleFullOut(geom,pdgCode,nPoints);
    particle.refitParticle();
    //particle.reconstructParticleRotate0(geom,pdgCode,nPoints);
    //particle.simulateParticle(geom, r,p,211, 250,161);
    //particle.reconstructParticle(geom,211,160);
    if (dumpStream==kFALSE) continue;
    if (tree) tree->Fill();
    else {
      (*pcstream) << "fastPart" <<
                  "i=" << i <<
                  "densScaling="<<matScaling<<
                  "geom.="<<&geom<<
                  "hasDecay="<<hasDecay<<
                  "isSecondary="<<isSecondary<<
                  "pidCode="<<pidCode<<
                  "pdgCode="<<pdgCode<<
                  "charge="<<charge<<
                  "part.=" << &particle0 <<
                  "partFull.=" << &particle <<
                  "\n";
      tree=  ((*pcstream) << "fastPart").GetTree();
    }
  }
  delete pcstream;
  timer.Print();
}

/// testAlice configuration ITS+TPC with material budget as in the Run1/2
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
    int    pidCode=int(gRandom->Rndm()*8);                   // PID code of particles - 0-electron ,1-muon, 2-pion, 3-Kaon

    //pidCode=2;
    long  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    long    pdgCode = AliPID::ParticleCode(pidCode)*charge;  // PID code covnerted to the PdgCode
    if (gRandom->Rndm()<kRandomPDGFraction) pdgCode=0;
     Float_t decayLength= (gRandom->Rndm()<kDecayFraction) ?gRandom->Rndm()*geom.fLayerRadius[geom.fLayerRadius.size()-1]:0;
    particle.fDecayLength=decayLength;
    particle.simulateParticle(geom, r,p,pdgCode, kMaxLength,nLayerTPC);
    particle.reconstructParticle(geom,pdgCode,nLayerTPC);
    particle.reconstructParticleRotate0(geom,pdgCode,nLayerTPC);

    if (dumpStream==kFALSE) continue;
    if (tree) tree->Fill();
    else {
      (*pcstream) << "fastPart" <<
                  "i=" << i <<
                  "isSecondary="<<isSecondary<<
                  "pidCode="<<pidCode<<
                  "pdgCode="<<pdgCode<<
                  "charge="<<charge<<
                  "part.=" << &particle <<
                  "\n";
      tree=  ((*pcstream) << "fastPart").GetTree();
    }
  }
  delete pcstream;
  timer.Print();
}


/// testAlice 3 configuration as proposed in https://github.com/preghenella/DelphesO2/blob/a058f94f6cb887edcf725fd991d16ca5f7b76e0b/src/lutWrite.werner.cc
/// \param nParticles
/// \param dumpStream
void testAlice3Werner(Int_t nParticles, bool dumpStream){
  // simulation setup parameters
  const Float_t smearR=10;
  const Float_t smearZ=10;

  const Float_t kMinPt=0.02;
  const Float_t kMax1Pt=1./100.;
  //const Float_t kMax1Pt=1./0.1;
  const Float_t kThetaMax=2.;
  const Float_t kMaxLength=300;
  //
  // new ideal Pixel properties?
  const Float_t x0IB     = 0.001;
  const Float_t x0OB     = 0.01;
  const Float_t xrhoIB     = 2.3292e-02; // 100 mum Si
  const Float_t xrhoOB     = 2.3292e-01; // 1000 mum Si
  const Float_t resRPhiIB     = 0.00025;
  const Float_t resZIB        = 0.00025;
  const Float_t resRPhiOB     = 0.00100;
  const Float_t resZOB        = 0.00100;
  const Float_t eff           = 0.98;
  //
  const Float_t bpipe0R       = 0.48;
  const Float_t x0bpipe0      = 0.00042;
  const Float_t xrhobpipe0    = 2.772e-02; // 150 mum Be
  const Float_t bpipe1R       = 3.7 ;
  const Float_t x0bpipe1      = 0.0014;
  const Float_t xrhobpipe1    = 9.24e-02; // 500 mum Be
  const Float_t xOuter        = 100.;
  const Int_t nLayersFull=30;
  //
  Float_t resolNo[2]={0,0};
  Float_t resolIB[2]={resRPhiIB,resRPhiIB};
  Float_t resolOB[2]={resRPhiOB,resRPhiOB};
  //
  // parameters
  float bz=0.5;
  Double_t scaleR=1;
  Int_t    nLayersIB=3;
  Float_t  powerIB=1;
  Int_t    nLayersOB=10;
  Float_t  powerOB=1.0;
  Int_t nLayersAll=2+(nLayersIB+nLayersOB);

  //
  TStopwatch timer;
  timer.Start();
  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticle.root","recreate");
  TTree * tree = 0;
  fastGeometry geom= fastGeometry(nLayersAll);
  geom.fBz=bz;
  fastParticle particle=fastParticle(2*nLayersAll);
  particle.fgStreamer=pcstream;
  //
  geom.setLayer(0,bpipe0R*scaleR, x0bpipe0,xrhobpipe0 , resolNo);
  geom.setLayerRadiusPower(1, nLayersIB+1,bpipe0R*1.02*scaleR, bpipe1R*0.98*scaleR, powerIB,x0IB,xrhoIB, resolIB);
  geom.setLayer(nLayersIB+1, bpipe1R*scaleR, x0bpipe1,xrhobpipe1 , resolNo);
  geom.setLayerRadiusPower(nLayersIB+2, nLayersIB+2+nLayersOB, bpipe1R*1.02*scaleR, xOuter *scaleR, 1, x0OB,xrhoOB , resolOB);
  //

  //geom.setLayerRadiusPower();

  //

  for (Int_t i=0; i<nParticles; i++){
    double r[]     = {0+gRandom->Rndm()*0.000000001,gRandom->Rndm()*0.000000001,0};
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
    int    pidCode=int(gRandom->Rndm()*8);                   // PID code of particles - 0-electron ,1-muon, 2-pion, 3-Kaon
    long  charge  = (gRandom->Rndm()<0.5) ? -1:1;
    long    pdgCode = AliPID::ParticleCode(pidCode)*charge;  // PID code covnerted to the PdgCode
    if (gRandom->Rndm()<kRandomPDGFraction) pdgCode=0;
     Float_t decayLength= (gRandom->Rndm()<kDecayFraction) ?gRandom->Rndm()*geom.fLayerRadius[geom.fLayerRadius.size()-1]:0;
    particle.fDecayLength=decayLength;
    particle.simulateParticle(geom, r,p,pdgCode, kMaxLength,nLayersAll*2);
    particle.reconstructParticle(geom,pdgCode,nLayersAll*2);
    particle.reconstructParticleRotate0(geom,pdgCode,nLayersAll*2);

    if (dumpStream==kFALSE) continue;
    if (tree) tree->Fill();
    else {
      (*pcstream) << "fastPart" <<
                  "i=" << i <<
                  "isSecondary="<<isSecondary<<
                  "pidCode="<<pidCode<<
                  "pdgCode="<<pdgCode<<
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

void initTreeFast(const char * inputList="fastParticle.list"){
  const char* inputListPath=gSystem->ExpandPathName(inputList);
  treeFast  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"fastPart",0,10000);
  treeTurn  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"turn",0,10000);
  treeUnit0  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"UnitTestDumpCorrectForMaterial",0,10000);
  treeSeed  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"seedDump",0,10000);
  treeFast->SetMarkerStyle(21);
  treeFast->SetMarkerSize(0.5);
   treeUnit0->SetMarkerStyle(21);
  treeUnit0->SetMarkerSize(0.5);
  treeSeed->BuildIndex("gid");
  treeFast->BuildIndex("gid");
  treeSeed->AddFriend(treeFast,"F");

  //AliDrawStyle::SetDefaults();
  //AliDrawStyle::ApplyStyle("figTemplate");
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
