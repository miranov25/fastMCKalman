/*
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I$fastMCKalman/fastMCKalman/MC/")
      gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
   .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/examples/fastSimulationRun3.C+
    AliPDG::AddParticlesToPdgDataBase();
    pgeom=createGeometry();
    makeSimulation(*pgeom,10000);
    #
     treeFast  = AliXRDPROOFtoolkit::MakeChainRandom("fastParticle.list","fastPart",0,10000);
     fastParticle::setAliases(*treeFast);
    .> a.log
    makeSimulation(*pgeom,10000);
     .>
 */
#include "fastSimulation.h"
#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
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


/// testAlice 3 configuration as proposed in https://github.com/preghenella/DelphesO2/blob/a058f94f6cb887edcf725fd991d16ca5f7b76e0b/src/lutWrite.werner.cc
fastGeometry * createGeometry() {
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
  fastGeometry *geom= new fastGeometry(nLayerTPC);
  geom->fBz=bz;
  // ITS emulation  0 6 layers wit  O(0.001 cm) resolution
  //     X0 is radiation length per layer and rho is the integrated density per layer - obtained from Geant simulation
  float resol[2]={0.001,0.001};
  geom->setLayerRadiusPower(0,5, 4,37, 1.0,0.006,0.16, resol);
  // TPC emulation
  //    X0 is radiation length per layer and rho is the integrated density per layer - obtained from Geant simulation
  resol[0]=0.1;
  resol[1]=0.1;
  //geom->setLayerRadiusPower(6,nLayerTPC,89,260,1.0,0.000025,0.0009,resol);
  geom->setLayerRadiusPower(6,nLayerTPC,89,260,1.0,xx0,xrho,resol);
  return geom;
}

void makeSimulation(fastGeometry &geom, int nParticles){
  // simulation setup parameters
  const Float_t smearR = 10;
  const Float_t smearZ = 10;
  const Float_t kMinPt = 0.02;
  const Float_t kMax1Pt = 1. / 100.;
  //const Float_t kMax1Pt=1./0.1;
  const Float_t kThetaMax = 2.;
  const Float_t kMaxLength = 300;
  //
  Int_t nLayersAll= geom.fLayerIndex.size();
  TStopwatch timer;
  timer.Start();
  TTreeSRedirector *pcstream = new TTreeSRedirector("fastParticle.root", "recreate");
  fastParticle particle = fastParticle(2 * nLayersAll);
  particle.fgStreamer = pcstream;
  TTree *tree = 0;
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
    long    pdgCode = AliPID::ParticleCode(pidCode)*charge;  // PID code converted to the PdgCode
    if (gRandom->Rndm()<kRandomPDGFraction) pdgCode=0;
     Float_t decayLength= (gRandom->Rndm()<kDecayFraction) ?gRandom->Rndm()*geom.fLayerRadius[geom.fLayerRadius.size()-1]:0;
    particle.fDecayLength=decayLength;
    particle.simulateParticle(geom, r,p,pdgCode, kMaxLength,nLayersAll*2);
    particle.reconstructParticle(geom,pdgCode,nLayersAll*2);
    particle.reconstructParticleRotate0(geom,pdgCode,nLayersAll*2);

    //if (dumpStream==kFALSE) continue;
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

