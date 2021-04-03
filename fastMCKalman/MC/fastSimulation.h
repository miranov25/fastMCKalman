/*
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.h
*/
#ifndef fastSimulation_H
#define fastSimulation_H

#include "TObject.h"
#include "AliExternalTrackParam.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "ROOT/RVec.hxx"
//using namespace ROOT;
using namespace ROOT::VecOps;
const Int_t kMaxLayers=10000;
#pragma link C++ class RVec<AliExternalTrackParam>+;


class TTreeSRedirector;

class fastGeometry:public TObject{
public:
  fastGeometry():TObject(){}
  ~fastGeometry(){}
  fastGeometry(Int_t maxLayer){
    fLayerRadius.resize(maxLayer); fLayerX0.resize(maxLayer);fLayerRho.resize(maxLayer);fLayerIndex.resize(maxLayer);fLayerResolRPhi.resize(maxLayer);fLayerResolZ.resize(maxLayer);
  }
  void setLayerRadiusPower(int layer0, int layerN, float r0, float rN, float power, float X0,float rho, float resol[2]);
  void Print();
  RVec<float> fLayerRadius;        //   barrel radius array
  RVec<int>   fLayerIndex;         //   barrel indices for fast navigation
  RVec<float> fLayerX0;            //   barrel X0 array
  RVec<float> fLayerRho;           //   barrel Rho array
  RVec<float> fLayerResolRPhi;     //   rphi resolution at layer
  RVec<float> fLayerResolZ;        //   z  resolution at layer
  RVec<float> fHitDensity;         //   hit density (hits per cm^2)
  float       fBz;                 //   magnetic field
  ClassDef(fastGeometry, 1)
};

class fastParticle : public TObject{
public:
  fastParticle():TObject(){}
  ~fastParticle(){}
  fastParticle(int nLayers){
    fLayerIndex.reserve(nLayers); fDirection.reserve(nLayers); fParamIn.reserve(nLayers); fParamInRot.reserve(nLayers);
    fParamOut.reserve(nLayers); fParamMC.reserve(nLayers);fStatus.reserve(nLayers);
    fChi2.resize(nLayers);
    fMaxLayer=0;
  }
  int simulateParticle(fastGeometry     &geom, double r[3], double p[3], int pdgCode, float maxLength, int maxPoints);
  int reconstructParticle(fastGeometry  &geom, int pdgCode, uint layerStart);
  int reconstructParticleRotate0(fastGeometry  &geom, int pdgCode, uint layerStart);
  double fR[3];                            //   initial position
  double fP[3];                            //   initial momentum
  int                         fPdgCodeMC;  //   PDG code used in simulation
  int                         fPdgCodeRec; //   PDG code as used in reconstruction
  int                         fMaxLayer;   //   maximal layer position
  RVec<int>                   fLayerIndex; //   layer index    - important for loopers
  RVec<float>                 fDirection;  //   particle direction - Out=1, In = -1
  RVec<AliExternalTrackParam> fParamMC;    //   "simulate"      Param MC
  RVec<AliExternalTrackParam> fParamOut;   //   "reconstructed" Param Out
  RVec<AliExternalTrackParam> fParamIn;    //   "reconstructed" Param In
  RVec<AliExternalTrackParam> fParamInRot;    //   "reconstructed" Param In - in rotated frame
  RVec<std::vector<int>>      fStatus;     //   propagation/update status
  RVec<float>                 fChi2;      //   chi2  at layer
  //                                       - local information - to be assigned to simulated track - if not set  taken symmetric from fastGeometry
  RVec<float> fLayerResolRPhi;             //   rphi resolution at layer
  RVec<float> fLayerResolZ;                //   z  resolution at layer
  RVec<float> fLayerDeltaRPhi;             //   rphi delta distortion at layer
  RVec<float> fLayerDeltaZ;                //   z  delta distortion at layer
  RVec<float> fLayerProb;                  //   probability  to assign point
  RVec<float> fHitDensity;                 //   hit density
  RVec<bool> fLayerInDead;                 //   flag in that region
  static TTreeSRedirector * fgStreamer;    //   debug streamer
  ClassDef(fastParticle, 1)
};
#endif

