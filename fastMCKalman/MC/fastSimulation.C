/*
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.C

    geom= fastGeometry(201.)
    particle=fastParticle(200);
    geom.setLayerRadiusPower(0,10,1,20,1.2,0);
    geom.setLayerRadiusPower(11,200,20,150,1.2,0);

*/
#include "TObject.h"
#include "AliExternalTrackParam.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "ROOT/RVec.hxx"
//using namespace ROOT;
using namespace ROOT::VecOps;
const Int_t kMaxLayers=10000;
#pragma link C++ class RVec<AliExternalTrackParam>+;

class fastGeometry:public TObject{
public:
  fastGeometry():TObject(){}
  ~fastGeometry(){}
  fastGeometry(Int_t maxLayer){
    fLayerRadius.resize(maxLayer); fLayerX0.resize(maxLayer);fLayerRho.resize(maxLayer);fLayerIndex.resize(maxLayer);
  }
  void setLayerRadiusPower(int layer0, int layerN, float r0, float rN, float power, float X0,float rho);
  void Print();
  RVec<float> fLayerRadius;  //   barrel radius array
  RVec<int>  fLayerIndex;   //   indices for fast navigation
  RVec<float> fLayerX0;      //   barrel X0 array
  RVec<float> fLayerRho;      //   barrel Rho array
  float              fBz;           //   magnetic field
  ClassDef(fastGeometry, 1)
};

class fastParticle : public TObject{
public:
    fastParticle():TObject(){}
    ~fastParticle(){}
    fastParticle(int nLayers){
        fPosition.reserve(nLayers); fLayerIndex.reserve(nLayers); fParamIn.reserve(nLayers); fParamOut.reserve(nLayers); fParamMC.reserve(nLayers);
    }
    int simulateParticle(fastGeometry  &geom, double r[3], double p[3], int pdgCode, float maxLength, int maxPoints);
    double fR[3];    // initial position
    double fP[3];    // initial momentum
    Int_t   fPdgCode;        // PDG code
    RVec<float>    fPosition;                //   layer position
    RVec<int>      fLayerIndex;              //   layer index
    RVec<AliExternalTrackParam> fParamMC;    //   "simulate"      Param MC
    RVec<AliExternalTrackParam> fParamOut;   //   "reconstructed" Param Out
    RVec<AliExternalTrackParam> fParamIn;    //   "reconstructed" Param In
    ClassDef(fastParticle, 1)
};
ClassImp(fastGeometry)
ClassImp(fastParticle)


/// make layer distribution in interval layer0(r0)->layerN(rN) following power low distance
/// \param layer0
/// \param layerN
/// \param r0
/// \param rN
/// \param power
/// \param X0
void fastGeometry::setLayerRadiusPower(int layer0, int layerN, float r0, float rN, float power, float X0,float rho){
    float dLayerN=(layerN-layer0);
    for (Int_t iLayer=layer0; iLayer<=layerN; iLayer++){
        fLayerX0[iLayer]=X0;
        fLayerRho[iLayer]=rho;
        float radius=r0+TMath::Power(float(iLayer-layer0)/dLayerN,power)*(rN-r0);
        fLayerRadius[iLayer]=radius;
    }
    fLayerIndex = Argsort(fLayerRadius);
}

int fastParticle::simulateParticle(fastGeometry  &geom, double r[3], double p[3], int pdgCode, float maxLength, int maxPoints){
  double covar[14]={};
  fPdgCode=pdgCode;
  TParticlePDG * particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
  if (particle== nullptr) {
    ::Error("fastParticle::simulateParticle","Invalid pdgCode %d",pdgCode);
    return -1;
  }
  float sign = particle->Charge()/3.;
  float mass = particle->Mass();
  AliExternalTrackParam param(r,p,covar,sign);
  float length=0, time=0;
  float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  uint indexR= uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
  int nPoint=0;
  fParamMC.resize(1);
  fParamMC[0]=param;
  for (nPoint=1; nPoint<maxPoints&&length<maxLength; nPoint++){
    //printf("%d\n",nPoint);
    //param.Print();
    if (indexR>geom.fLayerRadius.size()) break;
    double xyz[3];
    float radius  =geom.fLayerRadius[indexR];
    float xrho  =geom.fLayerRho[indexR];
    float xx0  =geom.fLayerX0[indexR];
    //Float_t x = param.GetXatLabR(r,localX,fBz,1);
    int status =  param.GetXYZatR(radius,geom.fBz,xyz);
    if (status==0) break;
    double alpha=TMath::ATan2(xyz[1],xyz[0]);
    param.Rotate(alpha);
    status = param.PropagateTo(radius,geom.fBz);
    param.CorrectForMeanMaterial(xx0,xrho,mass);
    fParamMC.resize(nPoint+1);
    fParamMC[nPoint]=param;
    indexR++;
  }

}