/*
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+
    geom= fastGeometry(201.)
    particle=fastParticle(200);
    geom.setLayerRadiusPower(0,10,1,20,1.2,0);
    geom.setLayerRadiusPower(11,200,20,150,1.2,0);
*/

#include "fastSimulation.h"
//#include "AliHelix.h"
#include "TTreeStream.h"

TTreeSRedirector* fastParticle::fgStreamer = nullptr;

ClassImp(fastGeometry)
ClassImp(fastParticle)

void getHelix(Double_t *fHelix,  AliExternalTrackParam t, float bz){
  // addapted code from AliHelix
  Double_t x,cs,sn;
  x=t.GetX();
  for (Int_t i=0;i<5; i++) fHelix[i]=t.GetParameter()[i];
  float alpha = t.GetAlpha();
  //
  //circle parameters
  //PH Sometimes fP4 and fHelix[4] are very big and the calculation
  //PH of the Sqrt cannot be done. To be investigated...
  fHelix[4]=fHelix[4]/(-1000/0.299792458/bz);    // C
  //  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  cs=cosf(alpha); sn=sinf(alpha); // RS use float versions: factor 2 in CPU speed
  Double_t xc, yc, rc;
  rc  =  1/fHelix[4];
  xc  =  x-fHelix[2]*rc;
  Double_t dummy = 1-(x-xc)*(x-xc)*fHelix[4]*fHelix[4];
  if (dummy<0) dummy = 0;
  yc  =  fHelix[0]+TMath::Sqrt(dummy)/fHelix[4];
  fHelix[6] = xc*cs - yc*sn;
  fHelix[7] = xc*sn + yc*cs;
  fHelix[8] =  TMath::Abs(rc);
  return;
}


/// make layer distribution in interval layer0(r0)->layerN(rN) following power low distance
/// \param layer0
/// \param layerN
/// \param r0
/// \param rN
/// \param power
/// \param X0
void fastGeometry::setLayerRadiusPower(int layer0, int layerN, float r0, float rN, float power, float X0,float rho, float resol[2]){
    float dLayerN=(layerN-layer0);
    for (Int_t iLayer=layer0; iLayer<=layerN; iLayer++){
        fLayerX0[iLayer]=X0;
        fLayerRho[iLayer]=rho;
        float radius=r0+TMath::Power(float(iLayer-layer0)/dLayerN,power)*(rN-r0);
        fLayerRadius[iLayer]=radius;
        fLayerResolRPhi[iLayer]=resol[0];
        fLayerResolZ[iLayer]=resol[1];
    }
    fLayerIndex = Argsort(fLayerRadius);
}

/// simulate particle    - barrel part, end cup and loopers not yet implemented
/// \param geom          - geometry and B field description
/// \param r             - initial position
/// \param p             - initial momenta
/// \param pdgCode       - PDG code of particle (pion, Kaon,proton,electron,muon)
/// \param maxLength     - max length to simulate
/// \param maxPoints     - maximal number of points to simulate
/// \return              - modify status of particles = create points along   - TODO status flags to be decides
int fastParticle::simulateParticle(fastGeometry  &geom, double r[3], double p[3], int pdgCode, float maxLength, int maxPoints){
    const float kMaxSnp=0.90;
   double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  fPdgCodeMC=pdgCode;
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
  float direction=r[0]*p[0]+r[1]*p[1]+r[2]*p[2];
  direction=(direction>0)? 1.:-1.;
  if (radius==0) direction=1;
  uint indexR= uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
  int nPoint=0;
  fParamMC.resize(1);
  fParamMC[0]=param;
  fLayerIndex[0]=indexR;
  double *par = (double*)param.GetParameter();
  for (nPoint=1; nPoint<maxPoints&&length<maxLength; nPoint++){
    //printf("%d\n",nPoint);
    //param.Print();
    if (indexR>geom.fLayerRadius.size()) {
      break;
    }
    double xyz[3],pxyz[3];
    float radius  = geom.fLayerRadius[indexR];
    //float xrho    = geom.fLayerRho[indexR];
    //float xx0     = geom.fLayerX0[indexR];
    //Float_t x = param.GetXatLabR(r,localX,fBz,1);
    int status =  param.GetXYZatR(radius,geom.fBz,xyz);
    if (status==0){   // if not possible to propagate to next radius - assume looper - change direction
      param.GetPxPyPz(pxyz);
      float C         = param.GetC(geom.fBz);
      float R         = TMath::Abs(1/C);
      float dca       = param.GetD(0,0,geom.fBz);         // signed distance of closest approach to origin (0,0) -- looks like  buggy in some cases - dca can be negative with dca<2*R
      double helix[9];
      getHelix(helix,param,geom.fBz);
      double x0=helix[6],y0=helix[7],r0=TMath::Sqrt(x0*x0+y0*y0);
      float dcaMI     = r0-R;
      float dlaMI     = r0+R;
      float dla       = dca+2*R;                          // distance of longest approach
      int  isUp       = abs(radius-dcaMI)>abs(radius-dlaMI)? 1:0;
      int        aSign= (param.GetSnp()>0)? 1:-1.;
      float alpha     = param.GetAlpha();
      float alphaC    = TMath::ATan2(y0,x0);
      float dAlpha    = alpha-alphaC;
      if (dAlpha>TMath::Pi()) dAlpha-=TMath::TwoPi();
      if (dAlpha<-TMath::Pi()) dAlpha+=TMath::TwoPi();
      float alpha2 = alphaC-dAlpha;
      double param2[5];
      param2[0]=0;
      param2[1]=param.GetParameter()[1];
      param2[2]=-param.GetParameter()[2];
      param2[3]=-param.GetParameter()[3];
      param2[4]=-param.GetParameter()[4];
      AliExternalTrackParam paramNew(param.GetX(),alpha2,param2,param.GetCovariance());
      if (fgStreamer){
        (*fgStreamer)<<"turn"<<
          "radius="<<radius<<     // radius to propagate
          "direction="<<direction<<
          "x0="<<x0<<
          "y0="<<y0<<
          "dca="<<dca<<
          "dcaMI="<<dcaMI<<
          "dlaMI="<<dlaMI<<
          "param.="<<&param<<
          "paramNew.="<<&paramNew<<
          "\n";
      }
      param=paramNew;
      direction*=-1;
    }else{
      double alpha  = TMath::ATan2(xyz[1],xyz[0]);
      if (fStatus.size()<=nPoint) fStatus.resize(nPoint+1);
      fStatus[nPoint].resize(6);                  // resize status array
      status = param.Rotate(alpha);
      fStatus[nPoint][0]=status;
      status = param.PropagateTo(radius,geom.fBz);
      fStatus[nPoint][1]=status;
    }
    //
    float xrho    = geom.fLayerRho[indexR];
    float xx0     = geom.fLayerX0[indexR];
    float tanPhi2 = par[2]*par[2];
    tanPhi2=(1-tanPhi2);
    float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);               /// geometrical path assuming crossing cylinder
    param.CorrectForMeanMaterial(crossLength*xx0,-crossLength*xrho,mass);
    fParamMC.resize(nPoint+1);
    fParamMC[nPoint]=param;
    fLayerIndex[nPoint]=indexR;
    indexR+=direction;
  }
  return 1;
}

/*

int fastParticle::simulateParticle(fastGeometry  &geom, double r[3], double p[3], int pdgCode, float maxLength, int maxPoints){
  const float kMaxSnp=0.90;
  double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  fPdgCodeMC=pdgCode;
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
  float direction=r[0]*p[0]+r[1]*p[1]+r[2]*p[2];
  direction=(direction>0)? 1.:-1.;
  if (radius==0) direction=1;
  uint indexR= uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
  int nPoint=0;
  fParamMC.resize(1);
  fParamMC[0]=param;
  fLayerIndex[0]=indexR;
  double *par = (double*)param.GetParameter();
  for (nPoint=1; nPoint<maxPoints&&length<maxLength; nPoint++){
    //printf("%d\n",nPoint);
    //param.Print();
    if (indexR>geom.fLayerRadius.size()) {
      break;
    }
    double xyz[3],pxyz[3];
    float radius  = geom.fLayerRadius[indexR];
    //float xrho    = geom.fLayerRho[indexR];
    //float xx0     = geom.fLayerX0[indexR];
    //Float_t x = param.GetXatLabR(r,localX,fBz,1);
    int status =  param.GetXYZatR(radius,geom.fBz,xyz);
    if (status==0){   // if not possible to propagate to next radius - assume looper - change direction
      float C         = param.GetC(geom.fBz);
      float R         = TMath::Abs(1/C);
      float dca       = param.GetD(0,0,geom.fBz);         // distance of closest approach to origin (0,0)
      float dla       = dca+2*R;                          // distance of longest approach
      bool  isUp      = abs(Radius-dca)>abs(Radius-dla);
      float dAlpha    = abs(TMath::ASin(param->GetSnp()))-TMath::Pi()*0.5;
      dAlpha *=(param->GetSnp()<0):-1.:1.;
      float alphaNew  = param.GetAlpha()+dAlpha;
      double paramNew[5]={par[0],par[1],par[2],par[3],par[4]};
      AliExternalTrackParam paramNew=param;
      double xyzNew[3],pxyzNew[3];
      paramNew.GetPxPyPz(pxyzNew);
      paramNew.GetXYZ(xyzNew);
      float stepDir=(direction*C*sign)>0? -1:1;
      double alphaDir  = TMath::ATan2(pxyzNew[1],pxyzNew[0]);
      if (stepDir<0) alphaDir+=TMath::Pi();
      double alpha     = TMath::ATan2(xyzNew[1],xyzNew[0]);
      float dPhi = TMath::Abs(TMath::ASin(paramNew.GetSnp()));
      float step = TMath::Cos(dPhi)*TMath::Abs(R);
      status = paramNew.Rotate(alphaDir);
      if (status==0) {
        param.Print();
        paramNew.Print();
        printf("1.) Can not rotate: \t%f\t%f\t%f\n\n",direction,C,sign);
        break;
      }

      float newX=paramNew.GetX()+stepDir*TMath::Abs(step*2);
      status = paramNew.PropagateTo(newX,geom.fBz);
      if (status==0) {
        printf("2.)\n");
        param.Print();
        paramNew.Print();
        AliExternalTrackParam paramNew2=param;
        status = paramNew2.Rotate(-alphaDir);
        float newX=paramNew2.GetX()+TMath::Abs(step*2);
        bool status2 = paramNew2.PropagateTo(newX,geom.fBz);
        paramNew2.Print();
        printf("2.) Can not propagate:\t%f\t%f\t%f\t%d\t%d\n\n\n",direction,C,sign, (direction*C*sign)>0,status2, isUp);
        break;
      }
      paramNew.GetXYZ(xyzNew);
      double alphaNew     = TMath::ATan2(xyzNew[1],xyzNew[0])+TMath::Pi();
      status = paramNew.Rotate(alphaNew);
      if (status==0) {
        printf("\n");
        param.Print();
        paramNew.Print();
        printf("3.) Can not rotate new: \t%f\t%f\t%f\t%d\n\n\n",direction,C,sign, (direction*C*sign)>0);
        break;
      }
      paramNew.SetParamOnly(-paramNew.GetX(),alphaNew-TMath::Pi(),paramNew.GetParameter());
      double *parNew = (double*)paramNew.GetParameter();
      parNew[4]*=-1; parNew[3]*=-1;
      printf("param print begin\n");
      param.Print();
      paramNew.Print();
      status = paramNew.PropagateTo(param.GetX(),geom.fBz);
      if (status==0) {
        printf("\n");
        param.Print();
        paramNew.Print();
        printf("4.) Can not propagate: \t%f\t%f\t%f\n\n",direction,C,sign);
      }
      paramNew.Print();
      param=paramNew;
      direction*=-1;
      indexR=fLayerIndex[nPoint-1];
      par = (double*)param.GetParameter();
    }else{
      double alpha  = TMath::ATan2(xyz[1],xyz[0]);
      status = param.Rotate(alpha);
      status = param.PropagateTo(radius,geom.fBz);
    }
    //
    float xrho    = geom.fLayerRho[indexR];
    float xx0     = geom.fLayerX0[indexR];
    float tanPhi2 = par[2]*par[2];
    tanPhi2=(1-tanPhi2);
    float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);               /// geometrical path assuming crossing cylinder
    param.CorrectForMeanMaterial(-crossLength*xx0,-crossLength*xrho,mass);
    fParamMC.resize(nPoint+1);
    fParamMC[nPoint]=param;
    fLayerIndex[nPoint]=indexR;
    indexR+=direction;
  }
  return 1;
}

*/

///
/// \param geom          - pointer to geomtery to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticle(fastGeometry  &geom, int pdgCode, uint layerStart){
  const Float_t chi2Cut=16;
  const float kMaxSnp=0.90;
  double covar0[5]={0.1,0.1,0.001,0.001,0.01};
  fPdgCodeRec   =pdgCode;
  TParticlePDG * particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
  if (particle== nullptr) {
    ::Error("fastParticle::simulateParticle","Invalid pdgCode %d",pdgCode);
    return -1;
  }
  float sign = particle->Charge()/3.;
  float mass = particle->Mass();
  uint layer1 = TMath::Min(layerStart,uint(fParamMC.size()-1));
  AliExternalTrackParam param=fParamMC[layer1];
  double *covar = (double*)param.GetCovariance();
  covar[0]=covar0[0];
  covar[2]=covar0[1];
  covar[5]=covar0[2];
  covar[9]=covar0[3];
  covar[14]=covar0[4];
  //
  float length=0, time=0;
  float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fParamIn.resize(layer1+1);
  fParamIn[layer1]=param;
  double xyz[3];
  int status=0;
  const double *par = param.GetParameter();
  for (int layer=layer1-1; layer>=1; layer--){   // dont propagate to vertex , will be done later ...
      double resol=0;
      AliExternalTrackParam & p = fParamMC[layer];
      p.GetXYZ(xyz);
      double alpha=TMath::ATan2(xyz[1],xyz[0]);
      double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
      status = param.Rotate(alpha);
      if (status==0){
        ::Error("reconstructParticle", "Rotation failed");
        break;
      }
      if (fStatus[layer].size()==0) fStatus[layer].resize(6);
      fStatus[layer][2]=status;
      status = param.PropagateTo(radius,geom.fBz);
      if (status==0){
        ::Error("reconstructParticle", "Propagation  failed");
        break;
      }
      fStatus[layer][3]=status;
      float xrho  =geom.fLayerRho[layer];
      float xx0  =geom.fLayerX0[layer];
      double pos[2]={0,xyz[2]};
      double cov[3]={geom.fLayerResolRPhi[layer]*geom.fLayerResolRPhi[layer],0, geom.fLayerResolZ[layer]*geom.fLayerResolZ[layer]};
      fParamIn[layer]=param;
      float chi2 =  param.GetPredictedChi2(pos, cov);
      if (chi2>chi2Cut){
        ::Error("reconstructParticle", "Too big chi2 %f", chi2);
        break;
      }

      if (TMath::Abs(param.GetSnp())<kMaxSnp) {
        param.Update(pos, cov);
      }
      float tanPhi2 = par[2]*par[2];
      tanPhi2=(1-tanPhi2);
      float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);                /// geometrical path assuming crossing cylinder
      param.CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
  }
  return 1;
}