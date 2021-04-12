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
ClassImp(AliExternalTrackParam4D)


AliExternalTrackParam4D::AliExternalTrackParam4D():AliExternalTrackParam(){}
AliExternalTrackParam4D::AliExternalTrackParam4D(const AliExternalTrackParam &t):AliExternalTrackParam(t){}

/// Runge-Kuta energy loss correction  - https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
/// \param xOverX0        - X/X0, the thickness in units of the radiation length.
/// \param xTimesRho      - is the product length*density (g/cm^2).
//                        - It should be passed as negative when propagating tracks
//                        - from the interaction point to the outside of the central barrel.
/// \param mass           - the mass of this particle (GeV/c^2). Negative mass means charge=2 particle
/// \param f              - dEdx formula
/// \param stepFraction   - step fraction
/// \return
Bool_t AliExternalTrackParam4D::CorrectForMeanMaterialRK(Double_t xOverX0, Double_t xTimesRho,Double_t mass, Double_t (*f)(Double_t), Float_t stepFraction){
  //  Runge Kuttta integral p(x)
  //k_{1} is the slope at the beginning of the interval, using {\displaystyle y}y (Euler's method);
  //k_{2} is the slope at the midpoint of the interval, using {\displaystyle y}y and {\displaystyle k_{1}}k_{1};
  //k_{3} is again the slope at the midpoint, but now using {\displaystyle y}y and {\displaystyle k_{2}}k_{2};
  //k_{4} is the slope at the end of the interval, using {\displaystyle y}y and {\displaystyle k_{3}}k_{3}.
  const Double_t kBGStop=0.02;
  Double_t p=GetP();
  Double_t q=(mass<0)?1.:2.;   // q=2 particle in ALICE convention
  mass=TMath::Abs(mass);
  Double_t mass2=mass*mass;
  p*=q;
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t p2=p*p;
  Double_t dEdxM=f(p/mass),dEdxMRK=0;
  Double_t Ein=TMath::Sqrt(p2+mass2);
  Double_t Eout=0;
  //Runge-Kuta
  if ( TMath::Abs(dEdxM*xTimesRho) > 0.3*Ein ) return kFALSE; //30% energy loss is too much!
  if (TMath::Abs(dEdxM*xTimesRho)>stepFraction*Ein || (Ein+dEdxM*xTimesRho)<mass){
    Double_t E1=Ein,E2=0,E3=0,E4=0;
    Double_t k1=0,k2=0,k3=0,k4=0;
    k1=dEdxM;
    E2=E1+k1*xTimesRho*0.5;
    if (E2*E2/mass2-1.< kBGStop*kBGStop) return kFALSE;
    k2=f(TMath::Sqrt(E2*E2/mass2-1.));
    E3=E1+k2*xTimesRho*0.5;
    if (E3*E3/mass2-1.< kBGStop*kBGStop) return kFALSE;
    k3=f(TMath::Sqrt(E3*E3/mass2-1.));
    E4=E1+k3*xTimesRho;
    if (E4*E4/mass2-1.< kBGStop*kBGStop) return kFALSE;
    k4 =f(TMath::Sqrt(E4*E4/mass2-1.));
    dEdxMRK=(k1+2.*k2+2.*k3+k4)/6.;
    Eout=E1+xTimesRho*dEdxMRK;
  }else{
    Eout=Ein+xTimesRho*dEdxM;
  }

  Double_t Emean=(Ein+Eout)*0.5;
  if (Emean<mass) {
    return kFALSE;
  }
  // Use mean values for p2, p and beta2
  p2=Emean*Emean-mass2;
  p =TMath::Sqrt(p2);
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t beta2=p2/(p2+mass2);
  //
  Double_t &fP2=fP[2];
  Double_t &fP3=fP[3];
  Double_t &fP4=fP[4];
  Double_t &fC22=fC[5];
  Double_t &fC33=fC[9];
  Double_t &fC43=fC[13];
  Double_t &fC44=fC[14];
  //
  //Calculating the multiple scattering corrections******************
  Double_t cC22 = 0.;
  Double_t cC33 = 0.;
  Double_t cC43 = 0.;
  Double_t cC44 = 0.;
  if (xOverX0 != 0) {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    if (GetUseLogTermMS()) {
      double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
      if (lt>0) theta2 *= lt*lt;
    }
    theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) return kFALSE;
    cC22 = theta2*((1.-fP2)*(1.+fP2))*(1. + fP3*fP3);
    cC33 = theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
    cC43 = theta2*fP3*fP4*(1. + fP3*fP3);
    cC44 = theta2*fP3*fP4*fP3*fP4;
  }

  //Calculating the energy loss corrections************************
  Double_t cP4=1.;
  if ((xTimesRho != 0.) && (beta2 < 1.)) {
    Double_t dE=Eout-Ein;
    if ( (1.+ dE/p2*(dE + 2*Ein)) < 0. ) return kFALSE;
    cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*Ein));  //A precise formula by Ruben !
    //if (TMath::Abs(fP4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c -dsiable controlled by the BG cut
    // Approximate energy loss fluctuation (M.Ivanov)
    const Double_t knst=0.07; // To be tuned.
    Double_t sigmadE=knst*TMath::Sqrt(TMath::Abs(dE));
    cC44 += ((sigmadE*Ein/p2*fP4)*(sigmadE*Ein/p2*fP4));

  }

  //Applying the corrections*****************************
  fC22 += cC22;
  fC33 += cC33;
  fC43 += cC43;
  fC44 += cC44;
  fP4  *= cP4;
  CheckCovariance();
  return kTRUE;
}



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
  fMaxLayer=0;
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
  AliExternalTrackParam4D param(AliExternalTrackParam(r,p,covar,sign));
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
    if (fStatusMaskMC.size()<=nPoint) fStatusMaskMC.resize(nPoint+1);
    //
    double xyz[3],pxyz[3];
    float radius  = geom.fLayerRadius[indexR];
    //float xrho    = geom.fLayerRho[indexR];
    //float xx0     = geom.fLayerX0[indexR];
    //Float_t x = param.GetXatLabR(r,localX,fBz,1);
    int status =  param.GetXYZatR(radius,geom.fBz,xyz);
    if (status==0){   // if not possible to propagate to next radius - assume looper - change direction
      break;    //this is temporary
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
      fStatusMaskMC[nPoint]=0;
      status = param.Rotate(alpha);
      if (status) {
        fStatusMaskMC[nPoint]|= kTrackRotate;
      }else{
        break;
      }
      status = param.PropagateTo(radius,geom.fBz);
      if (status) {
        fStatusMaskMC[nPoint]|=kTrackPropagate;
      }else{
        break;
      }
    }
    //
    float xrho    = geom.fLayerRho[indexR];
    float xx0     = geom.fLayerX0[indexR];
    float tanPhi2 = par[2]*par[2];
    tanPhi2=(1-tanPhi2);
    float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);               /// geometrical path assuming crossing cylinder
    status = param.CorrectForMeanMaterialRK(crossLength*xx0,-crossLength*xrho,mass);
    if (status) {
        fStatusMaskMC[nPoint]|=kTrackCorrectForMaterial;
      }else{
        break;
    }
    fParamMC.resize(nPoint+1);
    fParamMC[nPoint]=param;
    fLayerIndex[nPoint]=indexR;
    indexR+=direction;
    if (indexR>fMaxLayer) fMaxLayer=indexR;
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
/// \param geom          - pointer to geometry to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticle(fastGeometry  &geom, int pdgCode, uint layerStart){
  const Float_t chi2Cut=16;
  const float kMaxSnp=0.90;
  fLengthIn=0;
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
  double covar0[5]={1,1,0.001*(1+1./param.Pt()),0.001*(1+1./param.Pt()),0.01*(1.+1./param.Pt())};
  covar[0]=covar0[0];
  covar[2]=covar0[1];
  covar[5]=covar0[2];
  covar[9]=covar0[3];
  covar[14]=covar0[4];
  //
  float length=0, time=0;
  float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fParamIn.resize(layer1+1);
  fStatusMaskIn.resize(layer1+1);
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
      fStatusMaskIn[layer]=0;
      status = param.Rotate(alpha);
      if (status) {
        fStatusMaskIn[layer]|=kTrackRotate;
      }else{
        ::Error("reconstructParticle", "Rotation failed");
        break;
      }
      status = param.PropagateTo(radius,geom.fBz);
      if (status) {
        fStatusMaskIn[layer]|=kTrackPropagate;
      }else{
        ::Error("reconstructParticle", "Proapagation failed");
        break;
      }
      float xrho  =geom.fLayerRho[layer];
      float xx0  =geom.fLayerX0[layer];
      double pos[2]={0,xyz[2]};
      double cov[3]={geom.fLayerResolRPhi[layer]*geom.fLayerResolRPhi[layer],0, geom.fLayerResolZ[layer]*geom.fLayerResolZ[layer]};
      fParamIn[layer]=param;
      float chi2 =  param.GetPredictedChi2(pos, cov);
      fChi2[layer]=chi2;
      if (chi2<chi2Cut) {
        fStatusMaskIn[layer]|=kTrackChi2;
      }else{
        ::Error("reconstructParticle", "Too big chi2 %f", chi2);
        break;
      }
      if (TMath::Abs(param.GetSnp())<kMaxSnp) {
        status = param.Update(pos, cov);
        if (status) {
          fStatusMaskIn[layer]|=kTrackUpdate;
        }else{
          ::Error("reconstructParticle", "Update failed");
          break;
        }
      }

      float tanPhi2 = par[2]*par[2];
      tanPhi2=(1-tanPhi2);
      float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);                /// geometrical path assuming crossing cylinder
      status = param.CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
      if (status) {
        fStatusMaskIn[layer]|=kTrackCorrectForMaterial;
      }else{
        ::Error("reconstructParticle", "Correct for material failed");
        break;
      }
      fLengthIn++;
  }
  return 1;
}



///
/// \param geom          - pointer to geometry to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticleRotate0(fastGeometry  &geom, int pdgCode, uint layerStart){
  const Float_t chi2Cut=16;
  const float kMaxSnp=0.90;
  double covar0[5]={0.1,0.1,0.001,0.001,0.01};
  fLengthInRot=0;
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
  //
  Double_t sinT=TMath::Sin(param.GetAlpha());
  Double_t cosT=TMath::Cos(param.GetAlpha());
  //
  double *covar = (double*)param.GetCovariance();
  covar[0]=covar0[0];
  covar[2]=covar0[1];
  covar[5]=covar0[2];
  covar[9]=covar0[3];
  covar[14]=covar0[4];
  //
  float length=0, time=0;
  float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fParamInRot.resize(layer1+1);
    fStatusMaskInRot.resize(layer1+1);
  fParamInRot[layer1]=param;
  double xyz[3];
  int status=0;
  const double *par = param.GetParameter();
  for (int layer=layer1-1; layer>=1; layer--){   // dont propagate to vertex , will be done later ...
    double resol=0;
    AliExternalTrackParam & p = fParamMC[layer];
    p.GetXYZ(xyz);
    double alpha=TMath::ATan2(xyz[1],xyz[0]);
    double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    Double_t localX, localY;
    localX=cosT*xyz[0]+sinT*xyz[1];
    localY=-sinT*xyz[0]+cosT*xyz[1];
    fStatusMaskInRot[layer]=0;
    status = param.PropagateTo(localX,geom.fBz);
    if (status) {
      fStatusMaskInRot[layer]|=kTrackPropagate;
    }else{
      ::Error("reconstructParticle", "Propagation failed");
      break;
    }
    //

    float xrho  =geom.fLayerRho[layer];
    float xx0  =geom.fLayerX0[layer];
    double pos[2]={localY,xyz[2]};
    double cov[3]={geom.fLayerResolRPhi[layer]*geom.fLayerResolRPhi[layer],0, geom.fLayerResolZ[layer]*geom.fLayerResolZ[layer]};
    fParamInRot[layer]=param;    ///
    float chi2 =  param.GetPredictedChi2(pos, cov);
    fChi2[layer]=chi2;
    if (chi2<chi2Cut) {
      fStatusMaskInRot[layer] |= kTrackChi2;
    }else {
      ::Error("reconstructParticle", "Chi2 -  failed");
      break;
    }
    //
    if (TMath::Abs(param.GetSnp())<kMaxSnp) {
      status = param.Update(pos, cov);
      if (status) {
        fStatusMaskInRot[layer] |= kTrackUpdate;
      } else {
        ::Error("reconstructParticle", "Update failed");
        break;
      }
    }
    float tanPhi2 = par[2]*par[2];
    tanPhi2=(1-tanPhi2);
    float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);                /// geometrical path assuming crossing cylinder  = to be racalculated
    status = param.CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
    if (status) {
      fStatusMaskInRot[layer] |= kTrackCorrectForMaterial;
    } else {
      ::Error("reconstructParticle", "Correct for material failed");
      break;
    }
      fLengthInRot=0;
  }
  return 1;
}


/// set derived varaibles as alaiases to tree
/// \param tree
void fastParticle::setAliases(TTree & tree){
  tree.SetAlias("gxIn","cos(part.fParamIn[].fAlpha)*part.fParamIn[].fX");
  tree.SetAlias("gyIn","sin(part.fParamIn[].fAlpha)*part.fParamIn[].fX");
  tree.SetAlias("gxMC","cos(part.fParamMC[].fAlpha)*part.fParamMC[].fX");
  tree.SetAlias("gyMC","sin(part.fParamMC[].fAlpha)*part.fParamMC[].fX");
  tree.SetAlias("gzMC","part.fParamMC[].fP[1]");
  tree.SetAlias("rMC","part.fParamMC[].fX");
  //
  tree.SetAlias("ptMC","part.fParamMC[0].fData.Pt()");
  tree.SetAlias("ptMCL","part.fParamMC[].fData.Pt()");
  tree.SetAlias("pMC","part.fParamMC[0].fData.P()");
  tree.SetAlias("ptRec","part.fParamIn[0].fData.Pt()");
  tree.SetAlias("tglMC","part.fParamIn[0].fData.fP[3]");
  //
  //
  tree.SetAlias("layer","Iteration$");
  tree.SetAlias("c0MC","sqrt(part.fParamMC[].fC[0])");
  tree.SetAlias("c2MC","sqrt(part.fParamMC[].fC[2])");
  tree.SetAlias("c14MC","sqrt(part.fParamMC[].fC[14])");
  tree.SetAlias("c0In","sqrt(part.fParamIn[].fC[0])");
  tree.SetAlias("c2In","sqrt(part.fParamIn[].fC[2])");
  tree.SetAlias("c0InRot","sqrt(part.fParamInRot[].fC[0])");
  tree.SetAlias("c2InRot","sqrt(part.fParamInRot[].fC[2])");
  tree.SetAlias("dEdxExp","AliExternalTrackParam::BetheBlochAleph(pMC/AliPID::ParticleMass(pidCode))");
  tree.SetAlias("dEdxExpSolid","AliExternalTrackParam::BetheBlochSolid(pMC/AliPID::ParticleMass(pidCode))");
  tree.SetAlias("dEdxExpSolidL","AliExternalTrackParam::BetheBlochSolid(part.fParamMC[].fData.P()/AliPID::ParticleMass(pidCode))");
  tree.SetAlias("dEdxExpSolidL1","AliExternalTrackParam::BetheBlochSolid(part.fParamMC[Iteration$-1].fData.P()/AliPID::ParticleMass(pidCode))");
  tree.SetAlias("elossTPCIn","(part.fParamIn[159].fData.GetP()-part.fParamIn[7].fData.GetP())/part.fParamMC[1].fData.GetP()");
  tree.SetAlias("elossTPCMC","(part.fParamMC[159].fData.GetP()-part.fParamMC[7].fData.GetP())/part.fParamMC[1].fData.GetP()");
  //
  tree.SetAlias("sigmaY0","sqrt(part.fParamIn[1].fC[0])");
  tree.SetAlias("sigmaZ0","sqrt(part.fParamIn[1].fC[2])");
  tree.SetAlias("sigmaqPt0","sqrt(part.fParamIn[1].fC[14])");

  tree.SetAlias("sigmaY0Rot","sqrt(part.fParamInRot[1].fC[0])");
  tree.SetAlias("sigmaZ0Rot","sqrt(part.fParamInRot[1].fC[2])");
  tree.SetAlias("sigmaqPt0Rot","sqrt(part.fParamInRot[1].fC[14])");
}
