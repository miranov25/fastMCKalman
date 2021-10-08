#ifndef fastTracker_H
#define fastTracker_H
/*
 .L $fastMCKalman/fastMCKalman/MC/fastTracker.h
 code based on the AliTPCtracker.
*/
// basic tracking functionality extracked
// https://github.com/miranov25/AliRoot/blob/master/TPC/TPCrec/AliTPCtracker.cxx#L1325
// https://github.com/miranov25/AliRoot/blob/2f249b4d1cff4e54e0b17f2b23ddb4d9b5e4bd9e/STEER/ESD/AliTrackerBase.cxx#L673
#include "AliExternalTrackParam.h"

class fastTracker{
public:
  static AliExternalTrackParam* makeSeed(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz);
  static Double_t makeC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3);     // F1
  static Double_t makeSnp(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3);     // F2
  static Double_t makeTgln(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2,Double_t c);   // F3n
  static Double_t makeTgl(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2);   // F3
};


//_____________________________________________________________________________
Double_t fastTracker::makeC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u;
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t fastTracker::makeSnp(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approximation of the track snp at position x1
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u; 
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  x0+=x1;
  x0*=c2;  
  return x0;
}

//_____________________________________________________________________________
Double_t fastTracker::makeTgln(Double_t x1,Double_t y1, Double_t x2,Double_t y2,Double_t z1,Double_t z2,Double_t c){
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  Double_t d  =  TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  Double_t   angle2    = asinf(d*c*0.5);

  angle2  = (z1-z2)*c/(angle2*2.);
  return angle2;
  //return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

Double_t fastTracker::makeTgl(Double_t x1,Double_t y1, Double_t x2,Double_t y2,Double_t z1,Double_t z2){
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

AliExternalTrackParam * fastTracker::makeSeed(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz){
  Double_t sy2=sy*sy;
  Double_t sz2=sz*sz;
  Double_t param[5];
  Double_t c[15];
  // calculate initial param
  param[0]=xyz0[1];              
  param[1]=xyz0[2];
  param[2]=makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]);
  param[4]=makeC(xyz2[0],xyz2[1],xyz1[0],xyz1[1],xyz0[0],xyz0[1]);
  param[3]=makeTgln(xyz2[0],xyz2[1],xyz1[0],xyz1[1],xyz2[2],xyz1[2],param[4]);
  //
  Double_t f40=(makeC(xyz2[0],xyz2[1]+sy,xyz1[0],xyz1[1],xyz0[0],xyz0[1])-param[4])/sy;
  Double_t f42=(makeC(xyz2[0],xyz2[1],xyz1[0],xyz1[1]+sy,xyz0[0],xyz0[1])-param[4])/sy;
  Double_t f43=(makeC(xyz2[0],xyz2[1],xyz1[0],xyz1[1],xyz0[0],xyz0[1]+sy)-param[4])/sy;
  //
  Double_t f20=(makeSnp(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz2[0],xyz2[1])-param[2])/sy;
  Double_t f22=(makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz2[0],xyz2[1])-param[2])/sy;
  Double_t f23=(makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]+sy)-param[2])/sy;
  //
  //  makeTgln(xyz2[0],xyz2[1],xyz1[0],xyz1[1],xyz2[2],xyz1[2],param[4]);
  Double_t f30=(makeTgln(xyz2[0],xyz2[1]+sy,xyz1[0],xyz1[1],xyz2[2],xyz1[2],param[4])-param[3])/sy;
  Double_t f31=(makeTgln(xyz2[0],xyz2[1],xyz1[0],xyz1[1],xyz2[2]+sz,xyz1[2],param[4])-param[3])/sz;
  Double_t f32=(makeTgln(xyz2[0],xyz2[1],xyz1[0],xyz1[1]+sy,xyz2[2],xyz1[2],param[4])-param[3])/sy;
  Double_t f34=(makeTgln(xyz2[0],xyz2[1],xyz1[0],xyz1[1],xyz2[2],xyz1[2]+sz,param[4])-param[3])/sz;
  c[0]=sy2;
  c[1]=0.;       c[2]=sz2;
  c[3]=f20*sy2;   c[4]=0.;       c[5]=f20*sy2*f20+f22*sy2*f22+f23*sy2*f23;
  c[6]=f30*sy2;  c[7]=f31*sz2;  c[8]=f30*sy2*f20+f32*sy2*f22;
  c[9]=f30*sy2*f30+f31*sz2*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy2; c[11]=0.; c[12]=f40*sy2*f20+f42*sy2*f22+f43*sy2*f23;
  c[13]=f30*sy2*f40+f32*sy2*f42;
  c[14]=f40*sy2*f40+f42*sy2*f42+f43*sy2*f43;

  if (TMath::Abs(bz)>kAlmost0Field) {
    c[14]/=(bz*kB2C)*(bz*kB2C);
    param[4]/=(bz*kB2C); // transform to 1/pt
  }
  else { // assign 0.6 GeV pT
    const double kq2pt = 1./0.6;
    param[4] = kq2pt;
    c[14] = (0.5*0.5)*kq2pt;
  }

  //
  //AliExternalTrackParam *param=new AliExternalTrackParam(
  //AliExternalTrackParam AliExternalTrackParam(Double_t x, Double_t alpha, const Double_t[5] param, const Double_t[15] covar)
  AliExternalTrackParam *extParam=new AliExternalTrackParam(xyz0[0],0,param,c);
  return extParam;
}



#endif