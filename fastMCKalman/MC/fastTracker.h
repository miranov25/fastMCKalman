#ifndef fastTracker_H
#define fastTracker_H
/*
 .L $fastMCKalman/fastMCKalman/MC/fastTracker.h

*/
// basic tracking functionality extracked https://github.com/miranov25/AliRoot/blob/master/TPC/TPCrec/AliTPCtracker.cxx#L1325

class fastTracker{
  static AliExternalTrackParam* makeSeed(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz);
  static Double_t makeC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3);     // F1
  static Double_t makeSnp(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3);     // F2
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
  // Initial approximation of the track curvature0
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
Double_t fastTracker::makeTgl(Double_t x1,Double_t y1, Double_t x2,Double_t y2,Double_t z1,Double_t z2){
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

AliExternalTrackParam * fastTracker::makeSeed(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz){
  Double_t param[5];
  Double_t c[15];
  // calculate initial param
  param[0]=xyz0[1];              
  param[1]=xyz0[2];
  param[2]=makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]);
  param[4]=makeC(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]);
  param[3]=makeTgl(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2],xyz1[2]);
  //
  Double_t f40=(makeC(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz2[0],xyz2[1])-param[4])/sy;
  Double_t f42=(makeC(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz2[0],xyz2[1])-param[4])/sy;
  Double_t f43=(makeC(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]+sy)-param[4])/sy;
  Double_t f20=(makeSnp(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz2[0],xyz2[1])-param[2])/sy;
  Double_t f22=(makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz2[0],xyz2[1])-param[2])/sy;
  Double_t f23=(makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]+sy)-param[2])/sy;
	
  Double_t f30=(makeTgl(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz0[2],xyz1[2])-param[3])/sy;
  Double_t f31=(makeTgl(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2]+sz,xyz1[2])-param[3])/sz;
  Double_t f32=(makeTgl(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz0[2],xyz1[2])-param[3])/sy;
  Double_t f34=(makeTgl(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2],xyz1[2]+sz)-param[3])/sz;
  c[0]=sy;
  c[1]=0.;       c[2]=sz;
  c[3]=f20*sy;  c[4]=0.;       c[5]=f20*sy*f20+f22*sy*f22+f23*sy*f23;
  c[6]=f30*sy;  c[7]=f31*sz;  c[8]=f30*sy*f20+f32*sy*f22;
  c[9]=f30*sy*f30+f31*sz*f31+f32*sy*f32+f34*sz*f34;
  c[10]=f40*sy; c[11]=0.; c[12]=f40*sy*f20+f42*sy*f22+f43*sy*f23;
  c[13]=f30*sy*f40+f32*sy*f42;
  c[14]=f40*sy*f40+f42*sy*f42+f43*sy*f43;

  return 0;
}



#endif