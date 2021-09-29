#include "fastTracker.h"
/*
  .L fastTracker.h
 */

/// test if the track properties from seedign the same as from external track param
void testFasTracker(Float_t pt, Float_t bz, Float_t tgl){
  // Double_t pt=0.2,bz=0.5, tgl=0.1;
  Double_t pxpypz[3]={pt,0,pt*tgl};
  Double_t xyz[3]={1,0,0};
  Double_t cov[21]={0};
  //AliExternalTrackParam AliExternalTrackParam(Double_t[3] xyz, Double_t[3] pxpypz, Double_t[21] cv, Short_t sign)
  AliExternalTrackParam param(xyz,pxpypz,cov,1);
  Double_t xyz0[3], xyz1[3],xyz2[3];
  param.GetXYZatR(1, bz, xyz0);
  param.GetXYZatR(85, bz, xyz1);
  param.GetXYZatR(250, bz, xyz2);
  //
  AliExternalTrackParam * paramSeed = fastTracker::makeSeed(xyz0,xyz1,xyz2,0.1,0.1,bz);


}