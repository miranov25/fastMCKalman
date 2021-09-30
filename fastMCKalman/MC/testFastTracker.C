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

/// generate random tracks smear them and make kalman estimator
///    the error estimatro and the data should be within error
void testFasTracker(Int_t nPoints){
  // Double_t pt=0.2,bz=0.5, tgl=0.1;
   TTreeSRedirector *pcstream = new TTreeSRedirector("testSeed.root","recreate")

  Double_t pxpypz[3]={};
  Double_t xyz[3]={1,0,0};
  Double_t cov[21]={0};
  for (Int_t i=0; i<nPoints; i++){
     Float_t tgl=gRandom->Rndm();
     Float_t pt=(gRandom->Rndm()+0.05)*5;
     Float_t sy=(gRandom->Rndm()+0.1)*0.1;
     Float_t sz=(gRandom->Rndm()+0.1)*0.1;
     pxpypz[0]=pt;
     pxpypz[1]=pt*tgl;
     AliExternalTrackParam param(xyz,pxpypz,cov,1);
     param.Rotate(0);
     Double_t xyz0[3], xyz1[3],xyz2[3];
    param.GetXYZatR(1, bz, xyz0);
    param.GetXYZatR(85, bz, xyz1);
    param.GetXYZatR(250, bz, xyz2);
    xyz0[1]+=gRandom->Gaus()*sy;
    xyz1[1]+=gRandom->Gaus()*sy;
    xyz2[1]+=gRandom->Gaus()*sy;
    xyz0[2]+=gRandom->Gaus()*sz;
    xyz1[2]+=gRandom->Gaus()*sz;
    xyz2[2]+=gRandom->Gaus()*sz;
    //
    AliExternalTrackParam * paramSeed = fastTracker::makeSeed(xyz0,xyz1,xyz2,sy,sz,bz);
    (*pcstream)<<"seed"<<
      "param.="<<&param<<
       "paramSeed.="<<paramSeed<<
       "sy="<<sy<<
       "sz="<<sz<<
       "n";
  }
  delete pcstream;
}

