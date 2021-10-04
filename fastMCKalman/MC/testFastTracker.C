#include "fastTracker.h"
/*
  .L $fastMCKalman/fastMCKalman/MC/fastTracker.h
  .L $fastMCKalman/fastMCKalman/MC/testFastTracker.C
  testFasTrackerSimul(1000);

 */
TTree * tree = 0;

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
void testFasTrackerSimul(Int_t nPoints){
  // Double_t pt=0.2,bz=0.5, tgl=0.1;
  float bz=0.5;
   TTreeSRedirector *pcstream = new TTreeSRedirector("testSeed.root","recreate");

  Double_t pxpypz[3]={};
  Double_t xyz[3]={1,0,0};
  Double_t cov[21]={0};
  for (Int_t i=0; i<nPoints; i++){
     Float_t tgl=gRandom->Rndm();
     Float_t pt=(gRandom->Rndm()+0.05)*5;
     Float_t sy=(gRandom->Rndm()+0.1)*0.01;
     Float_t sz=(gRandom->Rndm()+0.1)*0.01;
     pxpypz[0]=pt;
     pxpypz[1]=0;
     pxpypz[2]=pt*tgl;
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
       "\n";
  }
  delete pcstream;
}
///
/// check pulls of the input parttice and reconstructed track
void testFastTrackerEval(){
  TFile *f = new TFile("testSeed.root");
  tree = (TTree *) f->Get("seed");
  Bool_t isOK=0;
  tree->SetAlias("errP3A","(sqrt(2.)*sz/(245-85))");
  // test 0
  tree->Draw("(paramSeed.fP[3]-param.fP[3])/(sqrt(2.)*sz/(245-85))","",""); // should be gaus with width ~ 1 - OK
  isOK= abs(1-tree->GetHistogram()->GetRMS())<4*tree->GetHistogram()->GetRMSError();
  if (isOK) {
    ::Info("testFastTracker P3 Test1","pullAnalytical");
  }else{
    ::Error("testFastTracker P3 Test1 ","pullAnalytical");
  }
  tree->Draw("(paramSeed.fP[3]-param.fP[3])/sqrt(paramSeed.fC[9])","","");

}
