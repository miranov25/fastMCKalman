#include "fastTracker.h"
/*
  .L $fastMCKalman/fastMCKalman/MC/fastTracker.h
  .L $fastMCKalman/fastMCKalman/MC/testFastTracker.C
  testFasTrackerSimul(100000);
  testFastTrackerEval();
  WDir:
      /home2/miranov/github/fastMCKalman/data/testSeed
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
  param.GetXYZAt(1, bz, xyz0);
  param.GetXYZAt(85, bz, xyz1);
  param.GetXYZAt(250, bz, xyz2);
  //
  AliExternalTrackParam * paramSeed = fastTracker::makeSeed(xyz0,xyz1,xyz2,0.1,0.1,bz);
}

/// generate random tracks smear them and make kalman estimator
///    the error estimatro and the data should be within error
void testFasTrackerSimul(Int_t nPoints) {
  // Double_t pt=0.2,bz=0.5, tgl=0.1;
  float bz = 0.5;
  TTreeSRedirector *pcstream = new TTreeSRedirector("testSeed.root", "recreate");

  Double_t pxpypz[3] = {};
  Double_t xyz[3] = {1, 0, 0};

  Double_t cov[21] = {0};
  for (Int_t i = 0; i < nPoints; i++) {
    Float_t tgl = gRandom->Rndm();
    Float_t pt = (gRandom->Rndm() + 0.2) * 5;
    Float_t sy = (gRandom->Rndm() + 0.001) * 0.001;
    Float_t sz = (gRandom->Rndm() + 0.001) * 0.001;
    xyz[1] = gRandom->Gaus() * 5;
    xyz[2] = gRandom->Gaus() * 5;
    pxpypz[0] = pt;
    pxpypz[1] = pt * (gRandom->Gaus() * 0.05);
    pxpypz[2] = pt * tgl;
    AliExternalTrackParam param(xyz, pxpypz, cov, 1);
    param.Rotate(0);
    param.PropagateTo(250, bz);
    //
    Double_t xyz0[3], xyz1[3], xyz2[3];
    param.GetXYZAt(250, bz, xyz0);
    param.GetXYZAt(85, bz, xyz1);
    param.GetXYZAt(1, bz, xyz2);
    xyz0[1] += gRandom->Gaus() * sy;
    xyz1[1] += gRandom->Gaus() * sy;
    xyz2[1] += gRandom->Gaus() * sy;
    xyz0[2] += gRandom->Gaus() * sz;
    xyz1[2] += gRandom->Gaus() * sz;
    xyz2[2] += gRandom->Gaus() * sz;
    //
    AliExternalTrackParam *paramSeed = fastTracker::makeSeed(xyz0, xyz1, xyz2, sy, sz, bz);
    param.Rotate(paramSeed->GetAlpha());
    param.PropagateTo(paramSeed->GetX(), bz);
    (*pcstream) << "seed" <<
                "param.=" << &param <<
                "paramSeed.=" << paramSeed <<
                "sy=" << sy <<
                "sz=" << sz <<
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
  //
  // test P[0]
  tree->Draw("(paramSeed.fP[0]-param.fP[0])/sy","",""); // should be gaus with width ~ 1 - OK
  isOK= abs(1-tree->GetHistogram()->GetRMS())<4*tree->GetHistogram()->GetRMSError();
  if (isOK) {::Info("testFastTracker pull test P0","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P0","pullAnalytical- FAILED");
  }
  // test P[1]
  tree->Draw("(paramSeed.fP[1]-param.fP[1])/sz","",""); // should be gaus with width ~ 1 - OK
  isOK= abs(1-tree->GetHistogram()->GetRMS())<4*tree->GetHistogram()->GetRMSError();
  if (isOK) {::Info("testFastTracker pull test P1","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P1","pullAnalytical- FAILED");
  }
  // test P[2]
  tree->Draw("(paramSeed.fP[2]-param.fP[2])/sqrt(paramSeed.fC[5])","",""); // should be gaus with width ~ 1 - OK
  isOK= abs(1-tree->GetHistogram()->GetRMS())<4*tree->GetHistogram()->GetRMSError();
  if (isOK) {::Info("testFastTracker pull test P2","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P2","pullAnalytical- FAILED");
  }
  // test P[3]
  tree->Draw("(paramSeed.fP[3]-param.fP[3])/sqrt(paramSeed.fC[9])","",""); // should be gaus with width ~ 1 - OK
  isOK= abs(1-tree->GetHistogram()->GetRMS())<4*tree->GetHistogram()->GetRMSError();
  if (isOK) {::Info("testFastTracker pull test P3","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P3","pullAnalytical- FAILED");
  }
  // test P[4]
  tree->Draw("(paramSeed.fP[4]-param.fP[4])/sqrt(paramSeed.fC[14])","",""); // should be gaus with width ~ 1 - OK
  isOK= abs(1-tree->GetHistogram()->GetRMS())<4*tree->GetHistogram()->GetRMSError();
  if (isOK) {::Info("testFastTracker pull test P4","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P4","pullAnalytical- FAILED");
  }




}
//
// 0
// 1  2
// 3  4   5
// 6  7   8  9
// 10 11  12 13 14