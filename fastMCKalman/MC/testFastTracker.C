
/*
     gSystem->Load("AliSkimmedDataAnalysisMaker.so");
  //.L $fastMCKalman/fastMCKalman/MC/fastTracker.h
  .L $fastMCKalman/fastMCKalman/MC/testFastTracker.C++g
  testFasTrackerSimul(100000);
  testFastTrackerEval();
  WDir:
      /home2/miranov/github/fastMCKalman/data/testSeed
 */

#include "TFile.h"
#include "TTree.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "fastTracker.h"


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
    Float_t pt = (gRandom->Rndm() + 0.05) * 5;
    Float_t sy = (gRandom->Rndm() + 0.01) * 0.01;
    Float_t sz = (gRandom->Rndm() + 0.01) * 0.01;
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
  TTree * tree =((*pcstream) << "seed").GetTree();
  tree->SetAlias("pull0","(paramSeed.fP[0]-param.fP[0])/sqrt(paramSeed.fC[0])");
  tree->SetAlias("pull1","(paramSeed.fP[1]-param.fP[1])/sqrt(paramSeed.fC[2])");
  tree->SetAlias("pull2","(paramSeed.fP[2]-param.fP[2])/sqrt(paramSeed.fC[5])");
  tree->SetAlias("pull3","(paramSeed.fP[3]-param.fP[3])/sqrt(paramSeed.fC[9])");
  tree->SetAlias("pull4","(paramSeed.fP[4]-param.fP[4])/sqrt(paramSeed.fC[14])");
  tree->SetAlias("c02","(paramSeed.fC[3])/sqrt(paramSeed.fC[0]*paramSeed.fC[5])");
  tree->SetAlias("c13","(paramSeed.fC[7])/sqrt(paramSeed.fC[2]*paramSeed.fC[9])");
  tree->SetAlias("c04","(paramSeed.fC[10])/sqrt(paramSeed.fC[0]*paramSeed.fC[14])");
  tree->SetAlias("c24","(paramSeed.fC[12])/sqrt(paramSeed.fC[5]*paramSeed.fC[14])");
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
  // corelation tests
  {
    TF1 *f1 = new TF1("f1", "[0]*x");
    tree->Draw("(paramSeed.fP[1]-param.fP[1])*(paramSeed.fP[3]-param.fP[3]):paramSeed.fC[7]", "sz>0.0005", "prof");
    tree->GetHistogram()->Fit("f1"); // ~ 1 - looks OK
    if (TMath::Abs(f1->GetParameter(0) - 1) < 0.05) {
      ::Info("testFastTracker C(1,3)", "pullAnalytical - OK");
    } else { ::Error("testFastTracker pull test C(1,2)", "pullAnalytical- FAILED"); };
    //
    tree->Draw("(paramSeed.fP[0]-param.fP[0])*(paramSeed.fP[2]-param.fP[2]):paramSeed.fC[3]", "sz>0.0005", "prof");
    tree->GetHistogram()->Fit("f1"); //  ~1 - looks OK
    if (TMath::Abs(f1->GetParameter(0) - 1) < 0.05) {
      ::Info("testFastTracker C(0,2)", "pullAnalytical - OK");
    } else { ::Error("testFastTracker pull test C(0,2)", "pullAnalytical- FAILED"); };
    //
    tree->Draw("(paramSeed.fP[2]-param.fP[2])*(paramSeed.fP[4]-param.fP[4]):paramSeed.fC[12]", "sz>0.0005", "prof");
    tree->GetHistogram()->Fit("f1"); //  ~1 - looks OK
    if (TMath::Abs(f1->GetParameter(0) - 1) < 0.05) {
      ::Info("testFastTracker C(2,4)", "pullAnalytical - OK");
    } else { ::Error("testFastTracker pull test C(2,4)", "pullAnalytical- FAILED"); };
    tree->Draw("(paramSeed.fP[0]-param.fP[0])*(paramSeed.fP[4]-param.fP[4]):paramSeed.fC[10]", "sz>0.0005", "prof");
    tree->GetHistogram()->Fit("f1"); //  ~  looks OK
    if (TMath::Abs(f1->GetParameter(0) - 1) < 0.05) {
      ::Info("testFastTracker C(0,4)", "pullAnalytical - OK");
    } else { ::Error("testFastTracker pull test C(0,4)", "pullAnalytical- FAILED"); };
  }



}
//
// 0
// 1  2
// 3  4   5
// 6  7   8  9
// 10 11  12 13 14