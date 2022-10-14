/*
  Here we would like to have tests - which checks internal consistence of the fastSimulation
  Dedicated debug streamaers to be used
      gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
.L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
     .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    .L $fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g
    initTreeFast();
 */

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "fastTracker.h"




extern TChain * treeSeed;
extern TChain * treeFast;

/// to check momentum bias in the seeds

/// Test pull of seeds
/// TODO add automatic alarm - and fix eventual problems
void testMB(){
  // Here we should make unit tests - it will work only up to some relative energy loss
  // treeSeed->Draw("log(paramSeed.P()/input.P()):log(input2.P()/input.P())","sign0>0");
  treeSeed->Draw("(seed.fP[4]-input.fP[4])/sqrt(seed.fC[14])","","");
  treeSeed->Draw("(seed.fP[4]-input.fP[4])/sqrt(seed.fC[14]):(input2.fP[4]-input.fP[4])/sqrt(seed.fC[14]):fMassMC","input.fP[4]<0&&sign0>0","colz");
}

void testLooperSmooth(){
  // checking the X position
  treeFast->Draw("gyMC:gxMC:(part.fParamMC[].fX==part.fParamMC[Iteration$-2].fX)","Iteration$>2","colz",10);
}

void testPulls(){
  /// define aliases  - to put later to the main code
  /// fit pulls  1D and as function of qPt
  //TFile *f = TFile::Open("fastParticle.root");
  //treeFast=(TTree*)f->Get("fastPart");
  Bool_t isOK=0;

  //Full Reco
  TH1F *myh = new TH1F("myh","myh",100,-6,6);
  TF1 *mygauss = new TF1("mygauss","gaus");
  treeFast->Draw("(partFull.fParamIn[0].fP[4]-partFull.fParamMC[0].fP[4])/sqrt(partFull.fParamIn[0].fC[14])>>myh","partFull.fStatusMaskIn[0].fData==31");
  myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
  if (isOK) {::Info("testFastTracker full reco pull test P4","pullAnalytical - OK");
  }else{::Error("testFastTracker full reco pull test P4","pullAnalytical- FAILED");
  }
  
  treeFast->Draw("(partFull.fParamIn[0].fP[3]-partFull.fParamMC[0].fP[3])/sqrt(partFull.fParamIn[0].fC[9])>>myh","partFull.fStatusMaskIn[0].fData==31");
  myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker full reco pull test P3","pullAnalytical - OK");
  }else{::Error("testFastTracker full reco pull test P3","pullAnalytical- FAILED");
  }
  
  treeFast->Draw("(partFull.fParamIn[0].fP[2]-partFull.fParamMC[0].fP[2])/sqrt(partFull.fParamIn[0].fC[5])>>myh","partFull.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker full reco pull test P2","pullAnalytical - OK");
  }else{::Error("testFastTracker full reco pull test P2","pullAnalytical- FAILED");
  }
  treeFast->Draw("(partFull.fParamIn[0].fP[1]-partFull.fParamMC[0].fP[1])/sqrt(partFull.fParamIn[0].fC[2])>>myh","partFull.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker full reco pull test P1","pullAnalytical - OK");
  }else{::Error("testFastTracker full reco pull test P1","pullAnalytical- FAILED");
  }
  treeFast->Draw("(partFull.fParamIn[0].fP[0]-partFull.fParamMC[0].fP[0])/sqrt(partFull.fParamIn[0].fC[0])>>myh","partFull.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker full reco pull test P0","pullAnalytical - OK");
  }else{::Error("testFastTracker full reco pull test P0","pullAnalytical- FAILED");
  }
 
  ////Leg Reco
  treeFast->Draw("(part.fParamIn[0].fP[4]-part.fParamMC[0].fP[4])/sqrt(part.fParamIn[0].fC[14])>>myh","part.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
  if (isOK) {::Info("testFastTracker pull test P4","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P4","pullAnalytical- FAILED");
  }

  treeFast->Draw("(part.fParamIn[0].fP[3]-part.fParamMC[0].fP[3])/sqrt(part.fParamIn[0].fC[9])>>myh","part.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker pull test P3","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P3","pullAnalytical- FAILED");
  }
  
  treeFast->Draw("(part.fParamIn[0].fP[2]-part.fParamMC[0].fP[2])/sqrt(part.fParamIn[0].fC[5])>>myh","part.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker pull test P2","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P2","pullAnalytical- FAILED");
  }
  treeFast->Draw("(part.fParamIn[0].fP[1]-part.fParamMC[0].fP[1])/sqrt(part.fParamIn[0].fC[2])>>myh","part.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker pull test P1","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P1","pullAnalytical- FAILED");
  }
  treeFast->Draw("(part.fParamIn[0].fP[0]-part.fParamMC[0].fP[0])/sqrt(part.fParamIn[0].fC[0])>>myh","part.fStatusMaskIn[0].fData==31");
   myh->Fit("mygauss","q");
  isOK= abs(1-mygauss->GetParameter(2))<4*mygauss->GetParError(2);
   if (isOK) {::Info("testFastTracker pull test P0","pullAnalytical - OK");
  }else{::Error("testFastTracker pull test P0","pullAnalytical- FAILED");
  
}

}

void testDebugFailure(){
  treeFast->Draw("partFull.fStatusMaskIn.fData","(partFull.fStatusMaskIn.fData)>0&&partFull.fStatusMaskIn.fData!=31");
  treeFast->Draw("partFull.fParamIn[Iteration$].fP[2]","(partFull.fStatusMaskIn.fData)==19&&partFull.fStatusMaskIn.fData!=31");

}

