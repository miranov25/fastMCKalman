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
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "fastTracker.h"
#include "fastSimulation.h"




extern TChain * treeSeed;
extern TChain * treeFast;

/// to check momentum bias in the seeds

/// Test pull of seeds
/// the seeding pulls at the big snp are underestimated - because of rotation - we call test only for smaller |fP[2]|<0.6 further from 1
void testPullsSeed() {
  // Here we should make unit tests - it will work only up to some relative energy loss
  // treeSeed->Draw("log(paramSeed.P()/input.P()):log(input2.P()/input.P())","sign0>0");

  //
  TF1 *mygauss = new TF1("mygauss", "gaus");
  for (int version=0; version<=1; version++) {
    for (int iPar = 0; iPar <= 4; iPar++) {
      treeSeed->Draw(Form("(seed.fP[%d]-input.fP[%d])/sqrt(seed.fC[%d])>>his(100,-6,6)", iPar, iPar, AliExternalTrackParam::GetIndex(iPar, iPar)),
                     Form("version==%d&&abs(seed.fP[2])<0.6",version), "");
      treeSeed->GetHistogram()->Fit("mygauss", "q");
      bool isOK = abs(1 - mygauss->GetParameter(2)) < 5 * mygauss->GetParError(2);
      float rms=treeSeed->GetHistogram()->GetRMS();
      if (isOK) {
        ::Info(Form("testFastTracker seed pull test P%d - version %d",iPar,version), "pullAnalytical - OK - %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      } else {
        ::Error(Form("testFastTracker seed pull test P%d - version %d",iPar,version), "pullAnalytical- FAILED- %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      }
    }
  }
}

void testLooperSmooth(){
  // checking the X position
  treeFast->Draw("gyMC:gxMC:(part.fParamMC[].fX==part.fParamMC[Iteration$-2].fX)","Iteration$>2","colz",10);
}


void testPulls() {
  TF1 *mygauss = new TF1("mygauss", "gaus");
  int isOK=fastParticle::kTrackUpdate |fastParticle::kTrackChi2;
  for (int version=0; version<=1; version++) {
    for (int iPar = 0; iPar <= 4; iPar++) {
      std::string sv = "";
      if (version==1) sv = "Full";
      treeFast->Draw(Form("(part%s.fParamIn[0].fP[%d]-part%s.fParamMC[0].fP[%d])/sqrt(part%s.fParamIn[0].fC[%d])>>his(100,-6,6)",sv.c_str(), iPar, sv.c_str(), iPar, sv.c_str(), AliExternalTrackParam::GetIndex(iPar, iPar)),
                    Form("(part%s.fStatusMaskIn[0].fData&%d)==%d&&abs(part%s.fParamIn[0].fP[2])<0.7",sv.c_str(),isOK,isOK,sv.c_str()), "");
      treeFast->GetHistogram()->Fit("mygauss", "q");
      bool isOK = abs(1 - mygauss->GetParameter(2)) < 5 * mygauss->GetParError(2);
      float rms=treeFast->GetHistogram()->GetRMS();
      if (isOK) {
        ::Info(Form("testFastTracker part%s reco pull test P%d ",sv.c_str(),iPar), "pullAnalytical - OK - %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      } else {
        ::Error(Form("testFastTracker part%s reco pull test P%d",sv.c_str(), iPar), "pullAnalytical- FAILED- %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      }
    }
  }
}

void test_compare(){

}

/*
/// this is pseudocode to be used in propagate to mirror to account or the material along the arc
/// hee we have the same logical problem as in correct for mean material - factor sqrt(2) in angular spread
/// Here we assume the parameters defining MS are not changing along trajectory
/// \param crossLength
/// \param xx0
/// \param beta2
/// \param p2
/// \param nSteps
/// \return
 std::vector<float>  makeIntegralCovar(float crossLength, float xx0, float beta2, float p2, float nSteps){
  float xOverX0=crossLength*xx0/nSteps;
  std::vector<float> covarLinear(10);   0-sy2, 2-sz2, 5-sphi2, 9-stheta2 as in AliExternalTrackParam
  Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0); // smearing angle per unit step
  float dL=crossLength/nSteps;
  //
  for (int i=0; i<nSteps;i++ ){
     covarLinear[0]+=(theta2*dL*dL);    //
     covarLinear[2]+=(theta2*dL*dL);
     covarLinear[5]+=(theta2);          //
     covarLinear[9]+=(theta2);          //
     //
     covarLinear[3]+=theta2*dL;
     covarLinear[7]+=theta2*dL;
  }
  return covarLinear;
}

void addCovariance(AliExternalTrackParam4D track){
  ///
  float sigmaPhi2= 0 ;  /// ???
  float sigmaTheta2=0 ; /// ???
  float sigmaCurv2=0;    /// ???
  fC[0]+=sigmaPhi2*crossLength;
  fC[2]+=sigmaTheta2*crossLength;
  fC[3]+=???
  fC[7]+=???
  dSigmaPhi=sigmaCurv2*crossLength;  /// to check
  //
}
*/

/// add canvas with views for MC and data
/// + proint particle properties
/// first
//void drawNextProblem(int counter){
//  treeFast->Draw("gyMC:gxMC:(partFull.fStatusMaskIn.fData==0x1)|Iteration$==0","","colz",1,counter++);
//}

void testDrawProblems(){
  // checking the X position
  // treeFast->Draw("gyMC:gxMC:(part.fParamMC[].fX==part.fParamMC[Iteration$-2].fX)","Iteration$>2","colz",10);
  //treeFast->Draw("gyMC:gxMC:partFull.fStatusMaskIn.fData==0x1","Iteration$>2","colz",100);
  //
  treeFast->Draw(">>ProblemRot","Sum$(partFull.fStatusMaskIn.fData==0x1)","entrylist");
  TEntryList* problemList0x1 =(TEntryList*)gDirectory->Get("ProblemRot");
  treeFast->SetEntryList(problemList0x1);
  int counter=0;
  treeFast->Draw("gyMC:gxMC:(partFull.fStatusMaskIn.fData==0x1)|((Iteration$==0)*2)|((Iteration$==partFull.fFirstIndex)*3)","","colz",1,counter++);

}