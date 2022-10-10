/*
  Here we would like to have tests - which checks internal consistence of the fastSimulation
  Dedicated debug streamaers to be used
      gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
.L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
 */

TTree * treeSeed=0;
/// here we would lik to test internal consistency of the seeding
void testInit(){
   TFile *f = TFile::Open("fastParticle.root");
   treeSeed=(TTree*)f->Get("seedDump");
     AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
}
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

  treeFast->Draw("(partFull.fParamIn[0].fP[4]-partFull.fParamMC[0].fP[4])/sqrt(partFull.fParamIn[0].fC[14])","partFull.fStatusMaskIn[0].fData==31");
  treeFast->Draw("(part.fParamIn[0].fP[4]-part.fParamMC[0].fP[4])/sqrt(part.fParamIn[0].fC[14])","part.fStatusMaskIn[0].fData==31");


}

