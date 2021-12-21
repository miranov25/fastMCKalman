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

void testMB(){
  // Here we should make unit tests - it will work only up to some relative energy loss
  // treeSeed->Draw("log(paramSeed.P()/input.P()):log(input2.P()/input.P())","sign0>0");
  treeSeed->Draw("(seed.fP[4]-input.fP[4])/sqrt(seed.fC[14])","","");
  treeSeed->Draw("(seed.fP[4]-input.fP[4])/sqrt(seed.fC[14]):(input2.fP[4]-input.fP[4])/sqrt(seed.fC[14]):fMassMC","input.fP[4]<0&&sign0>0","colz");

}