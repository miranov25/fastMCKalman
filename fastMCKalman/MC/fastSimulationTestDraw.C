/*
  .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
  .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
  .L  $fastMCKalman/fastMCKalman/MC/fastSimulationTestDraw.C
  initTreeFast()
 */

void drawTimeLength() {
  //
  TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  treeFast->Draw("fParamMC[4].fData.fLength/sqrt(1+fParamMC[0].fData.fP[3]**2):fParamMC[0].fData.fP[3]:abs(fParamMC[0].fData.fP[4])","","colz");
  gPad->SaveAs("fig/trackLengthLayer4_Normed.png");
  treeFast->Draw("fParamMC[4].fData.fTime/sqrt(1+fParamMC[0].fData.fP[3]**2):fParamMC[0].fData.fP[3]:fParamMC[0].Beta()","","colz");
  gPad->SaveAs("fig/trackTimeLayer4_Normed.png");

  treeFast->Draw("fParamMC[4].fData.fTime:(fParamMC[0].Beta()*fParamMC[4].fData.fLength/2.99792458e-2):fParamMC[0].Beta()", "fParamMC[4].fData.fX>20", "colz");
}
void drawdQ() {
  treeFast->Draw("eLossLog:fParamMC[0].fData.fP[3]:fParamMC[0].Beta()", "fParamMC[].Beta()>0.05&&Iteration$<7", "colz");
}