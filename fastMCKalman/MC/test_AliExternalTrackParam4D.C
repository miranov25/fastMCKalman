/*
  .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    .L $fastMCKalman/fastMCKalman/MC/test_AliExternalTrackParam4D.C
  initTreeFast()
*/

/// unit test dPdx formula
void testdPdx0(){
  const Float_t maxRMS=0.001;
  //treeUnit0->Draw("AliExternalTrackParam4D::dPdx(paramIn.P(),mass):(AliExternalTrackParam4D::dPdxEuler(paramIn.P(),mass,0.01*xTimesRho)/(0.01*xTimesRho))","statusStepRK==1&&mass>0.1&&paramIn.P()/mass<2","colz")
  treeUnit0->Draw("AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*(0.01*xTimesRho)-(AliExternalTrackParam4D::dPdxEuler(paramIn.P(),mass,0.01*xTimesRho))","statusStepRK==1","goff");
  if (treeUnit0->GetHistogram()->GetRMS()>maxRMS){
    ::Error("testdPdx0","FAILED RMS %f>%f", treeUnit0->GetHistogram()->GetRMS(),maxRMS);
  }else{
    ::Info("testdPdx0","OK RMS %f<%f", treeUnit0->GetHistogram()->GetRMS(),maxRMS);
  }
}
/// fit dPdxCorrection formula and test the content
void fitdPdxScaling(){
  TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  TF1 f1("f1","x*(1+[0]*x+[1]*x*x+[2]*x*x*x)");
  f1.SetParName(0, "dX2");
  f1.SetParName(1, "dX3");
  f1.SetParName(2, "dX4");
  gStyle->SetOptFit(1);
  treeUnit0->Draw("((paramStepRK.P()-paramIn.P())/paramIn.P()):(AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*xTimesRho)/paramIn.P()>>his(300,-0.3,0.2,300,-0.8,0.2)","statusStepRK==1&&status0==1","colz");
  treeUnit0->GetHistogram()->Fit("f1","","",-0.25,.2);
  gPad->SaveAs("fig/fitdPdxScaling.png");
  //
  treeUnit0->SetAlias("dPdxR0","((paramStepRK.P()-paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxR","(AliExternalTrackParam4D::dPdx(paramIn.P(),mass+0)*xTimesRho)/paramIn.P()");
  treeUnit0->SetAlias("dPdxRFit",Form("dPdxR*(1+%f*dPdxR+%f*dPdxR**2+%f*dPdxR**3)",f1.GetParameter(0), f1.GetParameter(1),f1.GetParameter(2)));
  treeUnit0->Draw("dPdxR0-dPdxRFit:dPdxR0:mass","statusStepRK==1&&status0==1&&dPdxR0>-0.5","colz");
  gPad->SaveAs("fig/fitdPdxScalingDiff.png");
  // make test of the AliExternalTrackParam4D::dPdxCorr(paramIn.P(),mass,xTimesRho)/paramIn.P()
  treeUnit0->Draw("dPdxR0-AliExternalTrackParam4D::dPdxCorr(paramIn.P(),mass,xTimesRho)/paramIn.P():dPdxR0:mass","statusStepRK==1&&status0==1&&dPdxR0>-0.5","colz");
  gPad->SaveAs("fig/fitdPdxCorrDiff.png");
}

void testRungeKutaDraw(){
    TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  // Euler propagation
  treeUnit0->Draw("log(paramStep.P()/param0.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&status0==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramEuler_PLoss.png");
  // Runge-Kuta in Energy
  treeUnit0->Draw("log(paramStepRK.P()/paramRK.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&status0==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramRK_PLoss.png");
  // Runge-Kuta in momenta
  treeUnit0->Draw("log(paramStepRK.P()/paramRKP.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&status0==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramRKP_PLoss.png");
  //
  treeUnit0->Draw("log(paramStepRK.P()/paramT4.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&status0==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramT4_PLoss.png");

}



void drawRungeKuttaTestMS() {
  //
     TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  treeUnit0->Draw("(paramStep.fC[9]-param0.fC[9])/(paramStep.fC[9]-paramIn.fC[9]):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
   gPad->SaveAs("fig/unitTest_paramEuler_CovTtgl.png");
  treeUnit0->Draw("(paramStep.fC[9]-paramRK.fC[9])/(paramStep.fC[9]-paramIn.fC[9]):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
  gPad->SaveAs("fig/unitTest_paramRK_CovTtgl.png");
  //

  treeUnit0->Draw("(paramStep.fC[9]-paramRK2.fC[9])/(paramStep.fC[9]-paramIn.fC[9]):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
  gPad->SaveAs("fig/unitTest_paramRK2_CovTtgl.png");
  //

  treeUnit0->Draw("log(paramStep.fC[9]-paramT4.fC[9]):log(paramStep.fC[9]-paramIn.fC[9]):mass","xOverX0>0.01&&statusStep==1&&statusT4==1&&paramIn.P()<1&&paramStep.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramT4_CovTtgl.png");

}