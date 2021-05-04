/*
  .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    .L $fastMCKalman/fastMCKalman/MC/test_AliExternalTrackParam4D.C
  initTreeFast()
  setFun();
  setUnitAlias();
*/
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TObject.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"

extern TChain * treeUnit0;

TF1 fT4("fT4","x*(1+[0]*x+[1]*x*x+[2]*x*x*x)");

/// taylor expansion  - expressing function using first derivative of function
///  f0=0
///  f'(0)=1
///  second, third and 4 derivative parameters
void setFun() {
  fT4.SetParName(0, "dX2");
  fT4.SetParName(1, "dX3");
  fT4.SetParName(2, "dX4");
}
//
void unitTestFormula(){
  TF2 fEulerPion02("fDeltaT42pion","log(AliExternalTrackParam4D::dPdxEulerStep(x,0.139,y,0.002))",0.02,0.2,0.001,0.1);
  TF2 fCorrT402Pion02("fCorr42pion","log(AliExternalTrackParam4D::dPdxCorrT42(x,0.139,y))",0.02,0.2,0.001,0.1);
  TF2 fDeltaT42Pion02("fDeltaT42pion","(AliExternalTrackParam4D::dPdxCorrT42(x,0.139,y)-AliExternalTrackParam4D::dPdxEulerStep(x,0.139,y,0.002))/(AliExternalTrackParam4D::dPdxCorrT42(x,0.139,y)+AliExternalTrackParam4D::dPdxEulerStep(x,0.139,y,0.002))",0.02,0.2,0.001,0.1);
  //
  TCanvas * c = new TCanvas("unitTestFormula","unitTestFormula",800,600);
  c->Divide(2,2);
  c->cd(1);
  fEulerPion02.Draw("colz");
  c->cd(2);
  fCorrT402Pion02.Draw("colz");
  c->cd(3);
  fDeltaT42Pion02.Draw("colz");
}

//
void setUnitAlias() {
  treeUnit0->SetAlias("stepOK", "statusStep==1&&statusStepRK==1");
  treeUnit0->SetAlias("sign", "-1.+2.*(Entry$%2)");
  treeUnit0->SetAlias("xTimesRhoS", "xTimesRho*sign");
  //
  treeUnit0->SetAlias("dPdxConstLR", "log((AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*xTimesRho+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxEuler01LR", "log((AliExternalTrackParam4D::dPdxEulerStep(paramIn.P(),mass,xTimesRho,0.01)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxEuler002LR", "log((AliExternalTrackParam4D::dPdxEulerStep(paramIn.P(),mass,xTimesRho,0.002)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxCorrT4LR", "log((AliExternalTrackParam4D::dPdxCorrT4(paramIn.P(),mass,xTimesRho)+paramIn.P())/paramIn.P())");
  // log (Pout/Pin) - based on different approximations
  /// constants, steps (-0.01,0.002,0.001), taylor 4 approximation
  treeUnit0->SetAlias("dPdxConstLRS", "log((AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*xTimesRho*sign+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxEuler01LRS", "log((AliExternalTrackParam4D::dPdxEulerStep(paramIn.P(),mass,xTimesRho*sign,0.01)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxEuler002LRS", "log((AliExternalTrackParam4D::dPdxEulerStep(paramIn.P(),mass,xTimesRho*sign,0.002)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxEuler001LRS", "log((AliExternalTrackParam4D::dPdxEulerStep(paramIn.P(),mass,xTimesRho*sign,0.001)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxCorrT4LRS", "log((AliExternalTrackParam4D::dPdxCorrT4(paramIn.P(),mass,xTimesRho*sign)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxCorrT42LRS", "log((AliExternalTrackParam4D::dPdxCorrT42(paramIn.P(),mass,xTimesRho*sign)+paramIn.P())/paramIn.P())");
  // define upper -lower band -assuming 1 % error
  treeUnit0->SetAlias("dPdxCorrT42LRSD01", "log((AliExternalTrackParam4D::dPdxCorrT42(paramIn.P()*1.01,mass,xTimesRho*sign)+paramIn.P())/paramIn.P())-log((AliExternalTrackParam4D::dPdxCorrT42(paramIn.P()*0.99,mass,xTimesRho*sign)+paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxCorrT42LRSD01", "abs(log((AliExternalTrackParam4D::dPdxCorrT42(paramIn.P()*(1+0.01/paramIn.Beta()),mass,xTimesRho*sign)+paramIn.P())/paramIn.P())-log((AliExternalTrackParam4D::dPdxCorrT42(paramIn.P()*(1-0.01/paramIn.Beta()),mass,xTimesRho*sign)+paramIn.P())/paramIn.P()))");
  //
  treeUnit0->SetAlias("dPdxStepLR", "log(paramStep.P()/paramIn.P())");
  treeUnit0->SetAlias("dPdxStepRKLR", "log(paramStepRK.P()/paramIn.P())");
}

/// show correlation and delta
///  TODO - this is obsolete
void drawUnitTest(){
  TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  treeUnit0->SetMarkerColor(2);
  treeUnit0->Draw("dPdxStepRKLR:dPdxStepRKLR","stepOK&&abs(dPdxStepLR)<0.6","");
  treeUnit0->SetMarkerColor(1);
  treeUnit0->Draw("dPdxConstLR:dPdxStepRKLR","stepOK&&abs(dPdxStepLR)<0.6","same");
  treeUnit0->SetMarkerColor(4);
  treeUnit0->Draw("dPdxEulerLR:dPdxStepRKLR","stepOK&&abs(dPdxStepLR)<0.6","same");
  gPad->Draw();
}



/// unit test dPdx formula
///
void test_dPdx0(){
  const Float_t maxRMS=0.001;
  //treeUnit0->Draw("AliExternalTrackParam4D::dPdx(paramIn.P(),mass):(AliExternalTrackParam4D::dPdxEuler(paramIn.P(),mass,0.01*xTimesRho)/(0.01*xTimesRho))","statusStepRK==1&&mass>0.1&&paramIn.P()/mass<2","colz")
  treeUnit0->Draw("AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*(0.01*xTimesRho)-(AliExternalTrackParam4D::dPdxEuler(paramIn.P(),mass,0.01*xTimesRho))","statusStepRK==1","goff");
  if (treeUnit0->GetHistogram()->GetRMS()>maxRMS){
    ::Error("testdPdx0","FAILED RMS %f>%f", treeUnit0->GetHistogram()->GetRMS(),maxRMS);
  }else{
    ::Info("testdPdx0","OK RMS %f<%f", treeUnit0->GetHistogram()->GetRMS(),maxRMS);
  }
}

/// fit taylor correction for derivative separatelly in pus and minus direction
///  1.) show qiality of the
void fitT4Correction(){
  TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  // correction vs 1 taylor expectation
  gSystem->GetFromPipe("mkdir -p fig/fitT4Correction");
  treeUnit0->Draw("dPdxEuler001LRS:dPdxConstLRS:abs(xTimesRho)","abs(xTimesRho)<0.2&&dPdxConstLRS>-0.5&&dPdxConstLRS<3","colz");
  // correction vs taylor expactation 4
  gPad->SaveAs("fig/fitT4Correction/dPdxEuler001LRS_dPdxConstLRS.png");
  treeUnit0->Draw("dPdxEuler001LRS-dPdxCorrT4LRS:dPdxConstLRS:abs(xTimesRho)","abs(xTimesRho)<0.2&&dPdxConstLRS>-0.5&&dPdxConstLRS<3","colz");
  gPad->SaveAs("fig/fitT4Correction/dPdxEuler001LRS_dPdxConstLRS.png");
  //correction vs taylor expactation 4 - 2 sides
  treeUnit0->Draw("dPdxEuler001LRS-dPdxCorrT42LRS:dPdxConstLRS:abs(xTimesRho)","abs(xTimesRho)<0.2&&dPdxConstLRS>-0.5&&dPdxConstLRS<3","colz");
  gPad->SaveAs("fig/fitT4Correction/dPdxEuler001LRS_dPdxCorrT42LRS.png");
  //
  // make fit of the taylow exansion
  treeUnit0->Draw("dPdxEuler001LRS:dPdxConstLRS>>his(200,-0.5,3,200,-2,2)","abs(xTimesRho)<0.2&&dPdxEuler001LRS!=0","colz");
  // fit both sides
  treeUnit0->GetHistogram()->Fit(&fT4,"","",-0.3,0.55);
  printf("Range both: %f,%f,%f\n", fT4.GetParameter(0), fT4.GetParameter(1),fT4.GetParameter(2));
  treeUnit0->SetAlias("dPdxRFitTB",Form("dPdxConstLRS*(1+%f*dPdxConstLRS+%f*dPdxConstLRS**2+%f*dPdxConstLRS**3)",fT4.GetParameter(0), fT4.GetParameter(1),fT4.GetParameter(2)));
  // fit left region
  treeUnit0->GetHistogram()->Fit(&fT4,"","",-0.3,0.1);
  printf("Range minus: %f,%f,%f\n", fT4.GetParameter(0), fT4.GetParameter(1),fT4.GetParameter(2));
  treeUnit0->SetAlias("dPdxRFitT4M",Form("dPdxConstLRS*(1+%f*dPdxConstLRS+%f*dPdxConstLRS**2+%f*dPdxConstLRS**3)",fT4.GetParameter(0), fT4.GetParameter(1),fT4.GetParameter(2)));
  // fit right region
  treeUnit0->GetHistogram()->Fit(&fT4,"","",-0.05,2);
  printf("Renage plus: %f,%f,%f\n", fT4.GetParameter(0), fT4.GetParameter(1),fT4.GetParameter(2));
  treeUnit0->SetAlias("dPdxRFitT4P",Form("dPdxConstLRS*(1+%f*dPdxConstLRS+%f*dPdxConstLRS**2+%f*dPdxConstLRS**3)",fT4.GetParameter(0), fT4.GetParameter(1),fT4.GetParameter(2)));
  treeUnit0->GetHistogram()->Fit(&fT4,"","",-0.05,2);
  treeUnit0->SetAlias("dPdxRFitT4","(xTimesRhoS<0)*dPdxRFitT4M+(xTimesRhoS>0)*dPdxRFitT4P");
  // Show residulas
  treeUnit0->Draw("dPdxEuler001LRS-dPdxRFitT4:dPdxConstLRS>>his(200,-0.5,3,200,-0.2,0.2)","abs(xTimesRho)<0.2&&dPdxEuler001LRS!=0","colz");
}



/// fit dPdxCorrection formula and test the content
/// TODO - this is obsolete
void fitdPdxScaling(){
  TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  // make T4 fit in range -0.25,0.6
  TF1 f1("f1","x*(1+[0]*x+[1]*x*x+[2]*x*x*x)");
  f1.SetParName(0, "dX2");
  f1.SetParName(1, "dX3");
  f1.SetParName(2, "dX4");
  gStyle->SetOptFit(1);
  treeUnit0->Draw("((paramStepRK.P()-paramIn.P())/paramIn.P()):(AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*xTimesRho)/paramIn.P()>>his(300,-0.4,0.6,300,-0.8,0.4)",
                  "statusStepRK==1&&status0==1","colz");
  treeUnit0->GetHistogram()->Fit("f1","","",-0.25,.6);
  gPad->SaveAs("fig/fitdPdxScaling.png");
  // Check fit  as fitted
  treeUnit0->SetAlias("dPdxConst","(AliExternalTrackParam4D::dPdx(paramIn.P(),mass)*xTimesRho)/paramIn.P()");
  treeUnit0->SetAlias("dPdxR0","((paramStepRK.P()-paramIn.P())/paramIn.P())");
  treeUnit0->SetAlias("dPdxR","(AliExternalTrackParam4D::dPdx(paramIn.P(),mass+0)*xTimesRho)/paramIn.P()");
  treeUnit0->SetAlias("dPdxRFit",Form("dPdxR*(1+%f*dPdxR+%f*dPdxR**2+%f*dPdxR**3)",f1.GetParameter(0), f1.GetParameter(1),f1.GetParameter(2)));
  treeUnit0->Draw("dPdxR0-dPdxRFit:dPdxConst>>hisDiff(300,-0.4,0.6,300,-0.8,0.4)","statusStepRK==1","colz");
  gPad->SaveAs("fig/fitdPdxScalingDiffTest.png");
  // Check fit  as as put into code
  treeUnit0->Draw("((paramStepRK.P()-paramIn.P()-AliExternalTrackParam4D::dPdxCorr(paramIn.P(),mass,xTimesRho))/paramIn.P()):dPdxConst>>hisDiffCorr(300,-0.4,0.6,300,-0.8,0.6)",
                  "statusStepRK==1&&status0==1","colz");
  treeUnit0->GetHistogram()->Fit("pol1","w","",-0.25,.6);
  gPad->SaveAs("fig/fitdPdxScalingDiffCorr.png");
  // make test of the AliExternalTrackParam4D::dPdxCorr(paramIn.P(),mass,xTimesRho)/paramIn.P()
  treeUnit0->Draw("dPdxR0-AliExternalTrackParam4D::dPdxCorr(paramIn.P(),mass,xTimesRho)/paramIn.P():dPdxConst:mass","statusStepRK==1&&statusT4==1","colz");
  gPad->SaveAs("fig/fitdPdxCorrDiff.png");
}

void test_RungeKutaDraw(){
  TCanvas *canvasfitdPdxScaling = new TCanvas("fitdPdxScaling","fitdPdxScaling",1200,400);
  canvasfitdPdxScaling->SetGrid(1,1);
  // Euler propagation
  treeUnit0->Draw("log(paramStep.P()/param0.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&status0==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramEuler_PLoss.png");
  // Runge-Kuta in Energy
  treeUnit0->Draw("log(paramStepRK.P()/paramRK.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&statusRK==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramRK_PLoss.png");
  // Runge-Kuta in momenta
  treeUnit0->Draw("log(paramStepRK.P()/paramRKP.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&statusRKP==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
  gPad->SaveAs("fig/unitTest_paramRKP_PLoss.png");
  //
  treeUnit0->Draw("log(paramStepRK.P()/paramT4.P()):paramStepRK.P()/paramIn.P():mass","statusStepRK==1&&statusT4==1&&paramStepRK.P()/paramIn.P()>0.5","colz");
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