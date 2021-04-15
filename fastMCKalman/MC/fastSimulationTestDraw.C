void drawRungeKuttaTest(){
  treeUnit0->Draw("log(paramStep.P()/paramRK2.P()):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
  gPad->SaveAs("unitTest_paramRK_PLoss.png");
  treeUnit0->Draw("log(paramStep.P()/param0.P()):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
  gPad->SaveAs("unitTest_param0_PLoss.png");
}

void drawRungeKuttaTest() {

  treeUnit0->Draw("(paramStep.fC[9]-param0.fC[9])/(paramStep.fC[9]-paramIn.fC[9]):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
  treeUnit0->Draw("(paramStep.fC[9]-paramRK.fC[9])/(paramStep.fC[9]-paramIn.fC[9]):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");
  treeUnit0->Draw("(paramStep.fC[9]-paramRK2.fC[9])/(paramStep.fC[9]-paramIn.fC[9]):paramRK2.P()/paramIn.P():mass","xOverX0>0.01&&statusStep==1&&paramRK2.P()/paramIn.P()>0.6","colz");


}