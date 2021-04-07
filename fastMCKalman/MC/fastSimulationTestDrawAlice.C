void drawCheckELoss05042020(){
  treeFast->Draw("log(part.fParamMC[Iteration$].fData.P()/part.fParamMC[Iteration$-1].fData.P()):dEdxExpSolidL1*sqrt(1+fP[2]**2+fP[3]**2)/part.fParamMC[Iteration$].fData.P():pidCode","part.fParamMC[].fData.P()<0.2&&Iteration$>1&&Iteration$<5&&Iteration$<fMaxLayer","colz",100000);
  // first version of drawing was wrong as fP[2] and fP[3] query was wrong
  // -->
  treeFast->Draw("log(part.fParamMC[Iteration$].fData.P()/part.fParamMC[Iteration$-1].fData.P()):dEdxExpSolidL1*sqrt(1+part.fParamMC[Iteration$].fData.fP[2]**2+part.fParamMC[Iteration$].fData.fP[3]**2)/part.fParamMC[Iteration$].fData.P():Entry$","part.fParamMC[].fData.P()<0.2&&Iteration$>1&&Iteration$<fMaxLayer&&pidCode==4","colz",100000);
  // this part is OK
}