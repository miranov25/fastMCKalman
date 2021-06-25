/*
  .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C
   .L $fastMCKalman/fastMCKalman/MC/fastSimulationFit.C
  initTreeFast();
   treeFast->SetAlias("meanptMC","Sum$(part.fParamMC.fData.Pt())/Sum$(part.fParamMC.fData.Pt()>0)");
   treeFast->SetAlias("sigmaPtRel0","sigmaqPt0*ptMC");

   makeFitResol(treeFast, "sigmaqPt0", "(10/sigmaqPt0)**2", "", "ptMC>0.01&&fMaxLayer>5&&ptMC<10",kFALSE);
    makeFitResolLArm(treeFast, "sigmaqPt0", "(10/sigmaqPt0)**2", "", "ptMC>0.01&&fMaxLayer>5&&ptMC<10",kFALSE);
     makeFitResolLArm(treeFast, "sigmaPtRel0", "(10/sigmaPtRel0)**2", "", "ptMC>0.01&&fMaxLayer>5&&ptMC<10",kFALSE);

*/
TF1 *likePseudoHuber = new TF1("likePseudoHuber", AliTMinuitToolkit::PseudoHuberLogLike,-10,10,1);
TStopwatch timer;
TStopwatch timerFunction;
AliTMinuitToolkit *fitter=0;

/// make a fit and set alas value for fit
///      fitter is in global scope parameters could be acessed via fitter
/// \param varFit            - variable to fit  sigmaqPt0
/// \param varWeight         - 1/errror
/// \param fitSelection      -
/// \param doCheck           -
/// Example usage:
void makeFitResol(TTree * treeFast, TString varFit, TString varWeight, TString varPrefix, TString fitSelection, bool doCheck=true){
  //  Int_t nPoints=1000; TString varFit="sigmaqPt0";  TString varPrefix=""; TString fitSelection="ptMC>0.01&&fMaxLayer>5"; TString varWeight="10/sigmaqPt0"
  //  makeFitResol(treeFast, "sigmaqPt0", "10/sigmaqPt0", "", "ptMC>0.01&&fMaxLayer>5",kFALSE)
  //
  //
  // variables:
  //     x[0]=q/pt
  //     x[1]=MIP/dEdx
  //     x[2]=tgl
  //     x[3]=alpha
  // parameters:
  //    [0] - resol 0
  //    [1] - qPt factor
  //    [2] - qPt power
  //    [3] - dEdx power
  ::Info("makeFitBG","BEGIN");
  timer.Start();
  timerFunction.Start();
  likePseudoHuber->SetParameter(0,3);
  const Int_t nPoints=1000000;
  TF1 *formula = new TF1("dcaResol0", "sqrt([0]**2+[1]*(abs(x[0])**[2])*(x[1]**[3]))");
  TMatrixD *initParam = new TMatrixD(4, 4), &p = *initParam;
  fitter = new AliTMinuitToolkit("fitterdEdxBG.root");
  fitter->SetFitFunction((TF1 *) formula, kTRUE);
  fitter->SetVerbose(0x7);
  // initial parameters
  // row    - corresponds to parameter
  // column - 0 - starting value, 1 - sigma, 2- lower bound , 3- upper bound
  p(0, 0)=0.1; p(0, 1)=0.1; p(0, 2)=0.;  p(0, 3)=1.;
  p(1, 0)=1;   p(1, 1)=0.1; p(1, 2)=0.;  p(1, 3)=1.;
  p(2, 0)=2;   p(2, 1)=0.1; p(2, 2)=-2.; p(2, 3)=4.;
  p(3, 0)=-1;  p(3, 1)=0.1; p(3, 2)=-4.; p(3, 3)=4.;
  //
  fitter->SetInitialParam(initParam);
  fitter->SetLogLikelihoodFunction(likePseudoHuber);
  //
  // const char *var = "LHC18r_p3.dcaRNtpcSignalPrim3Dist.rmsG";
  fitter->FillFitter(treeFast, (varFit+":"+varWeight).Data(), "1/meanptMC:1/dEdxExp", fitSelection.Data(), 0, nPoints, kTRUE);
  if (fitter->GetValues()->GetNrows()<=0){
    ::Error("","Invalid selection");
    return ;
  }
  fitter->SetTolerance(0.0001);
  fitter->SetMaxCalls(4000);
  timer.Start();
  fitter->Fit();
  timer.Print();
  treeFast->SetAlias(varFit+"Fit", fitter->GetFitFunctionAsAlias());
  //

  if (doCheck) {
    timer.Start();
    fitter->Bootstrap(10, "xxx");
    treeFast->SetAlias(varFit+"FitB10", fitter->GetFitFunctionAsAlias());
    fitter->GetRMSEstimator()->Print();
    //
    fitter->MISAC(100, 100, "");
    treeFast->SetAlias(varFit+"FitMI100", fitter->GetFitFunctionAsAlias());
    fitter->GetMISACEstimators()->Print();
    //
    //cons char *fdata
    Int_t entries= treeFast->Draw(Form("log(%sFitMI100/%sFit):log(%sFitB10/%sFit):log(%sFit/%s)",varFit.Data(), varFit.Data(),varFit.Data(),varFit.Data(),varFit.Data(),varFit.Data()),
                            fitSelection,"goff");
    if (entries>0) {
      ::Info("makeResolFitHighMult", "CHECK FITMI100/Fit \tMean:\t%0.4f\tRMS\t%0.4F", TMath::Median(entries, treeFast->GetV1()),
             TMath::RMS(entries, treeFast->GetV1()));
      ::Info("makeResolFitHighMult", "CHECK FITB19/Fit   \tMean:\t%0.4f\tRMS\t%0.4F", TMath::Median(entries, treeFast->GetV2()),
             TMath::RMS(entries, treeFast->GetV2()));
      ::Info("makeResolFitHighMult", "CHECK Fit/FitLocal \tMean%0.4f\tRMS\t%0.4F", TMath::Median(entries, treeFast->GetV3()),
             TMath::RMS(entries, treeFast->GetV3()));
       treeFast->GetListOfAliases()->Print("",(varFit+"*").Data());
       ::Info("makeResolFitHighMult","%s",fitter->GetFitFunctionAsAlias("latex").Data());
    }
  }
  ::Info("makeFitBG","END");
  timerFunction.Print();
}



/// make a fit and set alas value for fit
///      fitter is in global scope parameters could be acessed via fitter
/// \param varFit            - variable to fit  sigmaqPt0
/// \param varWeight         - 1/errror
/// \param fitSelection      -
/// \param doCheck           -
/// Example usage:
void makeFitResolLArm(TTree * treeFast, TString varFit, TString varWeight, TString varPrefix, TString fitSelection, bool doCheck=true){
  //  Int_t nPoints=1000; TString varFit="sigmaqPt0";  TString varPrefix=""; TString fitSelection="ptMC>0.01&&fMaxLayer>5"; TString varWeight="10/sigmaqPt0"
  //  makeFitResolLArm(treeFast, "sigmaqPt0", "(10/sigmaqPt0)**2", "", "ptMC>0.4&&fMaxLayer>5",kFALSE)
  //
  //
  // variables:
  //     x[0]=q/pt
  //     x[1]=MIP/dEdx
  //     x[2]=tgl
  //     x[3]=LArm
  // parameters:
  //    [0] - resol 0
  //    [1] - qPt factor
  //    [2] - qPt power
  //    [3] - dEdx power
  //    [4] - LArm power for resolution - should be -2
  //    [5] - LArm power for resolution - should be -1
  ::Info("makeFitBG","BEGIN");
  timer.Start();
  timerFunction.Start();
  likePseudoHuber->SetParameter(0,3);
  const Int_t nPoints=1000000;
  TF1 *formula = new TF1("dcaResol0", "sqrt(([0]**2)*(x[3]**[4])+[1]*(abs(x[0])**[2])*(x[1]**[3])*(x[3]**[5]))");
  TMatrixD *initParam = new TMatrixD(6, 4), &p = *initParam;
  fitter = new AliTMinuitToolkit("fitterdEdxBG.root");
  fitter->SetFitFunction((TF1 *) formula, kTRUE);
  fitter->SetVerbose(0x7);
  // initial parameters
  // row    - corresponds to parameter
  // column - 0 - starting value, 1 - sigma, 2- lower bound , 3- upper bound
  p(0, 0)=0.1; p(0, 1)=0.1; p(0, 2)=0.;  p(0, 3)=1.;
  p(1, 0)=1;   p(1, 1)=0.1; p(1, 2)=0.;  p(1, 3)=1.;
  p(2, 0)=2;   p(2, 1)=0.1; p(2, 2)=-2.; p(2, 3)=4.;
  p(3, 0)=-1;  p(3, 1)=0.1; p(3, 2)=-4.; p(3, 3)=4.;
  p(4, 0)=-4;  p(4, 1)=0.1; p(4, 2)=-5.; p(4, 3)=5.;
  p(5, 0)=-2;  p(5, 1)=0.1; p(5, 2)=-5.; p(5, 3)=5.;

  //
  fitter->SetInitialParam(initParam);
  fitter->SetLogLikelihoodFunction(likePseudoHuber);
  //
  // const char *var = "LHC18r_p3.dcaRNtpcSignalPrim3Dist.rmsG";
  fitter->FillFitter(treeFast, (varFit+":"+varWeight).Data(), "1/meanptMC:1/dEdxExp:tglMC:Larm", fitSelection.Data(), 0, nPoints, kTRUE);
  if (fitter->GetValues()->GetNrows()<=0){
    ::Error("","Invalid selection");
    return ;
  }
  fitter->SetTolerance(0.0001);
  fitter->SetMaxCalls(2000);
  timer.Start();
  fitter->Fit();
  timer.Print();
  treeFast->SetAlias(varFit+"FitL", fitter->GetFitFunctionAsAlias());
  //

  if (doCheck) {
    timer.Start();
    fitter->Bootstrap(10, "xxx");
    treeFast->SetAlias(varFit+"FitB10", fitter->GetFitFunctionAsAlias());
    fitter->GetRMSEstimator()->Print();
    //
    fitter->MISAC(100, 100, "");
    treeFast->SetAlias(varFit+"FitMI100", fitter->GetFitFunctionAsAlias());
    fitter->GetMISACEstimators()->Print();
    //
    //cons char *fdata
    Int_t entries= treeFast->Draw(Form("log(%sFitMI100/%sFit):log(%sFitB10/%sFit):log(%sFit/%s)",varFit.Data(), varFit.Data(),varFit.Data(),varFit.Data(),varFit.Data(),varFit.Data()),
                            fitSelection,"goff");
    if (entries>0) {
      ::Info("makeResolFitHighMult", "CHECK FITMI100/Fit \tMean:\t%0.4f\tRMS\t%0.4F", TMath::Median(entries, treeFast->GetV1()),
             TMath::RMS(entries, treeFast->GetV1()));
      ::Info("makeResolFitHighMult", "CHECK FITB19/Fit   \tMean:\t%0.4f\tRMS\t%0.4F", TMath::Median(entries, treeFast->GetV2()),
             TMath::RMS(entries, treeFast->GetV2()));
      ::Info("makeResolFitHighMult", "CHECK Fit/FitLocal \tMean%0.4f\tRMS\t%0.4F", TMath::Median(entries, treeFast->GetV3()),
             TMath::RMS(entries, treeFast->GetV3()));
       treeFast->GetListOfAliases()->Print("",(varFit+"*").Data());
       ::Info("makeResolFitHighMult","%s",fitter->GetFitFunctionAsAlias("latex").Data());
    }
  }
  ::Info("makeFitBG","END");
  timerFunction.Print();
}

/// TPC  - chenging lever arm using different points ... for non homogenous detector not working
void makeFitExample(){
   makeFitResolLArm(treeFast, "sigmaqPt0", "(10/sigmaqPt0)**2", "", "ptMC>0.01&&fMaxLayer>5&&ptMC<10",kFALSE);
   //
   treeFast->SetAlias("sigmaqPt","sqrt(part.fParamIn[].fC[14]+0)");
  treeFast->SetAlias("Larm","Max$(part.fParamMC[].fX)-Min$(part.fParamMC[Iteration$].fX)");
  makeFitResolLArm(treeFast, "sigmaqPt", "(10/sigmaqPt)**2", "", "ptMC>0.01&&sigmaqPt>0&&fMaxLayer>6&&ptMC<10&&rndm<0.02",kFALSE);
  //

}


