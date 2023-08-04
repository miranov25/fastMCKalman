/*
  Here we would like to have tests - which checks internal consistence of the fastSimulation
  Dedicated debug streamaers to be used
      gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
.L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
     .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    .L $fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g
    initTreeFast();

    //test RDF
    rdf1=makeDataFrame(treeFast);
    rdf1.Range(1000).Snapshot("xxx","xxx.root",{"deltaIn0","deltaIn1"});
    rdf1.Range(100).Snapshot("dump","dumpAll.root");
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
#include "TStyle.h"
//
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"




extern TChain * treeSeed;
extern TChain * treeFast;
int paramsDFGlobal[5]={0,1,2,3,4};

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


void testPulls(std::string sv="", std::string Id="In", std::string extra_condition="&&Iteration$==0") {
  TF1 *mygauss = new TF1("mygauss", "gaus");
  int isOK=fastParticle::kTrackisOK;
    for (int iPar = 0; iPar <= 4; iPar++) {
      treeFast->Draw(Form("(part%s.fParam%s[].fP[%d]-part%s.fParamMC[].fP[%d])/sqrt(part%s.fParam%s[].fC[%d])>>his(100,-6,6)",sv.c_str(),Id.c_str(), iPar, sv.c_str(), iPar, sv.c_str(),Id.c_str(), AliExternalTrackParam::GetIndex(iPar, iPar)),
                    Form("(part%s.fStatusMask%s[]&%d)==%d&&abs(part%s.fParam%s[].fP[2])<0.7%s",sv.c_str(),Id.c_str(),isOK,isOK,sv.c_str(),Id.c_str(),extra_condition.c_str()), "");
      treeFast->GetHistogram()->Fit("mygauss", "q");
      bool isOK = abs(1 - mygauss->GetParameter(2)) < 5 * mygauss->GetParError(2);
      float rms=treeFast->GetHistogram()->GetRMS();
      if (isOK) {
        ::Info(Form("testFastTracker part%s reco %s pull test P%d ",sv.c_str(),Id.c_str(),iPar), "pullAnalytical - OK - %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      } else {
        ::Error(Form("testFastTracker part%s reco %s pull test P%d",sv.c_str(),Id.c_str(), iPar), "pullAnalytical- FAILED- %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      }
    }
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
  float sigmaPhi2=????;
  float sigmaTheta2=????;
  //
  for (int i=0; i<nSteps;i++ ){
     covarLinear[0]+=(theta2*dL*dL) ;    // + add contribution frpm angle phi from before
     covarLinear[2]+=(theta2*dL*dL);    // + add contribution frpm angle phi from before
     covarLinear[5]+=(theta2);          //
     covarLinear[9]+=(theta2);          //
     //
     covarLinear[3]+=theta2*dL;         // +
     covarLinear[7]+=theta2*dL;         // +
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

void drawTrackStatus(int counter, std::string Id = "In", std::string Error = "0x1"){
  treeFast->SetMarkerColor(1);   /// all MC
  treeFast->Draw("gyMC:gxMC","","",1,counter);
  treeFast->SetMarkerColor(4);   /// all reco points
  treeFast->Draw("gyInF:gxInF",Form("partFull.fStatusMask%s>0",Id.c_str()),"same",1,counter);
  treeFast->SetMarkerColor(2);   /// first MC point
  treeFast->Draw("gyMC:gxMC","Iteration$==0","same",1,counter);
  treeFast->SetMarkerColor(3);   /// trigger problem
  treeFast->Draw("gyInF:gxInF",Form("partFull.fStatusMask%s==%s",Id.c_str(),Error.c_str()),"same",1,counter);
}

void drawTrackStatus3D(int counter, std::string Id = "In", std::string Error = "0x1"){
  treeFast->SetMarkerColor(1);   /// all MC
  treeFast->Draw("gyMC:gxMC:gzMC","","",1,counter);
  treeFast->SetMarkerColor(4);   /// all reco points
  treeFast->Draw("gyInF:gxInF:gzInF",Form("partFull.fStatusMask%s>0",Id.c_str()),"same",1,counter);
  treeFast->SetMarkerColor(2);   /// first MC point
  treeFast->Draw("gyMC:gxMC:gzMC","Iteration$==0","same",1,counter);
  treeFast->SetMarkerColor(3);   /// trigger problem
  treeFast->Draw("gyInF:gxInF:gzInF",Form("partFull.fStatusMask%s==%s",Id.c_str(),Error.c_str()),"same",1,counter);
}

void SetList(std::string Id = "In", std::string Error = "0x1"){

  treeFast->SetAlias("gyInF",Form("sin(partFull.fParam%s[].fAlpha)*partFull.fParam%s[].fX",Id.c_str(),Id.c_str()));
  treeFast->SetAlias("gxInF",Form("cos(partFull.fParam%s[].fAlpha)*partFull.fParam%s[].fX",Id.c_str(),Id.c_str()));
  treeFast->SetAlias("gzInF",Form("partFull.fParam%s[].fP[1]",Id.c_str()));

  treeFast->Draw(">>ProblemRot",Form("Sum$(partFull.fStatusMask%s==%s)",Id.c_str(),Error.c_str()),"entrylist");
  TEntryList* problemList0x1 =(TEntryList*)gDirectory->Get("ProblemRot");
  treeFast->SetEntryList(problemList0x1);
  int counter=0;
  treeFast->SetMarkerSize(1.5);
  gStyle->SetPalette(55);
}

ROOT::RVec<float>  dEdxRVec1(ROOT::RVec<AliExternalTrackParam4D>& track, float &mass) {
  ROOT::RVec<float> dEdx(track.size());
  for (size_t i = 0; i < track.size(); i++) dEdx[i] = AliExternalTrackParam::BetheBlochSolid(track[i].GetP() / mass);
  return dEdx;
}

ROOT::RVec<float>  Beta(ROOT::RVec<AliExternalTrackParam4D>& track, float &mass) {
  ROOT::RVec<float> beta(track.size());
  for (size_t i = 0; i < track.size(); i++) beta[i] = sqrt((track[i].GetP()*track[i].GetP()) / (mass*mass+track[i].GetP()*track[i].GetP()));
  return beta;
}

ROOT::RVec<float>  AvgSqrtdEdxOverpTOut(ROOT::RVec<AliExternalTrackParam4D>& track, float &mass, int FirstIndex) {
  ROOT::RVec<float> dEdxOverpT(track.size());
  for (size_t i = 0; i < track.size(); i++) 
  {
    float dEdxOverpTTot=0;
    int nPoints=0;
    if(i>=size_t(FirstIndex)){    
      for(size_t j = 0; j <= size_t(i-FirstIndex); j++)
        {
          dEdxOverpTTot+=sqrt(AliExternalTrackParam::BetheBlochSolid(track[j].GetP() / mass)) * abs(track[j].GetParameter()[4]);
          nPoints++;
        }
      }
    dEdxOverpT[i] = dEdxOverpTTot;
    if(nPoints!=0) dEdxOverpT[i]/=nPoints;
  }
  return dEdxOverpT;
}

ROOT::RVec<float>  AvgSqrtdEdxOverpTIn(ROOT::RVec<AliExternalTrackParam4D>& track, float &mass, int FirstIndex) {
  ROOT::RVec<float> dEdxOverpT(track.size());
  for (size_t i = 0; i < track.size(); i++) 
  {
    float dEdxOverpTTot=0;
    int nPoints=0;
    if(i<=size_t(FirstIndex)){    
      for(size_t j = 0; j <= size_t(FirstIndex); j++)
        {
          dEdxOverpTTot+=sqrt(AliExternalTrackParam::BetheBlochSolid(track[j].GetP() / mass)) * abs(track[j].GetParameter()[4]);
          nPoints++;
        }
      }
    dEdxOverpT[i] = dEdxOverpTTot;
    if(nPoints!=0) dEdxOverpT[i]/=nPoints;
  }
  return dEdxOverpT;
}

ROOT::RVec<float>  AvgInvBetapTOut(ROOT::RVec<AliExternalTrackParam4D>& track, float &mass, int FirstIndex) {
  ROOT::RVec<float> InvbetapT(track.size());
  for (size_t i = 0; i < track.size(); i++) 
  {
    float InvbetapTTot=0;
    int nPoints=0;
    if(i>=size_t(FirstIndex)){    
      for(size_t j = 0; j <= size_t(i-FirstIndex); j++)
        {
          float beta = sqrt((track[i].GetP()*track[i].GetP()) / (mass*mass+track[i].GetP()*track[i].GetP()));
          InvbetapTTot+=abs(track[j].GetParameter()[4])* 1 / beta;
          nPoints++;
        }
      }
    InvbetapT[i] = InvbetapTTot;
    if(nPoints!=0) InvbetapT[i]/=nPoints;
  }
  return InvbetapT;
}

ROOT::RVec<float>  AvgInvBetapTIn(ROOT::RVec<AliExternalTrackParam4D>& track, float &mass, int FirstIndex) {
  ROOT::RVec<float> InvbetapT(track.size());
  for (size_t i = 0; i < track.size(); i++) 
  {
    float InvbetapTTot=0;
    int nPoints=0;
    if(i<=size_t(FirstIndex)){    
      for(size_t j = 0; j <= size_t(FirstIndex); j++)
        {
          float beta = sqrt((track[i].GetP()*track[i].GetP()) / (mass*mass+track[i].GetP()*track[i].GetP()));
          InvbetapTTot += abs(track[j].GetParameter()[4])* 1 / beta;
          nPoints++;
        }
      }
    InvbetapT[i] = InvbetapTTot;
    if(nPoints!=0) InvbetapT[i]/=nPoints;
  }
  return InvbetapT;
}

ROOT::RVec<float>  deltaP(ROOT::RVec<AliExternalTrackParam4D>& track, ROOT::RVec<AliExternalTrackParam4D>& trackMC, int param) {
  ROOT::RVec<float> deltaP(track.size());
  for (size_t i = 0; i < track.size(); i++) deltaP[i] = track[i].GetParameter()[param]-trackMC[i].GetParameter()[param];
  return deltaP;
}

ROOT::RVec<float>  pullP(ROOT::RVec<AliExternalTrackParam4D>& track, ROOT::RVec<AliExternalTrackParam4D>& trackMC, int param) {
  ROOT::RVec<float> pullP(track.size());
  for (size_t i = 0; i < track.size(); i++) pullP[i] = (track[i].GetParameter()[param]-trackMC[i].GetParameter()[param])/sqrt(track[i].GetCovariance()[AliExternalTrackParam::GetIndex(param, param)]);
  return pullP;
}


ROOT::RVec<float>  covarP(ROOT::RVec<AliExternalTrackParam4D>& track, int param) {
  ROOT::RVec<float> covarP(track.size());
  for (size_t i = 0; i < track.size(); i++) covarP[i] = TMath::Sqrt(track[i].GetCovariance()[AliExternalTrackParam::GetIndex(param, param)]);
  return covarP;
}

ROOT::RVec<float>  paramP(ROOT::RVec<AliExternalTrackParam4D>& track, int param) {
  ROOT::RVec<float> paramP(track.size());
  for (size_t i = 0; i < track.size(); i++) paramP[i] = track[i].GetParameter()[param];
  return paramP;
}

ROOT::RVec<float>  LArm(ROOT::RVec<AliExternalTrackParam4D>& track, int FirstIndex) {
  ROOT::RVec<float> LArm(track.size());
  if(FirstIndex<0) FirstIndex=0;
  for (size_t i = 0; i < track.size(); i++) {
    Double_t xyzMC0[3];
    Double_t xyzi[3];
    track[FirstIndex].GetXYZ(xyzMC0);
    track[i].GetXYZ(xyzi);
    LArm[i] = TMath::Sqrt((xyzMC0[0]-xyzi[0])*(xyzMC0[0]-xyzi[0])+(xyzMC0[1]-xyzi[1])*(xyzMC0[1]-xyzi[1]));
  }
  return LArm;
}

ROOT::RVec<float>  LengthOut(ROOT::RVec<AliExternalTrackParam4D>& track, int FirstIndex) {
  ROOT::RVec<float> Length(track.size());
  for (size_t i = 0; i < track.size(); i++) {
    float Lengthtot = 0;
    if(i>size_t(FirstIndex)){    
      for(size_t j = 0; j < size_t(i-FirstIndex); j++)
        {
          Double_t xyz0[3];
          Double_t xyzi[3];
          track[FirstIndex+j].GetXYZ(xyz0);
          track[FirstIndex+j+1].GetXYZ(xyzi);
          Lengthtot+=TMath::Sqrt((xyz0[0]-xyzi[0])*(xyz0[0]-xyzi[0])+(xyz0[1]-xyzi[1])*(xyz0[1]-xyzi[1]));
        }
      }
    Length[i] = Lengthtot;
  }
  return Length;
}

ROOT::RVec<float>  LengthIn(ROOT::RVec<AliExternalTrackParam4D>& track, int FirstIndex) {
  ROOT::RVec<float> Length(track.size());
  for (size_t i = 0; i < track.size(); i++) {
    float Lengthtot = 0;
    if(i<size_t(FirstIndex)){ 
      for(size_t j = i; j < size_t(FirstIndex); j++)
        {
          Double_t xyzMC0[3];
          Double_t xyzi[3];
          track[j].GetXYZ(xyzMC0);
          track[j+1].GetXYZ(xyzi);
          Lengthtot+=TMath::Sqrt((xyzMC0[0]-xyzi[0])*(xyzMC0[0]-xyzi[0])+(xyzMC0[1]-xyzi[1])*(xyzMC0[1]-xyzi[1]));
        }
    }
    Length[i] = Lengthtot;
  }
  return Length;
}


ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>  makeDataFrame(TTree * treeFast){
  // problem defining parameters in scope
  std::vector<std::string> varList;
  std::string paramType[4]={"In","Out","Refit", "MC"};
  //int params[5]={0,1,2,3,4};
  //
  auto rdf=ROOT::RDataFrame(*treeFast);
  auto rdf1=rdf.Define("sigmaRPhi","geom.fLayerResolRPhi[0]");
  rdf1=rdf1.Define("sigmaZ","geom.fLayerResolZ[0]");
  rdf1=rdf1.Define("layerX0","1/geom.fLayerX0[0]");
  rdf1=rdf1.Define("massMC","partFull.fMassMC");
  rdf1=rdf1.Define("ptMC0",[](ROOT::RVec<AliExternalTrackParam4D> &track){return track[0].Pt();},{"partFull.fParamMC"});
  rdf1=rdf1.Define("ptMCEnd",[](ROOT::RVec<AliExternalTrackParam4D> &track){return track[track.size()-1].Pt();},{"partFull.fParamMC"});

  rdf1=rdf1.Define("dEdxMC",dEdxRVec1,{"partFull.fParamMC","partFull.fMassMC"});
  rdf1=rdf1.Define("dEdxIn",dEdxRVec1,{"partFull.fParamIn","partFull.fMassMC"});
  rdf1=rdf1.Define("dEdxOut",dEdxRVec1,{"partFull.fParamOut","partFull.fMassMC"});
  rdf1=rdf1.Define("dEdxRefit",dEdxRVec1,{"partFull.fParamRefit","partFull.fMassMC"});
  rdf1=rdf1.Define("dEdxInLeg",dEdxRVec1,{"part.fParamIn","part.fMassMC"});

  rdf1=rdf1.Define("BetaMC",Beta,{"partFull.fParamMC","partFull.fMassMC"});
  rdf1=rdf1.Define("BetaIn",Beta,{"partFull.fParamIn","partFull.fMassMC"});
  rdf1=rdf1.Define("BetaOut",Beta,{"partFull.fParamOut","partFull.fMassMC"});
  rdf1=rdf1.Define("BetaRefit",Beta,{"partFull.fParamRefit","partFull.fMassMC"});
  rdf1=rdf1.Define("BetaInLeg",Beta,{"part.fParamIn","part.fMassMC"});
  //
  //
  for (int i=0; i<5; i++){
    for (int iType=0; iType<4; iType++) {
      const char *type=paramType[iType].data();
      //
      rdf1 = rdf1.Define(Form("param%s%d", type, i),
                         [i](ROOT::RVec <AliExternalTrackParam4D> &track) {
                           return paramP(track, paramsDFGlobal[i]);
                         },
                         {Form("partFull.fParam%s",type)});
      if(iType==0){
        rdf1 = rdf1.Define(Form("param%sLeg%d", type, i),
                         [i](ROOT::RVec <AliExternalTrackParam4D> &track) {
                           return paramP(track, paramsDFGlobal[i]);
                         },
                         {Form("part.fParam%s",type)});
      }
      if (iType==3) continue;  /// skip delta and covar for the MC part
      rdf1 = rdf1.Define(Form("delta%s%d",type, i),
                         [i](ROOT::RVec <AliExternalTrackParam4D> &track, ROOT::RVec <AliExternalTrackParam4D> &trackMC) {
                           return deltaP(track, trackMC, paramsDFGlobal[i]);
                         },
                         {Form("partFull.fParam%s",type), "partFull.fParamMC"});
      rdf1 = rdf1.Define(Form("covar%s%d",type, i),
                         [i](ROOT::RVec <AliExternalTrackParam4D> &track) {
                           return covarP(track, paramsDFGlobal[i]);
                         },
                         {Form("partFull.fParam%s",type)});
      rdf1 = rdf1.Define(Form("pull%s%d",type, i),
                         [i](ROOT::RVec <AliExternalTrackParam4D> &track, ROOT::RVec <AliExternalTrackParam4D> &trackMC) {
                           return pullP(track, trackMC, paramsDFGlobal[i]);
                         },
                         {Form("partFull.fParam%s",type), "partFull.fParamMC"});

      if(iType==0){
        rdf1 = rdf1.Define(Form("delta%sLeg%d",type, i),
                         [i](ROOT::RVec <AliExternalTrackParam4D> &track, ROOT::RVec <AliExternalTrackParam4D> &trackMC) {
                           return deltaP(track, trackMC, paramsDFGlobal[i]);
                         },
                         {Form("part.fParam%s",type), "part.fParamMC"});
        rdf1 = rdf1.Define(Form("covar%sLeg%d",type, i),
                          [i](ROOT::RVec <AliExternalTrackParam4D> &track) {
                            return covarP(track, paramsDFGlobal[i]);
                          },
                          {Form("part.fParam%s",type)});
        rdf1 = rdf1.Define(Form("pull%sLeg%d",type, i),
                          [i](ROOT::RVec <AliExternalTrackParam4D> &track, ROOT::RVec <AliExternalTrackParam4D> &trackMC) {
                            return pullP(track, trackMC, paramsDFGlobal[i]);
                          },
                          {Form("part.fParam%s",type), "part.fParamMC"});
      }
    }
  }

   for (int iType=0; iType<4; iType++) {
     const char *type=paramType[iType].data();
     if(iType==2) continue;
     rdf1 = rdf1.Define(Form("LArm%s",type),LArm,{Form("partFull.fParam%s",type), Form("partFull.fFirstIndex%s",type)});
     if(iType==0)
     {
      rdf1 = rdf1.Define(Form("LengthXY%s",type),LengthIn,{Form("partFull.fParam%s",type), Form("partFull.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("AvgInvBetapT%s",type),AvgInvBetapTIn,{Form("partFull.fParam%s",type), "partFull.fMassMC", Form("partFull.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("AvgSqrtdEdxOverpT%s",type),AvgSqrtdEdxOverpTIn,{Form("partFull.fParam%s",type), "partFull.fMassMC", Form("partFull.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("LengthXY%sLeg",type),LengthIn,{Form("part.fParam%s",type), Form("part.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("LArm%sLeg",type),LArm,{Form("part.fParam%s",type), Form("part.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("AvgInvBetapT%sLeg",type),AvgInvBetapTIn,{Form("part.fParam%s",type), "part.fMassMC", Form("part.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("AvgSqrtdEdxOverpT%sLeg",type),AvgSqrtdEdxOverpTIn,{Form("part.fParam%s",type), "part.fMassMC", Form("part.fFirstIndex%s",type)});
     }
     else{
      rdf1 = rdf1.Define(Form("LengthXY%s",type),LengthOut,{Form("partFull.fParam%s",type), Form("partFull.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("AvgSqrtdEdxOverpT%s",type),AvgSqrtdEdxOverpTOut,{Form("partFull.fParam%s",type), "partFull.fMassMC", Form("partFull.fFirstIndex%s",type)});
      rdf1 = rdf1.Define(Form("AvgInvBetapT%s",type),AvgInvBetapTOut,{Form("partFull.fParam%s",type), "partFull.fMassMC", Form("partFull.fFirstIndex%s",type)});
     }
  }

  rdf1 = rdf1.Define("AvgInvBetapTRefit","(AvgInvBetapTIn+AvgInvBetapTOut)/2");
  
  rdf1 = rdf1.Define("LArmRefit","LArmIn+LArmOut");
  rdf1 = rdf1.Define("LengthXYRefit","LengthXYIn+LengthXYOut");
  rdf1 = rdf1.Define("AvgSqrtdEdxOverpTRefit","(AvgSqrtdEdxOverpTIn+AvgSqrtdEdxOverpTOut)/2");

  rdf1 = rdf1.Define("covarIn4ExpLPTMoliere","10.674*abs(paramIn4)/(BetaIn*sqrt(LengthXYIn*layerX0))");
  rdf1 = rdf1.Define("covarInLeg4ExpLPTMoliere","10.674*abs(paramInLeg4)/(BetaInLeg*sqrt(LengthXYInLeg*layerX0))");
  rdf1 = rdf1.Define("covarOut4ExpLPTMoliere","10.674*abs(paramOut4)/(BetaOut*sqrt(LengthXYOut*layerX0))");
  rdf1 = rdf1.Define("covarRefit4ExpLPTMoliere","10.674*abs(paramRefit4)/(BetaRefit*sqrt(LengthXYRefit*layerX0))");

  rdf1 = rdf1.Define("covarIn4ExpLPTMoliereAvg","10.674*AvgInvBetapTIn/(sqrt(LengthXYIn*layerX0))");
  rdf1 = rdf1.Define("covarInLeg4ExpLPTMoliereAvg","10.674*AvgInvBetapTInLeg/(sqrt(LengthXYInLeg*layerX0))");
  rdf1 = rdf1.Define("covarOut4ExpLPTMoliereAvg","10.674*AvgInvBetapTOut/(sqrt(LengthXYOut*layerX0))");
  rdf1 = rdf1.Define("covarRefit4ExpLPTMoliereAvg","10.674*AvgInvBetapTRefit/(sqrt(LengthXYRefit*layerX0))");

  rdf1 = rdf1.Define("covarIn4ExpHPT","17900.928*sigmaRPhi/((LArmIn*LArmIn)*(sqrt(partFull.fNPointsIn)))");
  rdf1 = rdf1.Define("covarInLeg4ExpHPT","17900.928*sigmaRPhi/((LArmInLeg*LArmInLeg)*(sqrt(part.fNPointsIn)))");
  rdf1 = rdf1.Define("covarOut4ExpHPT","17900.928*sigmaRPhi/((LArmOut*LArmOut)*(sqrt(partFull.fNPointsOut)))");
  rdf1 = rdf1.Define("covarRefit4ExpHPT","17900.928*sigmaRPhi/((LArmRefit*LArmRefit)*(sqrt(partFull.fNPointsRefit)))");
  return rdf1;
}


std::string testPullsSnapshot(std::string name="testVarRDF",std::string Id="In", std::string extra_condition="") {
  std::string filename = name+".root";
  TChain * tree = new TChain(name.c_str());
  tree->Add(filename.c_str());
  TF1 *mygauss = new TF1("mygauss", "gaus");
  int isOK=fastParticle::kTrackisOK;
  std::string results="";
    for (int iPar = 0; iPar <= 4; iPar++) {
      tree->Draw(Form("pull%s%d[]>>his(100,-6,6)",Id.c_str(), iPar),
                    Form("(StatusMask%s[]&%d)==%d&&abs(param%s2[])<0.7%s",Id.c_str(),isOK,isOK,Id.c_str(),extra_condition.c_str()), "");
      tree->GetHistogram()->Fit("mygauss", "q");
      bool pullisOK = abs(1 - mygauss->GetParameter(2)) < 5 * mygauss->GetParError(2);
      float rms=tree->GetHistogram()->GetRMS();
      if (pullisOK) {
        results += Form("testFastTracker snapShot reco %s pull test P%d - pullAnalytical - OK - %2.2f\t%2.2f \n",Id.c_str(),iPar, mygauss->GetParameter(2),rms);
      } else {
        results += Form("testFastTracker snapShot reco %s pull test P%d - pullAnalytical - FAILED - %2.2f\t%2.2f \n",Id.c_str(),iPar, mygauss->GetParameter(2),rms);
      }
    }
  return results;
}

void drawdEdxResolution(){
   //treeFast->Draw("partFull.@fParamMC.size():partFull.fParamMC[0].P():densScaling","hasDecay==0&&pidCode==1&&partFull.fStatusStopMC==0x20&&partFull.fParamMC[0].P()<0.1","colz")

  treeFast->SetAlias("sigmadEdxRel","0.06*(1/(densScaling**0.25))*sqrt(120./partFull.@fParamMC.size())*(1./(dEdxExp**0.25))");
  // sigmadEdx=k*((2/Beta^3)*sigmaBeta= 2*dedx/Beta*sigmaBeta --> sigmaBeta=Beta*(sigmadEdx/dEdx)/2.
  treeFast->SetAlias("betaResRel","sigmadEdxRel/2."); // 1/Beta^2 ~ dEdx  sigmadEdx=k*((2/Beta^3)*sigmaBeta= dedx/Beta*sigmaBeta
  treeFast->SetAlias("pResRel","betaResRel*2");    // in the region where Beta<<1
  treeFast->Draw("pResRel:1/ptMC:densScaling","ptMC<2&&ptMC>0.05&&pidCode==4&&hasDecay==0&&densScaling>5","colz",100000);


}