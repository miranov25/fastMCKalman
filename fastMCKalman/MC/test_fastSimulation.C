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
  for (size_t i = 0; i < track.size(); i++) {
    Double_t xyzMC0[3];
    Double_t xyzi[3];
    track[FirstIndex].GetXYZ(xyzMC0);
    track[i].GetXYZ(xyzi);
    LArm[i] = TMath::Sqrt((xyzMC0[0]-xyzi[0])*(xyzMC0[0]-xyzi[0])+(xyzMC0[1]-xyzi[1])*(xyzMC0[1]-xyzi[1]));
  }
  return LArm;
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
  rdf1=rdf1.Define("layerX0","geom.fLayerX0[0]");
  rdf1=rdf1.Define("dEdxMC",dEdxRVec1,{"partFull.fParamMC","partFull.fMassMC"});
  rdf1=rdf1.Define("dEdxIn",dEdxRVec1,{"partFull.fParamIn","partFull.fMassMC"});
  //
  varList.push_back("sigmaRPhi"); varList.push_back("sigmaZ");  varList.push_back("dEdxMC");
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
      varList.push_back(Form("param%s%d", type, i));
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
      varList.push_back(Form("delta%s%d", type, i));
      varList.push_back(Form("covar%s%d", type, i));
      varList.push_back(Form("pull%s%d", type, i));
    }
  }

   for (int iType=0; iType<4; iType++) {
     const char *type=paramType[iType].data();
     if(iType==2) continue;
     rdf1 = rdf1.Define(Form("LArm%s",type),LArm,{Form("partFull.fParam%s",type), Form("partFull.fFirstIndex%s",type)});
     varList.push_back(Form("LArm%s", type));
  }
  
  rdf1 = rdf1.Define("LArmRefit","LArmIn+LArmOut");
  varList.push_back(Form("LArmRefit"));
  //rdf1.Snapshot("xxx","xxx.root",varList);
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