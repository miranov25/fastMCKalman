/*
   gSystem->AddIncludePath("-I$NOTES/calibration/TPCSplineCreationFramework_Skimmed/code_treenSigmaParam/")
  // library to be created with makefile before -  with the same aliroot/root/container version
  gSystem->Load("$NOTES/calibration/TPCSplineCreationFramework_Skimmed/code_treenSigmaParam/AliSkimmedDataAnalysisMaker.so")
  AliSkimmedDataAnalysisMaker::doTest()

 */

#include <iostream>
#include <vector>
#include "V0FlatAna.h"
//#include "AliESDtrack.h"
#include "TrackFlatAna.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TVectorD.h"
#include "TSpline.h"
//#include "THnSparseF.h"
#include "AliSkimmedDataAnalysisMaker.h"

using namespace std;

ClassImp(AliSkimmedDataAnalysisMaker);

///
void AliSkimmedDataAnalysisMaker::doTest(){
  AliSkimmedDataAnalysisMaker m;
  m.bookHistogram();
  m.TreePrimary_BBFitAnalysis("CleanTrack.root");
  m.TreeV0_BBFitAnalysis("V0tree.root");
  m.WriteHistogram();

}


AliSkimmedDataAnalysisMaker::AliSkimmedDataAnalysisMaker() {

  mOutputFile_BBFitAnalysis = new TFile("TPCTree_BBFitAnalysis.root","RECREATE"); 
}
AliSkimmedDataAnalysisMaker::~AliSkimmedDataAnalysisMaker(){
  //..
  //..
}
void AliSkimmedDataAnalysisMaker::bookHistogram()
{
  //tree for BBfit analysis

for (Int_t i = 0; i<7; i++){
  tpcNsig[i]=-999;
  tofNsig[i]=-999;
}
  fTree2 = new TTree("fTree2", "Tree for analysis of BetheBloch fit");
  fTree2->Branch("maker",this);
  return;
//  fTree2->Branch("p", &p);
//  fTree2->Branch("oneoverpt",&oneoverpt);
//  fTree2->Branch("rawtpcsignal",&rawtpcsignal);
//  fTree2->Branch("rawtpcsignal_leg",&rawtpcsignal_leg);
//  fTree2->Branch("itssignal",&itssignal);
//  fTree2->Branch("pileupcor_tpcsignal",&pileupcor_tpcsignal);
//  fTree2->Branch("pileupcor_tpcsignal_leg",&pileupcor_tpcsignal_leg);
//  fTree2->Branch("BG", &BG);
//  fTree2->Branch("eta",&eta);
//  fTree2->Branch("tgl",&tgl);
//  fTree2->Branch("multiplicty",&multiplicty);
//  fTree2->Branch("centrality", &centrality);
//  fTree2->Branch("PDGcode", &PDGcode);
//  fTree2->Branch("isV0",&isV0);
//  fTree2->Branch("V0Pull",&V0Pull);
//  fTree2->Branch("V0PullEff",&V0PullEff);
//  fTree2->Branch("K0PullEff",&K0PullEff);
//  fTree2->Branch("LPullEff",&LPullEff);
//  fTree2->Branch("ALPullEff",&ALPullEff);
//  fTree2->Branch("EPullEff",&EPullEff);
//  fTree2->Branch("isPileUp",&isPileUp);
//  //fTree2->Branch("dEdxExpected_SatLund_woDeut",&dEdxExpected_SatLund_woDeut);
//  //fTree2->Branch("dEdxExpected_SatLund_wDeut",&dEdxExpected_SatLund_wDeut);
//  //fTree2->Branch("dEdxExpected_ALEPH_woDeut",&dEdxExpected_ALEPH_woDeut);
//  //fTree2->Branch("dEdxExpected_ALEPH_wDeut",&dEdxExpected_ALEPH_wDeut);
//  //fTree2->Branch("dEdxExpected_SatLund_woDeut_pileUpCut",&dEdxExpected_SatLund_woDeut_pileUpCut);
//  //fTree2->Branch("tpcnsigma",&tpcnsigma);
//  //fTree2->Branch("tofnsigma",&tofnsigma);
//  fTree2->Branch("tpcNsig",&tpcNsig);
//  fTree2->Branch("tofNsig",&tofNsig);
//  fTree2->Branch("shiftM",&shiftM);
//  fTree2->Branch("nPileUpPrim", &nPileUpPrim);
//  fTree2->Branch("nPileUpSum", &nPileUpSum);
//  fTree2->Branch("primMult",&primMult);
//  fTree2->Branch("tpcClusterMult", &tpcClusterMult);
//  fTree2->Branch("pileUpOK", &pileUpOK);
//  fTree2->Branch("multSSD", &multSSD);
//  fTree2->Branch("multSDD", &multSDD);
//  fTree2->Branch("SignalTot0", &SignalTot0);
//  fTree2->Branch("SignalTot1", &SignalTot1);
//  fTree2->Branch("SignalTot2", &SignalTot2);
//  fTree2->Branch("SignalTot3", &SignalTot3);
//  fTree2->Branch("SignalMax0", &SignalMax0);
//  fTree2->Branch("SignalMax1", &SignalMax1);
//  fTree2->Branch("SignalMax2", &SignalMax2);
//  fTree2->Branch("SignalMax3", &SignalMax3);
//  fTree2->Branch("SignalNcl0", &SignalNcl0);
//  fTree2->Branch("SignalNcl1", &SignalNcl1);
//  fTree2->Branch("SignalNcl2", &SignalNcl2);
//  fTree2->Branch("SignalNcl3", &SignalNcl3);
//  fTree2->Branch("SignalNcr0", &SignalNcr0);
//  fTree2->Branch("SignalNcr1", &SignalNcr1);
//  fTree2->Branch("SignalNcr2", &SignalNcr2);
//  fTree2->Branch("SignalNcr3", &SignalNcr3);
//  fTree2->Branch("tpc_ncls", &tpc_ncls);
//  fTree2->Branch("its_ncls", &its_ncls);
//  fTree2->Branch("tpc_chi2", &tpc_chi2);
//  fTree2->Branch("its_chi2", &its_chi2);
//  fTree2->Branch("dca_r",&dca_r);
//  fTree2->Branch("dca_z",&dca_z);
//  fTree2->Branch("dca_tpc_r",&dca_tpc_r);
//  fTree2->Branch("dca_tpc_z",&dca_tpc_z);
//  fTree2->Branch("fITSClusterMap",&fITSClusterMap);
//  fTree2->Branch("fITSClusterMap_leg",&fITSClusterMap_leg);
//  fTree2->Branch("fAlpha",&fAlpha);
//  fTree2->Branch("fSigned1Pt",&fSigned1Pt);
//  fTree2->Branch("selectionPtMask",&selectionPtMask);
//  fTree2->Branch("selectionPIDMask",&selectionPIDMask);
  //  THnSparseF* dEdx_THnF = new THnSparseF(name.Data(), title.Data(), kNdim, binsHistQA, xminHistQA, xmaxHistQA);
  // BinLogAxis(h, 0);
}
void AliSkimmedDataAnalysisMaker::TreeV0_BBFitAnalysis(TString filename_v0) {

  cout << "!!!!!!!!!!!!!!!!!read" << filename_v0 << endl;
  TFile *infile2 = new TFile(filename_v0.Data());
  TTree *Tree2 = (TTree *) infile2->Get("V0Flat");
  if (Tree2 == 0x0) {
  } else {
    V0FlatAna *V0ana = new V0FlatAna();
    V0ana->Init(Tree2);

    for (Int_t iV0 = 0; iV0 < Tree2->GetEntries(); iV0++) {
      if (iV0 % 10000 == 0)
        cout << "working on Event  " << iV0 << endl;
      V0ana->GetEntry(iV0);

      Int_t particleID = -1;
      Double_t tpccluster = -999;
      //Double_t dca_r = -999;
      //Double_t dca_z = -999;
      Double_t nCrossRows = -999;
      tgl = -999;
      tofnsigma = -999;
      tpcnsigma = -999;
      shiftM = -999;
      nPileUpPrim = -999;
      nPileUpSum = -999;
      primMult = -999;
      tpcClusterMult = -999;
      pileUpOK = -999;
      multSSD = -999;
      multSDD = -999;
      SignalTot0 = -999;
      SignalTot1 = -999;
      SignalTot2 = -999;
      SignalTot3 = -999;
      SignalMax0 = -999;
      SignalMax1 = -999;
      SignalMax2 = -999;
      SignalMax3 = -999;
      SignalNcl0 = -999;
      SignalNcl1 = -999;
      SignalNcl2 = -999;
      SignalNcl3 = -999;
      SignalNcr0 = -999;
      SignalNcr1 = -999;
      SignalNcr2 = -999;
      SignalNcr3 = -999;
      fSigned1Pt = -999;
      selectionPtMask = -999;
      selectionPIDMask = -999;

      const Int_t numCases = 8;
      Double_t tpcQA[numCases];
      Double_t tofQA[numCases];
      for (Int_t i = 0; i < numCases; i++) {
        tpcQA[i] = -999.;
        tofQA[i] = -999.;
      }

      for (Int_t track = 0; track < 2; track++) {


        bool isTreeK0 = false;
        bool isTreeGamma = false;
        bool isTreeLambda = false;
        bool isTreeALambda = false;
        tofnsigma = -999;
        tpcnsigma = -999;

        if (V0ana->K0Like > 0.7 && (V0ana->cleanK0 == 1)) {
          isTreeK0 = true;
          V0Pull = V0ana->K0Pull;
          V0PullEff = V0ana->K0PullEff;
          isV0 = 1.;
        }
        if (V0ana->ELike > 0.7 && (V0ana->cleanGamma == 1)) {
          isTreeGamma = true;
          V0Pull = V0ana->EPull;
          V0PullEff = V0ana->EPullEff;
          isV0 = 2.;
        }

        if (V0ana->LLike > 0.7 && (V0ana->cleanLambda == 1)) {
          isTreeLambda = true;
          V0Pull = V0ana->LPull;
          V0PullEff = V0ana->LPullEff;
          isV0 = 3;
        }

        if (V0ana->ALLike > 0.7 && (V0ana->cleanALambda == 1)) {
          isTreeALambda = true;
          V0Pull = V0ana->ALPull;
          V0PullEff = V0ana->ALPullEff;
          isV0 = 4;
        }

        if (!(isTreeK0 || isTreeGamma || isTreeLambda || isTreeALambda)) continue;

        if (track == 0) {

          tpccluster = V0ana->track0_fTPCsignalN;
          p = V0ana->track0P;
          if (V0ana->track0Pt == 0) {
            oneoverpt = 100000000;
          } else {
            oneoverpt = 1 / V0ana->track0Pt;
          }
          rawtpcsignal = V0ana->track0_fTPCsignal;
          rawtpcsignal_leg = V0ana->track1_fTPCsignal;
          itssignal = V0ana->track0_fITSsignal;
          pileupcor_tpcsignal = V0ana->track0_fTPCsignal - V0ana->track0GetPileupValue;
          pileupcor_tpcsignal_leg = V0ana->track1_fTPCsignal - V0ana->track1GetPileupValue;
          eta = V0ana->track0Eta;
          tgl = V0ana->track0Tgl;
          multiplicty = V0ana->tpcTrackBeforeClean;
          centrality = V0ana->centV0;
          //isV0 = 1.;
          isPileUp = V0ana->isPileUp;
          dca_r = V0ana->track0DCAxy;
          dca_z = V0ana->track0DCAz;
          nCrossRows = V0ana->track0ncrossrow;
          fSigned1Pt = V0ana->fSigned1Pt0;

          shiftM = V0ana->shiftM;
          nPileUpPrim = V0ana->nPileUpPrim;
          nPileUpSum = V0ana->nPileUpSum;
          primMult = V0ana->primMult;
          tpcClusterMult = V0ana->tpcClusterMult;
          pileUpOK = V0ana->pileUpOK;
          multSSD = V0ana->multSSD;
          multSDD = V0ana->multSDD;
          SignalTot0 = V0ana->logSignalTot_track0_0;
          SignalTot1 = V0ana->logSignalTot_track0_1;
          SignalTot2 = V0ana->logSignalTot_track0_2;
          SignalTot3 = V0ana->logSignalTot_track0_3;

          SignalMax0 = V0ana->logSignalMax_track0_0;
          SignalMax1 = V0ana->logSignalMax_track0_1;
          SignalMax2 = V0ana->logSignalMax_track0_2;
          SignalMax3 = V0ana->logSignalMax_track0_3;

          SignalNcl0 = V0ana->signalNcl_track0_0;
          SignalNcl1 = V0ana->signalNcl_track0_1;
          SignalNcl2 = V0ana->signalNcl_track0_2;
          SignalNcl3 = V0ana->signalNcl_track0_3;

          SignalNcr0 = V0ana->signalNcr_track0_0;
          SignalNcr1 = V0ana->signalNcr_track0_1;
          SignalNcr2 = V0ana->signalNcr_track0_2;
          SignalNcr3 = V0ana->signalNcr_track0_3;

          tpc_ncls = V0ana->track0TPC_cls;
          its_ncls = V0ana->track0ITS_cls;
          tpc_chi2 = V0ana->track0TPCchi2;
          its_chi2 = V0ana->track0ITSchi2;
          dca_r = V0ana->track0DCAxy;
          dca_z = V0ana->track0DCAz;
          dca_tpc_r = V0ana->track0DCAxy_TPC;
          dca_tpc_z = V0ana->track0DCAz_TPC;
          fITSClusterMap = V0ana->fITSClusterMap0;
          fITSClusterMap_leg = V0ana->fITSClusterMap1;
          fAlpha = V0ana->fAlpha0;
          K0PullEff = V0ana->K0PullEff;
          EPullEff = V0ana->EPullEff;
          LPullEff = V0ana->LPullEff;
          ALPullEff = V0ana->ALPullEff;
          // New information for track0
          trackPtComb=1/TMath::Abs(V0ana->track0CombP4);
          tglComb=V0ana->track0Tgl;
          trackPComb=trackPtComb*TMath::Sqrt(1+tglComb*tglComb);
          // local phi information
          phi_ROC0=V0ana->phi0_ROC0;
          phi_ROC1=V0ana->phi0_ROC1;
          phi_ROC2=V0ana->phi0_ROC2;
          // sector information if exist
          dSector0=0;
          dSector1=0;
          dSector2=0;
          //MC true information if exist
          ptTPCMC=0;
          ptMC=0;
          tglMC=0;
          fPdgCodeMC=0;


          tpcQA[0] = V0ana->track0tpcNsigma_el;
          tofQA[0] = V0ana->track0tofNsigmaElectron;
          tpcNsig[0] = V0ana->track0tpcNsigma_el;
          tofNsig[0] = V0ana->track0tofNsigmaElectron;

          tpcQA[1] = V0ana->track0tpcNsigma_pi;
          tofQA[1] = V0ana->track0tofNsigmaPion;
          tpcNsig[2] = V0ana->track0tpcNsigma_pi;
          tofNsig[2] = V0ana->track0tofNsigmaPion;

          tpcQA[2] = V0ana->track0tpcNsigma_ka;
          tofQA[2] = V0ana->track0tofNsigmaKaon;
          tpcNsig[3] = V0ana->track0tpcNsigma_ka;
          tofNsig[3] = V0ana->track0tofNsigmaKaon;

          tpcQA[3] = V0ana->track0tpcNsigma_pro;
          tofQA[3] = V0ana->track0tofNsigmaProton;
          tpcNsig[4] = V0ana->track0tpcNsigma_pro;
          tofNsig[4] = V0ana->track0tofNsigmaProton;

          tpcQA[4] = V0ana->track0tpcNsigma_ka;
          tofQA[4] = V0ana->track0tofNsigmaKaon;
          tpcNsig[3] = V0ana->track0tpcNsigma_ka;
          tofNsig[3] = V0ana->track0tofNsigmaKaon;

          tpcQA[5] = -999;
          tofQA[5] = -999;
          tpcNsig[5] = -999;
          tofNsig[5] = -999;

          tpcQA[6] = -999;
          tofQA[6] = -999;
          tpcNsig[6] = -999;
          tofNsig[6] = -999;

          tpcQA[7] = -999;
          tofQA[7] = -999;

          tofnsigma = -999;
          tpcnsigma = -999;

        } else if (track == 1) {

          tpccluster = V0ana->track1_fTPCsignalN;
          p = V0ana->track1P;
          if (V0ana->track1Pt == 0) {
            oneoverpt = 100000000;
          } else {
            oneoverpt = 1 / V0ana->track1Pt;
          }
          rawtpcsignal = V0ana->track1_fTPCsignal;
          rawtpcsignal_leg = V0ana->track0_fTPCsignal;
          itssignal = V0ana->track1_fITSsignal;
          pileupcor_tpcsignal = V0ana->track1_fTPCsignal - V0ana->track1GetPileupValue;
          pileupcor_tpcsignal_leg = V0ana->track0_fTPCsignal - V0ana->track0GetPileupValue;
          eta = V0ana->track1Eta;
          tgl = V0ana->track1Tgl;
          multiplicty = V0ana->tpcTrackBeforeClean;
          centrality = V0ana->centV0;
          //isV0 = 1.;
          isPileUp = V0ana->isPileUp;
          dca_r = V0ana->track1DCAxy;
          dca_z = V0ana->track1DCAz;
          nCrossRows = V0ana->track1ncrossrow;
          fSigned1Pt = V0ana->fSigned1Pt1;

          shiftM = V0ana->shiftM;
          nPileUpPrim = V0ana->nPileUpPrim;
          nPileUpSum = V0ana->nPileUpSum;
          primMult = V0ana->primMult;
          tpcClusterMult = V0ana->tpcClusterMult;
          pileUpOK = V0ana->pileUpOK;
          multSSD = V0ana->multSSD;
          multSDD = V0ana->multSDD;

          SignalTot0 = V0ana->logSignalTot_track1_0;
          SignalTot1 = V0ana->logSignalTot_track1_1;
          SignalTot2 = V0ana->logSignalTot_track1_2;
          SignalTot3 = V0ana->logSignalTot_track1_3;

          SignalMax0 = V0ana->logSignalMax_track1_0;
          SignalMax1 = V0ana->logSignalMax_track1_1;
          SignalMax2 = V0ana->logSignalMax_track1_2;
          SignalMax3 = V0ana->logSignalMax_track1_3;

          SignalNcl0 = V0ana->signalNcl_track1_0;
          SignalNcl1 = V0ana->signalNcl_track1_1;
          SignalNcl2 = V0ana->signalNcl_track1_2;
          SignalNcl3 = V0ana->signalNcl_track1_3;

          SignalNcr0 = V0ana->signalNcr_track1_0;
          SignalNcr1 = V0ana->signalNcr_track1_1;
          SignalNcr2 = V0ana->signalNcr_track1_2;
          SignalNcr3 = V0ana->signalNcr_track1_3;

          tpc_ncls = V0ana->track1TPC_cls;
          its_ncls = V0ana->track1ITS_cls;
          tpc_chi2 = V0ana->track1TPCchi2;
          its_chi2 = V0ana->track1ITSchi2;
          dca_r = V0ana->track1DCAxy;
          dca_z = V0ana->track1DCAz;
          dca_tpc_r = V0ana->track1DCAxy_TPC;
          dca_tpc_z = V0ana->track1DCAxy_TPC;
          fITSClusterMap = V0ana->fITSClusterMap1;
          fITSClusterMap_leg = V0ana->fITSClusterMap0;
          fAlpha = V0ana->fAlpha1;
          K0PullEff = V0ana->K0PullEff;
          EPullEff = V0ana->EPullEff;
          LPullEff = V0ana->LPullEff;
          ALPullEff = V0ana->ALPullEff;
           // New information for track1
          trackPtComb=1/TMath::Abs(V0ana->track1CombP4);
          tglComb=V0ana->track1Tgl;
          trackPComb=trackPtComb*TMath::Sqrt(1+tglComb*tglComb);
          // local phi information
          phi_ROC0=V0ana->phi1_ROC0;
          phi_ROC1=V0ana->phi1_ROC1;
          phi_ROC2=V0ana->phi1_ROC2;
          // sector information if exist
          dSector0=0;
          dSector1=0;
          dSector2=0;
          //MC true information if exist
          ptTPCMC=0;
          ptMC=0;
          tglMC=0;
          fPdgCodeMC=0;
          //
          dZ_Row0=V0ana->dZ_Row0;
          dRPhi_Row0=V0ana->dRPhi_Row0;
          dZ_Row62=V0ana->dZ_Row62;
          dRPhi_Row62=V0ana->dRPhi_Row62;
          dZ_Row125=V0ana->dZ_Row125;
          dRPhi_Row125=V0ana->dRPhi_Row125;
          dZ_Row158=V0ana->dZ_Row158;
          dRPhi_Row158=V0ana->dRPhi_Row158;
          v0_fRr=V0ana->v0_fRr;

          tpcQA[0] = V0ana->track1tpcNsigma_el;
          tofQA[0] = V0ana->track1tofNsigmaElectron;
          tpcNsig[0] = V0ana->track1tpcNsigma_el;
          tofNsig[0] = V0ana->track1tofNsigmaElectron;

          tpcQA[1] = V0ana->track1tpcNsigma_pi;
          tofQA[1] = V0ana->track1tofNsigmaPion;
          tpcNsig[2] = V0ana->track1tpcNsigma_pi;
          tofNsig[2] = V0ana->track1tofNsigmaPion;

          tpcQA[2] = V0ana->track1tpcNsigma_ka;
          tofQA[2] = V0ana->track1tofNsigmaKaon;
          tpcNsig[3] = V0ana->track1tpcNsigma_ka;
          tofNsig[3] = V0ana->track1tofNsigmaKaon;

          tpcQA[3] = V0ana->track1tpcNsigma_pro;
          tofQA[3] = V0ana->track1tofNsigmaProton;
          tpcNsig[4] = V0ana->track1tpcNsigma_pro;
          tofNsig[4] = V0ana->track1tofNsigmaProton;

          tpcQA[4] = V0ana->track1tpcNsigma_ka;
          tofQA[4] = V0ana->track1tofNsigmaKaon;
          tpcNsig[3] = V0ana->track1tpcNsigma_ka;
          tofNsig[3] = V0ana->track1tofNsigmaKaon;

          tpcQA[5] = -999;
          tofQA[5] = -999;
          tpcNsig[5] = -999;
          tofNsig[5] = -999;

          tpcQA[6] = -999;
          tofQA[6] = -999;
          tpcNsig[6] = -999;
          tofNsig[6] = -999;

          tpcQA[7] = -999;
          tofQA[7] = -999;

          tofnsigma = -999;
          tpcnsigma = -999;

        }


        if (tpccluster < 50) continue;
        //if(TMath::Abs(dca_r) > 3.0) continue;
        //if(TMath::Abs(dca_z) > 3.0) continue;
        if (nCrossRows < 70) continue;

        if (isTreeK0) {
          particleID = 2; //pion
          PDGcode = 2;
          BG = p / 1.39569997787475586e-01;
          tofnsigma = tofQA[1];
          tpcnsigma = tpcQA[1];

        } else if (isTreeGamma) {
          particleID = 1; // electron
          PDGcode = 0;//AliPID::kElectron;
          BG = p / 5.10998885147273540e-04;
          tofnsigma = tofQA[0];
          tpcnsigma = tpcQA[0];

        } else if (isTreeLambda) {
          if (track == 0) {
            particleID = 3; // proton
            PDGcode = 4.;
            BG = p / 9.38271999359130859e-01;
            tofnsigma = tofQA[3];
            tpcnsigma = tpcQA[3];
          } else {
            particleID = 2; // pion
            PDGcode = 2;//AliPID::kPion;
            BG = p / 1.39569997787475586e-01;
            tofnsigma = tofQA[1];
            tpcnsigma = tpcQA[1];
          }
        } else if (isTreeALambda) {
          if (track == 0) {
            particleID = 2; // pion
            PDGcode = 2.;//AliPID::kPion;
            BG = p / 1.39569997787475586e-01;
            tofnsigma = tofQA[1];
            tpcnsigma = tpcQA[1];

          } else {
            particleID = 3; // proton
            PDGcode = 4;;
            BG = p / 9.38271999359130859e-01;
            tofnsigma = tofQA[3];
            tpcnsigma = tpcQA[3];
          }
        }

        fTree2->Fill();

      }    // end track loop

    } // end treeloop
  }

} //end function


void AliSkimmedDataAnalysisMaker::TreePrimary_BBFitAnalysis(TString filename_track) {

  cout << "!!!!!!!!!!!!!!!!!read!!!!!!!" << filename_track << endl;
  TFile *infile = new TFile(filename_track.Data());
  TTree *Tree = (TTree *) infile->Get("CleanTrackFlat");
  if (Tree == 0x0) {
  } else {
    TrackFlatAna *TrackAna = new TrackFlatAna();
    TrackAna->Init(Tree);

    for (Int_t iV0 = 0; iV0 < Tree->GetEntries(); iV0++) {
      if (iV0 % 10000 == 0)
        cout << "working on Event  " << iV0 << endl;
      TrackAna->GetEntry(iV0);

      const Int_t numCases = 8;
      Double_t tpcQA[numCases];
      Double_t tofQA[numCases];
      for (Int_t i = 0; i < numCases; i++) {
        tpcQA[i] = -999.;
        tofQA[i] = -999.;
      }

      Double_t tpccluster = -999;
      p = -999;
      oneoverpt = -999;
      rawtpcsignal = -999;
      itssignal = -999;
      pileupcor_tpcsignal = -999;
      eta = -999;
      multiplicty = -999;
      centrality = -999;
      isV0 = 0.;
      PDGcode = -999;
      tgl = -999;
      BG = -999;
      tofnsigma = -999;
      tpcnsigma = -999;
      //Double_t dca_r = -999;
      //Double_t dca_z = -999;
      Double_t nCrossRows = -999;
      fSigned1Pt = -999;
      selectionPtMask = -999;
      selectionPIDMask = -999;

      shiftM = -999;
      nPileUpPrim = -999;
      nPileUpSum = -999;
      primMult = -999;
      tpcClusterMult = -999;
      pileUpOK = -999;
      multSSD = -999;
      multSDD = -999;
      SignalTot0 = -999;
      SignalTot1 = -999;
      SignalTot2 = -999;
      SignalTot3 = -999;
      SignalMax0 = -999;
      SignalMax1 = -999;
      SignalMax2 = -999;
      SignalMax3 = -999;
      SignalNcl0 = -999;
      SignalNcl1 = -999;
      SignalNcl2 = -999;
      SignalNcl3 = -999;
      SignalNcr0 = -999;
      SignalNcr1 = -999;
      SignalNcr2 = -999;
      SignalNcr3 = -999;
      tpc_ncls = -999;
      its_ncls = -999;
      tpc_chi2 = -999;
      its_chi2 = -999;
      V0Pull = -999;
      V0PullEff = -999;
      K0PullEff = -999;
      EPullEff = -999;
      LPullEff = -999;
      ALPullEff = -999;

      tpccluster = TrackAna->esdTrack_fTPCsignalN;
      eta = TrackAna->trackEta;
      //dca_r = TrackAna->dca_r;
      //dca_z = TrackAna->dca_z;
      nCrossRows = TrackAna->nCrossRows;

      if (tpccluster < 50) continue;
      if (TMath::Abs(eta) > 0.9) continue;
      //if(TMath::Abs(dca_r) > 3.0) continue;
      //if(TMath::Abs(dca_z) > 3.0) continue;
      if (nCrossRows < 70) continue;

      p = TrackAna->trackP;
      oneoverpt = 1 / TrackAna->trackPt;
      rawtpcsignal = TrackAna->fTPCsignal;
      itssignal = TrackAna->esdTrack_fITSsignal;
      pileupcor_tpcsignal = TrackAna->fTPCsignal - TrackAna->track_GetPileupValue;
      multiplicty = TrackAna->tpcTrackBeforeClean;
      centrality = TrackAna->centV0;
      isPileUp = TrackAna->isPileUp;
      tgl = TrackAna->tgl;
      fSigned1Pt = TrackAna->fSigned1Pt;
      selectionPtMask = TrackAna->selectionPtMask;
      selectionPIDMask = TrackAna->selectionPIDMask;

      shiftM = TrackAna->shiftM;
      nPileUpPrim = TrackAna->nPileUpPrim;
      nPileUpSum = TrackAna->nPileUpSum;
      primMult = TrackAna->primMult;
      tpcClusterMult = TrackAna->tpcClusterMult;
      pileUpOK = TrackAna->pileUpOK;
      multSSD = TrackAna->multSSD;
      multSDD = TrackAna->multSDD;
      SignalTot0 = TrackAna->logSignalTot0;
      SignalTot1 = TrackAna->logSignalTot1;
      SignalTot2 = TrackAna->logSignalTot2;
      SignalTot3 = TrackAna->logSignalTot3;

      SignalMax0 = TrackAna->logSignalMax0;
      SignalMax1 = TrackAna->logSignalMax1;
      SignalMax2 = TrackAna->logSignalMax2;
      SignalMax3 = TrackAna->logSignalMax3;

      SignalNcl0 = TrackAna->signalNcl0;
      SignalNcl1 = TrackAna->signalNcl1;
      SignalNcl2 = TrackAna->signalNcl2;
      SignalNcl3 = TrackAna->signalNcl3;

      SignalNcr0 = TrackAna->signalNcr0;
      SignalNcr1 = TrackAna->signalNcr1;
      SignalNcr2 = TrackAna->signalNcr2;
      SignalNcr3 = TrackAna->signalNcr3;

      tpc_ncls = TrackAna->trackTPC_cls;
      its_ncls = TrackAna->trackITS_cls;
      tpc_chi2 = TrackAna->trackTPCchi2;
      its_chi2 = TrackAna->trackITSchi2;
      dca_r = TrackAna->trackDCAxy;
      dca_z = TrackAna->trackDCAz;
      dca_tpc_r = TrackAna->trackDCAxy_TPC;
      dca_tpc_z = TrackAna->trackDCAz_TPC;
      fITSClusterMap = TrackAna->fITSClusterMap;
      fAlpha = TrackAna->fAlpha;

      /// PDG/PID code is assigned to the tracks based on some varaibles - thi segment was done originally for QA - not well readable
      tpcQA[0] = TrackAna->track_tpcNsigma_el;
      tofQA[0] = TrackAna->tracktofNsigmaElectron;
      tpcNsig[0] = TrackAna->track_tpcNsigma_el;
      tofNsig[0] = TrackAna->tracktofNsigmaElectron;

      tpcQA[1] = TrackAna->track_tpcNsigma_pi;
      tofQA[1] = TrackAna->tracktofNsigmaPion;
      tpcNsig[2] = TrackAna->track_tpcNsigma_pi;
      tofNsig[2] = TrackAna->tracktofNsigmaPion;

      tpcQA[2] = TrackAna->track_tpcNsigma_ka;
      tofQA[2] = TrackAna->tracktofNsigmaKaon;
      tpcNsig[3] = TrackAna->track_tpcNsigma_ka;
      tofNsig[3] = TrackAna->tracktofNsigmaKaon;

      tpcQA[3] = TrackAna->track_tpcNsigma_pro;
      tofQA[3] = TrackAna->tracktofNsigmaProton;
      tpcNsig[4] = TrackAna->track_tpcNsigma_pro;
      tofNsig[4] = TrackAna->tracktofNsigmaProton;

      tpcQA[4] = TrackAna->track_tpcNsigma_ka;
      tofQA[4] = TrackAna->tracktofNsigmaKaon;
      tpcNsig[3] = TrackAna->track_tpcNsigma_ka;
      tofNsig[3] = TrackAna->tracktofNsigmaKaon;

      tpcQA[5] = TrackAna->track_tpcNsigma_deut;
      tofQA[5] = TrackAna->tracktofNsigmaDeuteron;
      tpcNsig[5] = TrackAna->track_tpcNsigma_deut;
      tofNsig[5] = TrackAna->tracktofNsigmaDeuteron;

      tpcQA[6] = TrackAna->track_tpcNsigma_tri;
      tofQA[6] = TrackAna->tracktofNsigmaTrition;
      tpcNsig[6] = TrackAna->track_tpcNsigma_tri;
      tofNsig[6] = TrackAna->tracktofNsigmaTrition;

      //assign PDG code to different particles:


      // pion - (TOF 3 sigma & TPC4 sigma in TOF region  p(0.5-1.) or TPC 4 sigma
      if ((p > 0.5 && p < 1.1 && TMath::Abs(tofQA[1]) < 3.0 && TMath::Abs(tpcQA[1]) < 4.0) || p < 0.5 && TMath::Abs(tpcQA[1]) < 4.0) {
        PDGcode = 2;
        BG = p / 1.39569997787475586e-01;
        tofnsigma = tofQA[1];
        tpcnsigma = tpcQA[1];
      }



      //proton
      if ((p < 3.0 && p > 0.6 && TMath::Abs(tofQA[3]) < 3.0 && (rawtpcsignal >= 50. / TMath::Power(p, 1.3))) || (p < 0.6) && (rawtpcsignal >= 50. / TMath::Power(p, 1.3))) {
        if (p > 1.0) {
          if (TrackAna->Nucleitrigger_OFF == 0) continue;
        }
        PDGcode = 4;
        BG = p / 9.38271999359130859e-01;
        tofnsigma = tofQA[3];
        tpcnsigma = tpcQA[3];
      }


      //Deuteron  (TPC only below 1.5 GeV ) || TOF on p (1.5-2 GeV)  -
      if ((p < 1.5 && TMath::Abs(tpcQA[5]) < 5.0) || (p < 2.0 && p > 1.5 && TMath::Abs(tofQA[5]) < 3.0 && TMath::Abs(tpcQA[5]) < 5.0)) {
        PDGcode = 5;
        BG = p / 1.87561297416687012e+00;
        tofnsigma = tofQA[5];
        tpcnsigma = tpcQA[5];
      }

//Trition
      if ((p < 1.5 && TMath::Abs(tpcQA[6]) < 5.0) || (p < 2.0 && p > 1.5 && TMath::Abs(tofQA[6]) < 3.0 && TMath::Abs(tpcQA[6]) < 5.0)) {
        PDGcode = 6;
        BG = p / 2.80925011634826660e+00;
        tofnsigma = tofQA[6];
        tpcnsigma = tpcQA[6];
      }


      //kaon
      if ((p > 0.4 && p < 1.5 && TMath::Abs(tofQA[2]) < 3.0 && TMath::Abs(tpcQA[2]) < 5.0) || (p < 0.4 && TMath::Abs(tpcQA[2]) < 5.0)) {
        PDGcode = 3;
        BG = p / 4.93676990270614624e-01;
        tofnsigma = tofQA[2];
        tpcnsigma = tpcQA[2];
      }

      if (PDGcode == -999) continue;
      // new information
      // comined information
      trackPComb = TrackAna->trackPComb;
      trackPtComb = TrackAna->trackPtComb;
      tglComb = TrackAna->tglComb;
      // local phi information
      phi_ROC0 = TrackAna->phi_ROC0;
      phi_ROC1 = TrackAna->phi_ROC1;
      phi_ROC2 = TrackAna->phi_ROC2;
      // sector information if exist
      dSector0 = TrackAna->dSector0;
      dSector1 = TrackAna->dSector1;
      dSector2 = TrackAna->dSector2;
      //MC true information if exist
      ptTPCMC = TrackAna->ptTPCMC;
      ptMC = TrackAna->ptMC;
      tglMC = TrackAna->tglMC;
      fPdgCodeMC = TrackAna->fPdgCode;


      fTree2->Fill();
    }
  }

}

void AliSkimmedDataAnalysisMaker::WriteHistogram()
{

      cout<<"  write"<<endl;

    mOutputFile_BBFitAnalysis->cd();
    fTree2->Write();

    mOutputFile_BBFitAnalysis->Close();
}
  
void AliSkimmedDataAnalysisMaker::BinLogAxis(THnSparseF *h, Int_t axisNumber)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];

  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;

}
void AliSkimmedDataAnalysisMaker::SetAxisNamesFromTitle(const THnSparseF *h)
{
  // set the histogram axis names from the axis titles
  for (Int_t i=0; i<h->GetNdimensions(); ++i) {
    h->GetAxis(i)->SetName(h->GetAxis(i)->GetTitle());
  }
}

Double_t AliSkimmedDataAnalysisMaker::Lund(Double_t* xx, Double_t* par)
{
  // bg is beta*gamma
  const Double_t bg = xx[0];

  const Double_t beta2 = bg*bg / (1.0 + bg*bg);

  const Double_t a = par[0];
  const Double_t b = par[1];
  const Double_t c = par[2];
  const Double_t e = par[3];
  const Double_t f = par[4];

  const Double_t d = TMath::Exp(c*(a - f) / b);

  const Double_t powbg = TMath::Power(1.0 + bg, c);

  const Double_t value = a / TMath::Power(beta2,e) +
    b/c * TMath::Log(powbg / (1.0 + d*powbg));

  return value;
}


//________________________________________________________________________
Double_t AliSkimmedDataAnalysisMaker::SaturatedLund(Double_t* xx, Double_t* par)
{
  const Double_t qq = Lund(xx, par);
  return qq * TMath::Exp(par[5] / qq);
}

Double_t AliSkimmedDataAnalysisMaker::Aleph(Double_t* xx, Double_t* par) {
  const Double_t bg = xx[0];

  const Double_t a0 = par[0];
  const Double_t a1 = par[1];
  const Double_t a2 = par[2];
  const Double_t a3 = par[3];
  const Double_t a4 = par[4];
  const Double_t a5 = par[5];

  const Double_t beta = TMath::Sqrt(bg*bg / (1.0 + bg*bg));
  const Double_t powbetaa3 = TMath::Power(beta,a3);

  const Double_t value = a0/powbetaa3 * (a1 - a2 - a5 * powbetaa3 - TMath::Log(1.0 + TMath::Power(bg, -a4)*TMath::Exp(-a2)));

  return value;
}
