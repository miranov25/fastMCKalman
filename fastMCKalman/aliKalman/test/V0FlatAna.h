//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 21 11:14:00 2021 by ROOT version 6.20/08
// from TTree V0Flat/V0Flat
// found on file: V0tree.root
//////////////////////////////////////////////////////////

#ifndef V0FlatAna_h
#define V0FlatAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class V0FlatAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        run;
   Double_t        intrate;
   Double_t        timeStampS;
   Double_t        timestamp;
   Double_t        bField;
   Double_t        selectionPtMask;
   Double_t        gid;
   Double_t        shiftA;
   Double_t        shiftC;
   Double_t        shiftM;
   Double_t        nPileUpPrim;
   Double_t        nPileUpSum;
   Double_t        nPileUpSumCorr;
   Double_t        primMult;
   Double_t        tpcClusterMult;
   Double_t        pileUpOK;
   Double_t        v0_fPointAngle;
   Double_t        kf_fChi2;
   Double_t        K0Like;
   Double_t        ELike;
   Double_t        LLike;
   Double_t        ALLike;
   Double_t        cleanK0;
   Double_t        cleanGamma;
   Double_t        cleanLambda;
   Double_t        cleanALambda;
   Double_t        K0Pull;
   Double_t        LPull;
   Double_t        ALPull;
   Double_t        EPull;
   Double_t        K0PullEff;
   Double_t        LPullEff;
   Double_t        ALPullEff;
   Double_t        EPullEff;
   Double_t        centV0;
   Double_t        centITS0;
   Double_t        centITS1;
   Double_t        multSSD;
   Double_t        multSDD;
   Double_t        tpcTrackBeforeClean;
   Double_t        triggerMask;
   Double_t        pairDCA_InerTPC;
   Double_t        isMinBias;
   Double_t        isPileUp;
   Double_t        track0_fTPCsignal;
   Double_t        track1_fTPCsignal;
   Double_t        track0_fTPCsignalN;
   Double_t        track1_fTPCsignalN;
   Double_t        type;
   Double_t        track0_fITSsignal;
   Double_t        track1_fITSsignal;
   Double_t        track0status;
   Double_t        track1status;
   Double_t        track0_hasTOF;
   Double_t        track1_hasTOF;
   Double_t        track0P;
   Double_t        track0Pt;
   Double_t        track0Eta;
   Double_t        track0Phi;
   Double_t        track0Px;
   Double_t        track0Py;
   Double_t        track0Pz;
   Double_t        track0Tgl;
   Double_t        track0DCAxy;
   Double_t        track0DCAz;
   Double_t        track0ncrossrow;
   Double_t        track0TPCchi2;
   Double_t        track0ITSchi2;
   Double_t        track0TPC_cls;
   Double_t        track0ITS_cls;
   Double_t        track1P;
   Double_t        track1Pt;
   Double_t        track1Eta;
   Double_t        track1Phi;
   Double_t        track1Px;
   Double_t        track1Py;
   Double_t        track1Pz;
   Double_t        track1Tgl;
   Double_t        track1DCAxy;
   Double_t        track1DCAz;
   Double_t        track1ncrossrow;
   Double_t        track1TPCchi2;
   Double_t        track1ITSchi2;
   Double_t        track1TPC_cls;
   Double_t        track1ITS_cls;
   Double_t        track0DCAxy_TPC;
   Double_t        track0DCAz_TPC;
   Double_t        fAlpha0;
   Double_t        fITSClusterMap0;
   Double_t        fSigned1Pt0;
   Double_t        track1DCAxy_TPC;
   Double_t        track1DCAz_TPC;
   Double_t        fAlpha1;
   Double_t        fITSClusterMap1;
   Double_t        fSigned1Pt1;
   Double_t        track0tofNsigmaElectron;
   Double_t        track0tofNsigmaPion;
   Double_t        track0tofNsigmaKaon;
   Double_t        track0tofNsigmaProton;
   Double_t        track1tofNsigmaElectron;
   Double_t        track1tofNsigmaPion;
   Double_t        track1tofNsigmaKaon;
   Double_t        track1tofNsigmaProton;
   Double_t        track0tpcNsigma_el;
   Double_t        track0tpcNsigma_pi;
   Double_t        track0tpcNsigma_ka;
   Double_t        track0tpcNsigma_pro;
   Double_t        track1tpcNsigma_el;
   Double_t        track1tpcNsigma_pi;
   Double_t        track1tpcNsigma_ka;
   Double_t        track1tpcNsigma_pro;
   Double_t        track0ExpectedTPCSignalV0_el;
   Double_t        track0ExpectedTPCSignalV0_pi;
   Double_t        track0ExpectedTPCSignalV0_ka;
   Double_t        track0ExpectedTPCSignalV0_pro;
   Double_t        track1ExpectedTPCSignalV0_el;
   Double_t        track1ExpectedTPCSignalV0_pi;
   Double_t        track1ExpectedTPCSignalV0_ka;
   Double_t        track1ExpectedTPCSignalV0_pro;
   Double_t        track0CorrectedTPCSignalV0_el;
   Double_t        track0CorrectedTPCSignalV0_pi;
   Double_t        track0CorrectedTPCSignalV0_ka;
   Double_t        track0CorrectedTPCSignalV0_pro;
   Double_t        track1CorrectedTPCSignalV0_el;
   Double_t        track1CorrectedTPCSignalV0_pi;
   Double_t        track1CorrectedTPCSignalV0_ka;
   Double_t        track1CorrectedTPCSignalV0_pro;
   Double_t        track0GetPileupValue;
   Double_t        track1GetPileupValue;
   Double_t        track0tpcNsigma_no_corrected_el;
   Double_t        track0tpcNsigma_no_corrected_pi;
   Double_t        track0tpcNsigma_no_corrected_pr;
   Double_t        track1tpcNsigma_no_corrected_el;
   Double_t        track1tpcNsigma_no_corrected_pi;
   Double_t        track1tpcNsigma_no_corrected_pr;
   Double_t        track0ExpectedTPCSignalV0_no_corrected_el;
   Double_t        track0ExpectedTPCSignalV0_no_corrected_pi;
   Double_t        track0ExpectedTPCSignalV0_no_corrected_pr;
   Double_t        track1ExpectedTPCSignalV0_no_corrected_el;
   Double_t        track1ExpectedTPCSignalV0_no_corrected_pi;
   Double_t        track1ExpectedTPCSignalV0_no_corrected_pr;
   Double_t        track0tpcNsigma_corrected_eta_el;
   Double_t        track0tpcNsigma_corrected_eta_pi;
   Double_t        track0tpcNsigma_corrected_eta_pr;
   Double_t        track1tpcNsigma_corrected_eta_el;
   Double_t        track1tpcNsigma_corrected_eta_pi;
   Double_t        track1tpcNsigma_corrected_eta_pr;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_el;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_pi;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_pr;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_el;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_pi;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_pr;
   Double_t        track0tpcNsigma_corrected_multip_el;
   Double_t        track0tpcNsigma_corrected_multip_pi;
   Double_t        track0tpcNsigma_corrected_multip_pr;
   Double_t        track1tpcNsigma_corrected_multip_el;
   Double_t        track1tpcNsigma_corrected_multip_pi;
   Double_t        track1tpcNsigma_corrected_multip_pr;
   Double_t        track0ExpectedTPCSignalV0_corrected_multip_el;
   Double_t        track0ExpectedTPCSignalV0_corrected_multip_pi;
   Double_t        track0ExpectedTPCSignalV0_corrected_multip_pr;
   Double_t        track1ExpectedTPCSignalV0_corrected_multip_el;
   Double_t        track1ExpectedTPCSignalV0_corrected_multip_pi;
   Double_t        track1ExpectedTPCSignalV0_corrected_multip_pr;
   Double_t        track0tpcNsigma_corrected_pileup_el;
   Double_t        track0tpcNsigma_corrected_pileup_pi;
   Double_t        track0tpcNsigma_corrected_pileup_pr;
   Double_t        track1tpcNsigma_corrected_pileup_el;
   Double_t        track1tpcNsigma_corrected_pileup_pi;
   Double_t        track1tpcNsigma_corrected_pileup_pr;
   Double_t        track0ExpectedTPCSignalV0_corrected_pileup_el;
   Double_t        track0ExpectedTPCSignalV0_corrected_pileup_pi;
   Double_t        track0ExpectedTPCSignalV0_corrected_pileup_pr;
   Double_t        track1ExpectedTPCSignalV0_corrected_pileup_el;
   Double_t        track1ExpectedTPCSignalV0_corrected_pileup_pi;
   Double_t        track1ExpectedTPCSignalV0_corrected_pileup_pr;
   Double_t        track0tpcNsigma_corrected_eta_multip_el;
   Double_t        track0tpcNsigma_corrected_eta_multip_pi;
   Double_t        track0tpcNsigma_corrected_eta_multip_pr;
   Double_t        track1tpcNsigma_corrected_eta_multip_el;
   Double_t        track1tpcNsigma_corrected_eta_multip_pi;
   Double_t        track1tpcNsigma_corrected_eta_multip_pr;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_multip_el;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_multip_pi;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_multip_pr;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_multip_el;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_multip_pi;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_multip_pr;
   Double_t        track0tpcNsigma_corrected_eta_pileup_el;
   Double_t        track0tpcNsigma_corrected_eta_pileup_pi;
   Double_t        track0tpcNsigma_corrected_eta_pileup_pr;
   Double_t        track1tpcNsigma_corrected_eta_pileup_el;
   Double_t        track1tpcNsigma_corrected_eta_pileup_pi;
   Double_t        track1tpcNsigma_corrected_eta_pileup_pr;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_pileup_el;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_pileup_pi;
   Double_t        track0ExpectedTPCSignalV0_corrected_eta_pileup_pr;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_pileup_el;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_pileup_pi;
   Double_t        track1ExpectedTPCSignalV0_corrected_eta_pileup_pr;
   Double_t        logSignalTot_track0_0;
   Double_t        logSignalTot_track0_1;
   Double_t        logSignalTot_track0_2;
   Double_t        logSignalTot_track0_3;
   Double_t        logSignalMax_track0_0;
   Double_t        logSignalMax_track0_1;
   Double_t        logSignalMax_track0_2;
   Double_t        logSignalMax_track0_3;
   Double_t        signalNcl_track0_0;
   Double_t        signalNcl_track0_1;
   Double_t        signalNcl_track0_2;
   Double_t        signalNcl_track0_3;
   Double_t        signalNcr_track0_0;
   Double_t        signalNcr_track0_1;
   Double_t        signalNcr_track0_2;
   Double_t        signalNcr_track0_3;
   Double_t        logSignalTot_track1_0;
   Double_t        logSignalTot_track1_1;
   Double_t        logSignalTot_track1_2;
   Double_t        logSignalTot_track1_3;
   Double_t        logSignalMax_track1_0;
   Double_t        logSignalMax_track1_1;
   Double_t        logSignalMax_track1_2;
   Double_t        logSignalMax_track1_3;
   Double_t        signalNcl_track1_0;
   Double_t        signalNcl_track1_1;
   Double_t        signalNcl_track1_2;
   Double_t        signalNcl_track1_3;
   Double_t        signalNcr_track1_0;
   Double_t        signalNcr_track1_1;
   Double_t        signalNcr_track1_2;
   Double_t        signalNcr_track1_3;
   Double_t        track0P4;
   Double_t        track1P4;
   Double_t        track0CombP4;
   Double_t        track1CombP4;
   Double_t        track0C14;
   Double_t        track1C14;
   Double_t        track0CombC14;
   Double_t        track1CombC14;
   Double_t        dZ_Row0;
   Double_t        dRPhi_Row0;
   Double_t        dZ_Row62;
   Double_t        dRPhi_Row62;
   Double_t        dZ_Row125;
   Double_t        dRPhi_Row125;
   Double_t        dZ_Row158;
   Double_t        dRPhi_Row158;
   Double_t        v0_fRr;
   Double_t        phi0_ROC0;
   Double_t        phi1_ROC0;
   Double_t        phi0_ROC1;
   Double_t        phi1_ROC1;
   Double_t        phi0_ROC2;
   Double_t        phi1_ROC2;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_intrate;   //!
   TBranch        *b_timeStampS;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_bField;   //!
   TBranch        *b_selectionPtMask;   //!
   TBranch        *b_gid;   //!
   TBranch        *b_shiftA;   //!
   TBranch        *b_shiftC;   //!
   TBranch        *b_shiftM;   //!
   TBranch        *b_nPileUpPrim;   //!
   TBranch        *b_nPileUpSum;   //!
   TBranch        *b_nPileUpSumCorr;   //!
   TBranch        *b_primMult;   //!
   TBranch        *b_tpcClusterMult;   //!
   TBranch        *b_pileUpOK;   //!
   TBranch        *b_v0_fPointAngle;   //!
   TBranch        *b_kf_fChi2;   //!
   TBranch        *b_K0Like;   //!
   TBranch        *b_ELike;   //!
   TBranch        *b_LLike;   //!
   TBranch        *b_ALLike;   //!
   TBranch        *b_cleanK0;   //!
   TBranch        *b_cleanGamma;   //!
   TBranch        *b_cleanLambda;   //!
   TBranch        *b_cleanALambda;   //!
   TBranch        *b_K0Pull;   //!
   TBranch        *b_LPull;   //!
   TBranch        *b_ALPull;   //!
   TBranch        *b_EPull;   //!
   TBranch        *b_K0PullEff;   //!
   TBranch        *b_LPullEff;   //!
   TBranch        *b_ALPullEff;   //!
   TBranch        *b_EPullEff;   //!
   TBranch        *b_centV0;   //!
   TBranch        *b_centITS0;   //!
   TBranch        *b_centITS1;   //!
   TBranch        *b_multSSD;   //!
   TBranch        *b_multSDD;   //!
   TBranch        *b_tpcTrackBeforeClean;   //!
   TBranch        *b_triggerMask;   //!
   TBranch        *b_pairDCA_InerTPC;   //!
   TBranch        *b_isMinBias;   //!
   TBranch        *b_isPileUp;   //!
   TBranch        *b_track0_fTPCsignal;   //!
   TBranch        *b_track1_fTPCsignal;   //!
   TBranch        *b_track0_fTPCsignalN;   //!
   TBranch        *b_track1_fTPCsignalN;   //!
   TBranch        *b_type;   //!
   TBranch        *b_track0_fITSsignal;   //!
   TBranch        *b_track1_fITSsignal;   //!
   TBranch        *b_track0status;   //!
   TBranch        *b_track1status;   //!
   TBranch        *b_track0_hasTOF;   //!
   TBranch        *b_track1_hasTOF;   //!
   TBranch        *b_track0P;   //!
   TBranch        *b_track0Pt;   //!
   TBranch        *b_track0Eta;   //!
   TBranch        *b_track0Phi;   //!
   TBranch        *b_track0Px;   //!
   TBranch        *b_track0Py;   //!
   TBranch        *b_track0Pz;   //!
   TBranch        *b_track0Tgl;   //!
   TBranch        *b_track0DCAxy;   //!
   TBranch        *b_track0DCAz;   //!
   TBranch        *b_track0ncrossrow;   //!
   TBranch        *b_track0TPCchi2;   //!
   TBranch        *b_track0ITSchi2;   //!
   TBranch        *b_track0TPC_cls;   //!
   TBranch        *b_track0ITS_cls;   //!
   TBranch        *b_track1P;   //!
   TBranch        *b_track1Pt;   //!
   TBranch        *b_track1Eta;   //!
   TBranch        *b_track1Phi;   //!
   TBranch        *b_track1Px;   //!
   TBranch        *b_track1Py;   //!
   TBranch        *b_track1Pz;   //!
   TBranch        *b_track1Tgl;   //!
   TBranch        *b_track1DCAxy;   //!
   TBranch        *b_track1DCAz;   //!
   TBranch        *b_track1ncrossrow;   //!
   TBranch        *b_track1TPCchi2;   //!
   TBranch        *b_track1ITSchi2;   //!
   TBranch        *b_track1TPC_cls;   //!
   TBranch        *b_track1ITS_cls;   //!
   TBranch        *b_track0DCAxy_TPC;   //!
   TBranch        *b_track0DCAz_TPC;   //!
   TBranch        *b_fAlpha0;   //!
   TBranch        *b_fITSClusterMap0;   //!
   TBranch        *b_fSigned1Pt0;   //!
   TBranch        *b_track1DCAxy_TPC;   //!
   TBranch        *b_track1DCAz_TPC;   //!
   TBranch        *b_fAlpha1;   //!
   TBranch        *b_fITSClusterMap1;   //!
   TBranch        *b_fSigned1Pt1;   //!
   TBranch        *b_track0tofNsigmaElectron;   //!
   TBranch        *b_track0tofNsigmaPion;   //!
   TBranch        *b_track0tofNsigmaKaon;   //!
   TBranch        *b_track0tofNsigmaProton;   //!
   TBranch        *b_track1tofNsigmaElectron;   //!
   TBranch        *b_track1tofNsigmaPion;   //!
   TBranch        *b_track1tofNsigmaKaon;   //!
   TBranch        *b_track1tofNsigmaProton;   //!
   TBranch        *b_track0tpcNsigma_el;   //!
   TBranch        *b_track0tpcNsigma_pi;   //!
   TBranch        *b_track0tpcNsigma_ka;   //!
   TBranch        *b_track0tpcNsigma_pro;   //!
   TBranch        *b_track1tpcNsigma_el;   //!
   TBranch        *b_track1tpcNsigma_pi;   //!
   TBranch        *b_track1tpcNsigma_ka;   //!
   TBranch        *b_track1tpcNsigma_pro;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_ka;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_pro;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_ka;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_pro;   //!
   TBranch        *b_track0CorrectedTPCSignalV0_el;   //!
   TBranch        *b_track0CorrectedTPCSignalV0_pi;   //!
   TBranch        *b_track0CorrectedTPCSignalV0_ka;   //!
   TBranch        *b_track0CorrectedTPCSignalV0_pro;   //!
   TBranch        *b_track1CorrectedTPCSignalV0_el;   //!
   TBranch        *b_track1CorrectedTPCSignalV0_pi;   //!
   TBranch        *b_track1CorrectedTPCSignalV0_ka;   //!
   TBranch        *b_track1CorrectedTPCSignalV0_pro;   //!
   TBranch        *b_track0GetPileupValue;   //!
   TBranch        *b_track1GetPileupValue;   //!
   TBranch        *b_track0tpcNsigma_no_corrected_el;   //!
   TBranch        *b_track0tpcNsigma_no_corrected_pi;   //!
   TBranch        *b_track0tpcNsigma_no_corrected_pr;   //!
   TBranch        *b_track1tpcNsigma_no_corrected_el;   //!
   TBranch        *b_track1tpcNsigma_no_corrected_pi;   //!
   TBranch        *b_track1tpcNsigma_no_corrected_pr;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_no_corrected_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_no_corrected_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_no_corrected_pr;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_no_corrected_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_no_corrected_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_no_corrected_pr;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_el;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_pi;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_pr;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_el;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_pi;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_pr;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_pr;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_pr;   //!
   TBranch        *b_track0tpcNsigma_corrected_multip_el;   //!
   TBranch        *b_track0tpcNsigma_corrected_multip_pi;   //!
   TBranch        *b_track0tpcNsigma_corrected_multip_pr;   //!
   TBranch        *b_track1tpcNsigma_corrected_multip_el;   //!
   TBranch        *b_track1tpcNsigma_corrected_multip_pi;   //!
   TBranch        *b_track1tpcNsigma_corrected_multip_pr;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_multip_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_multip_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_multip_pr;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_multip_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_multip_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_multip_pr;   //!
   TBranch        *b_track0tpcNsigma_corrected_pileup_el;   //!
   TBranch        *b_track0tpcNsigma_corrected_pileup_pi;   //!
   TBranch        *b_track0tpcNsigma_corrected_pileup_pr;   //!
   TBranch        *b_track1tpcNsigma_corrected_pileup_el;   //!
   TBranch        *b_track1tpcNsigma_corrected_pileup_pi;   //!
   TBranch        *b_track1tpcNsigma_corrected_pileup_pr;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_pileup_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_pileup_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_pileup_pr;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_pileup_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_pileup_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_pileup_pr;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_multip_el;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_multip_pi;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_multip_pr;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_multip_el;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_multip_pi;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_multip_pr;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_multip_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_multip_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_multip_pr;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_multip_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_multip_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_multip_pr;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_pileup_el;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_pileup_pi;   //!
   TBranch        *b_track0tpcNsigma_corrected_eta_pileup_pr;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_pileup_el;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_pileup_pi;   //!
   TBranch        *b_track1tpcNsigma_corrected_eta_pileup_pr;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_pileup_el;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_pileup_pi;   //!
   TBranch        *b_track0ExpectedTPCSignalV0_corrected_eta_pileup_pr;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_pileup_el;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_pileup_pi;   //!
   TBranch        *b_track1ExpectedTPCSignalV0_corrected_eta_pileup_pr;   //!
   TBranch        *b_logSignalTot_track0_0;   //!
   TBranch        *b_logSignalTot_track0_1;   //!
   TBranch        *b_logSignalTot_track0_2;   //!
   TBranch        *b_logSignalTot_track0_3;   //!
   TBranch        *b_logSignalMax_track0_0;   //!
   TBranch        *b_logSignalMax_track0_1;   //!
   TBranch        *b_logSignalMax_track0_2;   //!
   TBranch        *b_logSignalMax_track0_3;   //!
   TBranch        *b_signalNcl_track0_0;   //!
   TBranch        *b_signalNcl_track0_1;   //!
   TBranch        *b_signalNcl_track0_2;   //!
   TBranch        *b_signalNcl_track0_3;   //!
   TBranch        *b_signalNcr_track0_0;   //!
   TBranch        *b_signalNcr_track0_1;   //!
   TBranch        *b_signalNcr_track0_2;   //!
   TBranch        *b_signalNcr_track0_3;   //!
   TBranch        *b_logSignalTot_track1_0;   //!
   TBranch        *b_logSignalTot_track1_1;   //!
   TBranch        *b_logSignalTot_track1_2;   //!
   TBranch        *b_logSignalTot_track1_3;   //!
   TBranch        *b_logSignalMax_track1_0;   //!
   TBranch        *b_logSignalMax_track1_1;   //!
   TBranch        *b_logSignalMax_track1_2;   //!
   TBranch        *b_logSignalMax_track1_3;   //!
   TBranch        *b_signalNcl_track1_0;   //!
   TBranch        *b_signalNcl_track1_1;   //!
   TBranch        *b_signalNcl_track1_2;   //!
   TBranch        *b_signalNcl_track1_3;   //!
   TBranch        *b_signalNcr_track1_0;   //!
   TBranch        *b_signalNcr_track1_1;   //!
   TBranch        *b_signalNcr_track1_2;   //!
   TBranch        *b_signalNcr_track1_3;   //!
   TBranch        *b_track0P4;   //!
   TBranch        *b_track1P4;   //!
   TBranch        *b_track0CombP4;   //!
   TBranch        *b_track1CombP4;   //!
   TBranch        *b_track0C14;   //!
   TBranch        *b_track1C14;   //!
   TBranch        *b_track0CombC14;   //!
   TBranch        *b_track1CombC14;   //!
   TBranch        *b_dZ_Row0;   //!
   TBranch        *b_dRPhi_Row0;   //!
   TBranch        *b_dZ_Row62;   //!
   TBranch        *b_dRPhi_Row62;   //!
   TBranch        *b_dZ_Row125;   //!
   TBranch        *b_dRPhi_Row125;   //!
   TBranch        *b_dZ_Row158;   //!
   TBranch        *b_dRPhi_Row158;   //!
   TBranch        *b_v0_fRr;   //!
   TBranch        *b_phi0_ROC0;   //!
   TBranch        *b_phi1_ROC0;   //!
   TBranch        *b_phi0_ROC1;   //!
   TBranch        *b_phi1_ROC1;   //!
   TBranch        *b_phi0_ROC2;   //!
   TBranch        *b_phi1_ROC2;   //!

   V0FlatAna(TTree *tree=0);
   virtual ~V0FlatAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef V0FlatAna_cxx
V0FlatAna::V0FlatAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("V0tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("V0tree.root");
      }
      f->GetObject("V0Flat",tree);

   }
   Init(tree);
}

V0FlatAna::~V0FlatAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t V0FlatAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t V0FlatAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void V0FlatAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("intrate", &intrate, &b_intrate);
   fChain->SetBranchAddress("timeStampS", &timeStampS, &b_timeStampS);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("bField", &bField, &b_bField);
   fChain->SetBranchAddress("selectionPtMask", &selectionPtMask, &b_selectionPtMask);
   fChain->SetBranchAddress("gid", &gid, &b_gid);
   fChain->SetBranchAddress("shiftA", &shiftA, &b_shiftA);
   fChain->SetBranchAddress("shiftC", &shiftC, &b_shiftC);
   fChain->SetBranchAddress("shiftM", &shiftM, &b_shiftM);
   fChain->SetBranchAddress("nPileUpPrim", &nPileUpPrim, &b_nPileUpPrim);
   fChain->SetBranchAddress("nPileUpSum", &nPileUpSum, &b_nPileUpSum);
   fChain->SetBranchAddress("nPileUpSumCorr", &nPileUpSumCorr, &b_nPileUpSumCorr);
   fChain->SetBranchAddress("primMult", &primMult, &b_primMult);
   fChain->SetBranchAddress("tpcClusterMult", &tpcClusterMult, &b_tpcClusterMult);
   fChain->SetBranchAddress("pileUpOK", &pileUpOK, &b_pileUpOK);
   fChain->SetBranchAddress("v0.fPointAngle", &v0_fPointAngle, &b_v0_fPointAngle);
   fChain->SetBranchAddress("kf.fChi2", &kf_fChi2, &b_kf_fChi2);
   fChain->SetBranchAddress("K0Like", &K0Like, &b_K0Like);
   fChain->SetBranchAddress("ELike", &ELike, &b_ELike);
   fChain->SetBranchAddress("LLike", &LLike, &b_LLike);
   fChain->SetBranchAddress("ALLike", &ALLike, &b_ALLike);
   fChain->SetBranchAddress("cleanK0", &cleanK0, &b_cleanK0);
   fChain->SetBranchAddress("cleanGamma", &cleanGamma, &b_cleanGamma);
   fChain->SetBranchAddress("cleanLambda", &cleanLambda, &b_cleanLambda);
   fChain->SetBranchAddress("cleanALambda", &cleanALambda, &b_cleanALambda);
   fChain->SetBranchAddress("K0Pull", &K0Pull, &b_K0Pull);
   fChain->SetBranchAddress("LPull", &LPull, &b_LPull);
   fChain->SetBranchAddress("ALPull", &ALPull, &b_ALPull);
   fChain->SetBranchAddress("EPull", &EPull, &b_EPull);
   fChain->SetBranchAddress("K0PullEff", &K0PullEff, &b_K0PullEff);
   fChain->SetBranchAddress("LPullEff", &LPullEff, &b_LPullEff);
   fChain->SetBranchAddress("ALPullEff", &ALPullEff, &b_ALPullEff);
   fChain->SetBranchAddress("EPullEff", &EPullEff, &b_EPullEff);
   fChain->SetBranchAddress("centV0", &centV0, &b_centV0);
   fChain->SetBranchAddress("centITS0", &centITS0, &b_centITS0);
   fChain->SetBranchAddress("centITS1", &centITS1, &b_centITS1);
   fChain->SetBranchAddress("multSSD", &multSSD, &b_multSSD);
   fChain->SetBranchAddress("multSDD", &multSDD, &b_multSDD);
   fChain->SetBranchAddress("tpcTrackBeforeClean", &tpcTrackBeforeClean, &b_tpcTrackBeforeClean);
   fChain->SetBranchAddress("triggerMask", &triggerMask, &b_triggerMask);
   fChain->SetBranchAddress("pairDCA_InerTPC", &pairDCA_InerTPC, &b_pairDCA_InerTPC);
   fChain->SetBranchAddress("isMinBias", &isMinBias, &b_isMinBias);
   fChain->SetBranchAddress("isPileUp", &isPileUp, &b_isPileUp);
   fChain->SetBranchAddress("track0.fTPCsignal", &track0_fTPCsignal, &b_track0_fTPCsignal);
   fChain->SetBranchAddress("track1.fTPCsignal", &track1_fTPCsignal, &b_track1_fTPCsignal);
   fChain->SetBranchAddress("track0.fTPCsignalN", &track0_fTPCsignalN, &b_track0_fTPCsignalN);
   fChain->SetBranchAddress("track1.fTPCsignalN", &track1_fTPCsignalN, &b_track1_fTPCsignalN);
   fChain->SetBranchAddress("type", &type, &b_type);
   fChain->SetBranchAddress("track0.fITSsignal", &track0_fITSsignal, &b_track0_fITSsignal);
   fChain->SetBranchAddress("track1.fITSsignal", &track1_fITSsignal, &b_track1_fITSsignal);
   fChain->SetBranchAddress("track0status", &track0status, &b_track0status);
   fChain->SetBranchAddress("track1status", &track1status, &b_track1status);
   fChain->SetBranchAddress("track0_hasTOF", &track0_hasTOF, &b_track0_hasTOF);
   fChain->SetBranchAddress("track1_hasTOF", &track1_hasTOF, &b_track1_hasTOF);
   fChain->SetBranchAddress("track0P", &track0P, &b_track0P);
   fChain->SetBranchAddress("track0Pt", &track0Pt, &b_track0Pt);
   fChain->SetBranchAddress("track0Eta", &track0Eta, &b_track0Eta);
   fChain->SetBranchAddress("track0Phi", &track0Phi, &b_track0Phi);
   fChain->SetBranchAddress("track0Px", &track0Px, &b_track0Px);
   fChain->SetBranchAddress("track0Py", &track0Py, &b_track0Py);
   fChain->SetBranchAddress("track0Pz", &track0Pz, &b_track0Pz);
   fChain->SetBranchAddress("track0Tgl", &track0Tgl, &b_track0Tgl);
   fChain->SetBranchAddress("track0DCAxy", &track0DCAxy, &b_track0DCAxy);
   fChain->SetBranchAddress("track0DCAz", &track0DCAz, &b_track0DCAz);
   fChain->SetBranchAddress("track0ncrossrow", &track0ncrossrow, &b_track0ncrossrow);
   fChain->SetBranchAddress("track0TPCchi2", &track0TPCchi2, &b_track0TPCchi2);
   fChain->SetBranchAddress("track0ITSchi2", &track0ITSchi2, &b_track0ITSchi2);
   fChain->SetBranchAddress("track0TPC_cls", &track0TPC_cls, &b_track0TPC_cls);
   fChain->SetBranchAddress("track0ITS_cls", &track0ITS_cls, &b_track0ITS_cls);
   fChain->SetBranchAddress("track1P", &track1P, &b_track1P);
   fChain->SetBranchAddress("track1Pt", &track1Pt, &b_track1Pt);
   fChain->SetBranchAddress("track1Eta", &track1Eta, &b_track1Eta);
   fChain->SetBranchAddress("track1Phi", &track1Phi, &b_track1Phi);
   fChain->SetBranchAddress("track1Px", &track1Px, &b_track1Px);
   fChain->SetBranchAddress("track1Py", &track1Py, &b_track1Py);
   fChain->SetBranchAddress("track1Pz", &track1Pz, &b_track1Pz);
   fChain->SetBranchAddress("track1Tgl", &track1Tgl, &b_track1Tgl);
   fChain->SetBranchAddress("track1DCAxy", &track1DCAxy, &b_track1DCAxy);
   fChain->SetBranchAddress("track1DCAz", &track1DCAz, &b_track1DCAz);
   fChain->SetBranchAddress("track1ncrossrow", &track1ncrossrow, &b_track1ncrossrow);
   fChain->SetBranchAddress("track1TPCchi2", &track1TPCchi2, &b_track1TPCchi2);
   fChain->SetBranchAddress("track1ITSchi2", &track1ITSchi2, &b_track1ITSchi2);
   fChain->SetBranchAddress("track1TPC_cls", &track1TPC_cls, &b_track1TPC_cls);
   fChain->SetBranchAddress("track1ITS_cls", &track1ITS_cls, &b_track1ITS_cls);
   fChain->SetBranchAddress("track0DCAxy_TPC", &track0DCAxy_TPC, &b_track0DCAxy_TPC);
   fChain->SetBranchAddress("track0DCAz_TPC", &track0DCAz_TPC, &b_track0DCAz_TPC);
   fChain->SetBranchAddress("fAlpha0", &fAlpha0, &b_fAlpha0);
   fChain->SetBranchAddress("fITSClusterMap0", &fITSClusterMap0, &b_fITSClusterMap0);
   fChain->SetBranchAddress("fSigned1Pt0", &fSigned1Pt0, &b_fSigned1Pt0);
   fChain->SetBranchAddress("track1DCAxy_TPC", &track1DCAxy_TPC, &b_track1DCAxy_TPC);
   fChain->SetBranchAddress("track1DCAz_TPC", &track1DCAz_TPC, &b_track1DCAz_TPC);
   fChain->SetBranchAddress("fAlpha1", &fAlpha1, &b_fAlpha1);
   fChain->SetBranchAddress("fITSClusterMap1", &fITSClusterMap1, &b_fITSClusterMap1);
   fChain->SetBranchAddress("fSigned1Pt1", &fSigned1Pt1, &b_fSigned1Pt1);
   fChain->SetBranchAddress("track0tofNsigmaElectron", &track0tofNsigmaElectron, &b_track0tofNsigmaElectron);
   fChain->SetBranchAddress("track0tofNsigmaPion", &track0tofNsigmaPion, &b_track0tofNsigmaPion);
   fChain->SetBranchAddress("track0tofNsigmaKaon", &track0tofNsigmaKaon, &b_track0tofNsigmaKaon);
   fChain->SetBranchAddress("track0tofNsigmaProton", &track0tofNsigmaProton, &b_track0tofNsigmaProton);
   fChain->SetBranchAddress("track1tofNsigmaElectron", &track1tofNsigmaElectron, &b_track1tofNsigmaElectron);
   fChain->SetBranchAddress("track1tofNsigmaPion", &track1tofNsigmaPion, &b_track1tofNsigmaPion);
   fChain->SetBranchAddress("track1tofNsigmaKaon", &track1tofNsigmaKaon, &b_track1tofNsigmaKaon);
   fChain->SetBranchAddress("track1tofNsigmaProton", &track1tofNsigmaProton, &b_track1tofNsigmaProton);
   fChain->SetBranchAddress("track0tpcNsigma_el", &track0tpcNsigma_el, &b_track0tpcNsigma_el);
   fChain->SetBranchAddress("track0tpcNsigma_pi", &track0tpcNsigma_pi, &b_track0tpcNsigma_pi);
   fChain->SetBranchAddress("track0tpcNsigma_ka", &track0tpcNsigma_ka, &b_track0tpcNsigma_ka);
   fChain->SetBranchAddress("track0tpcNsigma_pro", &track0tpcNsigma_pro, &b_track0tpcNsigma_pro);
   fChain->SetBranchAddress("track1tpcNsigma_el", &track1tpcNsigma_el, &b_track1tpcNsigma_el);
   fChain->SetBranchAddress("track1tpcNsigma_pi", &track1tpcNsigma_pi, &b_track1tpcNsigma_pi);
   fChain->SetBranchAddress("track1tpcNsigma_ka", &track1tpcNsigma_ka, &b_track1tpcNsigma_ka);
   fChain->SetBranchAddress("track1tpcNsigma_pro", &track1tpcNsigma_pro, &b_track1tpcNsigma_pro);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_el", &track0ExpectedTPCSignalV0_el, &b_track0ExpectedTPCSignalV0_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_pi", &track0ExpectedTPCSignalV0_pi, &b_track0ExpectedTPCSignalV0_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_ka", &track0ExpectedTPCSignalV0_ka, &b_track0ExpectedTPCSignalV0_ka);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_pro", &track0ExpectedTPCSignalV0_pro, &b_track0ExpectedTPCSignalV0_pro);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_el", &track1ExpectedTPCSignalV0_el, &b_track1ExpectedTPCSignalV0_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_pi", &track1ExpectedTPCSignalV0_pi, &b_track1ExpectedTPCSignalV0_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_ka", &track1ExpectedTPCSignalV0_ka, &b_track1ExpectedTPCSignalV0_ka);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_pro", &track1ExpectedTPCSignalV0_pro, &b_track1ExpectedTPCSignalV0_pro);
   fChain->SetBranchAddress("track0CorrectedTPCSignalV0_el", &track0CorrectedTPCSignalV0_el, &b_track0CorrectedTPCSignalV0_el);
   fChain->SetBranchAddress("track0CorrectedTPCSignalV0_pi", &track0CorrectedTPCSignalV0_pi, &b_track0CorrectedTPCSignalV0_pi);
   fChain->SetBranchAddress("track0CorrectedTPCSignalV0_ka", &track0CorrectedTPCSignalV0_ka, &b_track0CorrectedTPCSignalV0_ka);
   fChain->SetBranchAddress("track0CorrectedTPCSignalV0_pro", &track0CorrectedTPCSignalV0_pro, &b_track0CorrectedTPCSignalV0_pro);
   fChain->SetBranchAddress("track1CorrectedTPCSignalV0_el", &track1CorrectedTPCSignalV0_el, &b_track1CorrectedTPCSignalV0_el);
   fChain->SetBranchAddress("track1CorrectedTPCSignalV0_pi", &track1CorrectedTPCSignalV0_pi, &b_track1CorrectedTPCSignalV0_pi);
   fChain->SetBranchAddress("track1CorrectedTPCSignalV0_ka", &track1CorrectedTPCSignalV0_ka, &b_track1CorrectedTPCSignalV0_ka);
   fChain->SetBranchAddress("track1CorrectedTPCSignalV0_pro", &track1CorrectedTPCSignalV0_pro, &b_track1CorrectedTPCSignalV0_pro);
   fChain->SetBranchAddress("track0GetPileupValue", &track0GetPileupValue, &b_track0GetPileupValue);
   fChain->SetBranchAddress("track1GetPileupValue", &track1GetPileupValue, &b_track1GetPileupValue);
   fChain->SetBranchAddress("track0tpcNsigma_no_corrected_el", &track0tpcNsigma_no_corrected_el, &b_track0tpcNsigma_no_corrected_el);
   fChain->SetBranchAddress("track0tpcNsigma_no_corrected_pi", &track0tpcNsigma_no_corrected_pi, &b_track0tpcNsigma_no_corrected_pi);
   fChain->SetBranchAddress("track0tpcNsigma_no_corrected_pr", &track0tpcNsigma_no_corrected_pr, &b_track0tpcNsigma_no_corrected_pr);
   fChain->SetBranchAddress("track1tpcNsigma_no_corrected_el", &track1tpcNsigma_no_corrected_el, &b_track1tpcNsigma_no_corrected_el);
   fChain->SetBranchAddress("track1tpcNsigma_no_corrected_pi", &track1tpcNsigma_no_corrected_pi, &b_track1tpcNsigma_no_corrected_pi);
   fChain->SetBranchAddress("track1tpcNsigma_no_corrected_pr", &track1tpcNsigma_no_corrected_pr, &b_track1tpcNsigma_no_corrected_pr);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_no_corrected_el", &track0ExpectedTPCSignalV0_no_corrected_el, &b_track0ExpectedTPCSignalV0_no_corrected_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_no_corrected_pi", &track0ExpectedTPCSignalV0_no_corrected_pi, &b_track0ExpectedTPCSignalV0_no_corrected_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_no_corrected_pr", &track0ExpectedTPCSignalV0_no_corrected_pr, &b_track0ExpectedTPCSignalV0_no_corrected_pr);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_no_corrected_el", &track1ExpectedTPCSignalV0_no_corrected_el, &b_track1ExpectedTPCSignalV0_no_corrected_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_no_corrected_pi", &track1ExpectedTPCSignalV0_no_corrected_pi, &b_track1ExpectedTPCSignalV0_no_corrected_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_no_corrected_pr", &track1ExpectedTPCSignalV0_no_corrected_pr, &b_track1ExpectedTPCSignalV0_no_corrected_pr);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_el", &track0tpcNsigma_corrected_eta_el, &b_track0tpcNsigma_corrected_eta_el);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_pi", &track0tpcNsigma_corrected_eta_pi, &b_track0tpcNsigma_corrected_eta_pi);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_pr", &track0tpcNsigma_corrected_eta_pr, &b_track0tpcNsigma_corrected_eta_pr);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_el", &track1tpcNsigma_corrected_eta_el, &b_track1tpcNsigma_corrected_eta_el);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_pi", &track1tpcNsigma_corrected_eta_pi, &b_track1tpcNsigma_corrected_eta_pi);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_pr", &track1tpcNsigma_corrected_eta_pr, &b_track1tpcNsigma_corrected_eta_pr);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_el", &track0ExpectedTPCSignalV0_corrected_eta_el, &b_track0ExpectedTPCSignalV0_corrected_eta_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_pi", &track0ExpectedTPCSignalV0_corrected_eta_pi, &b_track0ExpectedTPCSignalV0_corrected_eta_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_pr", &track0ExpectedTPCSignalV0_corrected_eta_pr, &b_track0ExpectedTPCSignalV0_corrected_eta_pr);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_el", &track1ExpectedTPCSignalV0_corrected_eta_el, &b_track1ExpectedTPCSignalV0_corrected_eta_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_pi", &track1ExpectedTPCSignalV0_corrected_eta_pi, &b_track1ExpectedTPCSignalV0_corrected_eta_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_pr", &track1ExpectedTPCSignalV0_corrected_eta_pr, &b_track1ExpectedTPCSignalV0_corrected_eta_pr);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_multip_el", &track0tpcNsigma_corrected_multip_el, &b_track0tpcNsigma_corrected_multip_el);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_multip_pi", &track0tpcNsigma_corrected_multip_pi, &b_track0tpcNsigma_corrected_multip_pi);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_multip_pr", &track0tpcNsigma_corrected_multip_pr, &b_track0tpcNsigma_corrected_multip_pr);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_multip_el", &track1tpcNsigma_corrected_multip_el, &b_track1tpcNsigma_corrected_multip_el);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_multip_pi", &track1tpcNsigma_corrected_multip_pi, &b_track1tpcNsigma_corrected_multip_pi);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_multip_pr", &track1tpcNsigma_corrected_multip_pr, &b_track1tpcNsigma_corrected_multip_pr);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_multip_el", &track0ExpectedTPCSignalV0_corrected_multip_el, &b_track0ExpectedTPCSignalV0_corrected_multip_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_multip_pi", &track0ExpectedTPCSignalV0_corrected_multip_pi, &b_track0ExpectedTPCSignalV0_corrected_multip_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_multip_pr", &track0ExpectedTPCSignalV0_corrected_multip_pr, &b_track0ExpectedTPCSignalV0_corrected_multip_pr);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_multip_el", &track1ExpectedTPCSignalV0_corrected_multip_el, &b_track1ExpectedTPCSignalV0_corrected_multip_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_multip_pi", &track1ExpectedTPCSignalV0_corrected_multip_pi, &b_track1ExpectedTPCSignalV0_corrected_multip_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_multip_pr", &track1ExpectedTPCSignalV0_corrected_multip_pr, &b_track1ExpectedTPCSignalV0_corrected_multip_pr);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_pileup_el", &track0tpcNsigma_corrected_pileup_el, &b_track0tpcNsigma_corrected_pileup_el);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_pileup_pi", &track0tpcNsigma_corrected_pileup_pi, &b_track0tpcNsigma_corrected_pileup_pi);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_pileup_pr", &track0tpcNsigma_corrected_pileup_pr, &b_track0tpcNsigma_corrected_pileup_pr);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_pileup_el", &track1tpcNsigma_corrected_pileup_el, &b_track1tpcNsigma_corrected_pileup_el);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_pileup_pi", &track1tpcNsigma_corrected_pileup_pi, &b_track1tpcNsigma_corrected_pileup_pi);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_pileup_pr", &track1tpcNsigma_corrected_pileup_pr, &b_track1tpcNsigma_corrected_pileup_pr);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_pileup_el", &track0ExpectedTPCSignalV0_corrected_pileup_el, &b_track0ExpectedTPCSignalV0_corrected_pileup_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_pileup_pi", &track0ExpectedTPCSignalV0_corrected_pileup_pi, &b_track0ExpectedTPCSignalV0_corrected_pileup_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_pileup_pr", &track0ExpectedTPCSignalV0_corrected_pileup_pr, &b_track0ExpectedTPCSignalV0_corrected_pileup_pr);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_pileup_el", &track1ExpectedTPCSignalV0_corrected_pileup_el, &b_track1ExpectedTPCSignalV0_corrected_pileup_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_pileup_pi", &track1ExpectedTPCSignalV0_corrected_pileup_pi, &b_track1ExpectedTPCSignalV0_corrected_pileup_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_pileup_pr", &track1ExpectedTPCSignalV0_corrected_pileup_pr, &b_track1ExpectedTPCSignalV0_corrected_pileup_pr);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_multip_el", &track0tpcNsigma_corrected_eta_multip_el, &b_track0tpcNsigma_corrected_eta_multip_el);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_multip_pi", &track0tpcNsigma_corrected_eta_multip_pi, &b_track0tpcNsigma_corrected_eta_multip_pi);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_multip_pr", &track0tpcNsigma_corrected_eta_multip_pr, &b_track0tpcNsigma_corrected_eta_multip_pr);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_multip_el", &track1tpcNsigma_corrected_eta_multip_el, &b_track1tpcNsigma_corrected_eta_multip_el);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_multip_pi", &track1tpcNsigma_corrected_eta_multip_pi, &b_track1tpcNsigma_corrected_eta_multip_pi);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_multip_pr", &track1tpcNsigma_corrected_eta_multip_pr, &b_track1tpcNsigma_corrected_eta_multip_pr);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_multip_el", &track0ExpectedTPCSignalV0_corrected_eta_multip_el, &b_track0ExpectedTPCSignalV0_corrected_eta_multip_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_multip_pi", &track0ExpectedTPCSignalV0_corrected_eta_multip_pi, &b_track0ExpectedTPCSignalV0_corrected_eta_multip_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_multip_pr", &track0ExpectedTPCSignalV0_corrected_eta_multip_pr, &b_track0ExpectedTPCSignalV0_corrected_eta_multip_pr);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_multip_el", &track1ExpectedTPCSignalV0_corrected_eta_multip_el, &b_track1ExpectedTPCSignalV0_corrected_eta_multip_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_multip_pi", &track1ExpectedTPCSignalV0_corrected_eta_multip_pi, &b_track1ExpectedTPCSignalV0_corrected_eta_multip_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_multip_pr", &track1ExpectedTPCSignalV0_corrected_eta_multip_pr, &b_track1ExpectedTPCSignalV0_corrected_eta_multip_pr);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_pileup_el", &track0tpcNsigma_corrected_eta_pileup_el, &b_track0tpcNsigma_corrected_eta_pileup_el);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_pileup_pi", &track0tpcNsigma_corrected_eta_pileup_pi, &b_track0tpcNsigma_corrected_eta_pileup_pi);
   fChain->SetBranchAddress("track0tpcNsigma_corrected_eta_pileup_pr", &track0tpcNsigma_corrected_eta_pileup_pr, &b_track0tpcNsigma_corrected_eta_pileup_pr);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_pileup_el", &track1tpcNsigma_corrected_eta_pileup_el, &b_track1tpcNsigma_corrected_eta_pileup_el);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_pileup_pi", &track1tpcNsigma_corrected_eta_pileup_pi, &b_track1tpcNsigma_corrected_eta_pileup_pi);
   fChain->SetBranchAddress("track1tpcNsigma_corrected_eta_pileup_pr", &track1tpcNsigma_corrected_eta_pileup_pr, &b_track1tpcNsigma_corrected_eta_pileup_pr);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_pileup_el", &track0ExpectedTPCSignalV0_corrected_eta_pileup_el, &b_track0ExpectedTPCSignalV0_corrected_eta_pileup_el);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_pileup_pi", &track0ExpectedTPCSignalV0_corrected_eta_pileup_pi, &b_track0ExpectedTPCSignalV0_corrected_eta_pileup_pi);
   fChain->SetBranchAddress("track0ExpectedTPCSignalV0_corrected_eta_pileup_pr", &track0ExpectedTPCSignalV0_corrected_eta_pileup_pr, &b_track0ExpectedTPCSignalV0_corrected_eta_pileup_pr);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_pileup_el", &track1ExpectedTPCSignalV0_corrected_eta_pileup_el, &b_track1ExpectedTPCSignalV0_corrected_eta_pileup_el);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_pileup_pi", &track1ExpectedTPCSignalV0_corrected_eta_pileup_pi, &b_track1ExpectedTPCSignalV0_corrected_eta_pileup_pi);
   fChain->SetBranchAddress("track1ExpectedTPCSignalV0_corrected_eta_pileup_pr", &track1ExpectedTPCSignalV0_corrected_eta_pileup_pr, &b_track1ExpectedTPCSignalV0_corrected_eta_pileup_pr);
   fChain->SetBranchAddress("logSignalTot_track0_0", &logSignalTot_track0_0, &b_logSignalTot_track0_0);
   fChain->SetBranchAddress("logSignalTot_track0_1", &logSignalTot_track0_1, &b_logSignalTot_track0_1);
   fChain->SetBranchAddress("logSignalTot_track0_2", &logSignalTot_track0_2, &b_logSignalTot_track0_2);
   fChain->SetBranchAddress("logSignalTot_track0_3", &logSignalTot_track0_3, &b_logSignalTot_track0_3);
   fChain->SetBranchAddress("logSignalMax_track0_0", &logSignalMax_track0_0, &b_logSignalMax_track0_0);
   fChain->SetBranchAddress("logSignalMax_track0_1", &logSignalMax_track0_1, &b_logSignalMax_track0_1);
   fChain->SetBranchAddress("logSignalMax_track0_2", &logSignalMax_track0_2, &b_logSignalMax_track0_2);
   fChain->SetBranchAddress("logSignalMax_track0_3", &logSignalMax_track0_3, &b_logSignalMax_track0_3);
   fChain->SetBranchAddress("signalNcl_track0_0", &signalNcl_track0_0, &b_signalNcl_track0_0);
   fChain->SetBranchAddress("signalNcl_track0_1", &signalNcl_track0_1, &b_signalNcl_track0_1);
   fChain->SetBranchAddress("signalNcl_track0_2", &signalNcl_track0_2, &b_signalNcl_track0_2);
   fChain->SetBranchAddress("signalNcl_track0_3", &signalNcl_track0_3, &b_signalNcl_track0_3);
   fChain->SetBranchAddress("signalNcr_track0_0", &signalNcr_track0_0, &b_signalNcr_track0_0);
   fChain->SetBranchAddress("signalNcr_track0_1", &signalNcr_track0_1, &b_signalNcr_track0_1);
   fChain->SetBranchAddress("signalNcr_track0_2", &signalNcr_track0_2, &b_signalNcr_track0_2);
   fChain->SetBranchAddress("signalNcr_track0_3", &signalNcr_track0_3, &b_signalNcr_track0_3);
   fChain->SetBranchAddress("logSignalTot_track1_0", &logSignalTot_track1_0, &b_logSignalTot_track1_0);
   fChain->SetBranchAddress("logSignalTot_track1_1", &logSignalTot_track1_1, &b_logSignalTot_track1_1);
   fChain->SetBranchAddress("logSignalTot_track1_2", &logSignalTot_track1_2, &b_logSignalTot_track1_2);
   fChain->SetBranchAddress("logSignalTot_track1_3", &logSignalTot_track1_3, &b_logSignalTot_track1_3);
   fChain->SetBranchAddress("logSignalMax_track1_0", &logSignalMax_track1_0, &b_logSignalMax_track1_0);
   fChain->SetBranchAddress("logSignalMax_track1_1", &logSignalMax_track1_1, &b_logSignalMax_track1_1);
   fChain->SetBranchAddress("logSignalMax_track1_2", &logSignalMax_track1_2, &b_logSignalMax_track1_2);
   fChain->SetBranchAddress("logSignalMax_track1_3", &logSignalMax_track1_3, &b_logSignalMax_track1_3);
   fChain->SetBranchAddress("signalNcl_track1_0", &signalNcl_track1_0, &b_signalNcl_track1_0);
   fChain->SetBranchAddress("signalNcl_track1_1", &signalNcl_track1_1, &b_signalNcl_track1_1);
   fChain->SetBranchAddress("signalNcl_track1_2", &signalNcl_track1_2, &b_signalNcl_track1_2);
   fChain->SetBranchAddress("signalNcl_track1_3", &signalNcl_track1_3, &b_signalNcl_track1_3);
   fChain->SetBranchAddress("signalNcr_track1_0", &signalNcr_track1_0, &b_signalNcr_track1_0);
   fChain->SetBranchAddress("signalNcr_track1_1", &signalNcr_track1_1, &b_signalNcr_track1_1);
   fChain->SetBranchAddress("signalNcr_track1_2", &signalNcr_track1_2, &b_signalNcr_track1_2);
   fChain->SetBranchAddress("signalNcr_track1_3", &signalNcr_track1_3, &b_signalNcr_track1_3);
   fChain->SetBranchAddress("track0P4", &track0P4, &b_track0P4);
   fChain->SetBranchAddress("track1P4", &track1P4, &b_track1P4);
   fChain->SetBranchAddress("track0CombP4", &track0CombP4, &b_track0CombP4);
   fChain->SetBranchAddress("track1CombP4", &track1CombP4, &b_track1CombP4);
   fChain->SetBranchAddress("track0C14", &track0C14, &b_track0C14);
   fChain->SetBranchAddress("track1C14", &track1C14, &b_track1C14);
   fChain->SetBranchAddress("track0CombC14", &track0CombC14, &b_track0CombC14);
   fChain->SetBranchAddress("track1CombC14", &track1CombC14, &b_track1CombC14);
   fChain->SetBranchAddress("dZ_Row0", &dZ_Row0, &b_dZ_Row0);
   fChain->SetBranchAddress("dRPhi_Row0", &dRPhi_Row0, &b_dRPhi_Row0);
   fChain->SetBranchAddress("dZ_Row62", &dZ_Row62, &b_dZ_Row62);
   fChain->SetBranchAddress("dRPhi_Row62", &dRPhi_Row62, &b_dRPhi_Row62);
   fChain->SetBranchAddress("dZ_Row125", &dZ_Row125, &b_dZ_Row125);
   fChain->SetBranchAddress("dRPhi_Row125", &dRPhi_Row125, &b_dRPhi_Row125);
   fChain->SetBranchAddress("dZ_Row158", &dZ_Row158, &b_dZ_Row158);
   fChain->SetBranchAddress("dRPhi_Row158", &dRPhi_Row158, &b_dRPhi_Row158);
   fChain->SetBranchAddress("v0.fRr", &v0_fRr, &b_v0_fRr);
   fChain->SetBranchAddress("phi0_ROC0", &phi0_ROC0, &b_phi0_ROC0);
   fChain->SetBranchAddress("phi1_ROC0", &phi1_ROC0, &b_phi1_ROC0);
   fChain->SetBranchAddress("phi0_ROC1", &phi0_ROC1, &b_phi0_ROC1);
   fChain->SetBranchAddress("phi1_ROC1", &phi1_ROC1, &b_phi1_ROC1);
   fChain->SetBranchAddress("phi0_ROC2", &phi0_ROC2, &b_phi0_ROC2);
   fChain->SetBranchAddress("phi1_ROC2", &phi1_ROC2, &b_phi1_ROC2);
   Notify();
}

Bool_t V0FlatAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void V0FlatAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t V0FlatAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef V0FlatAna_cxx
