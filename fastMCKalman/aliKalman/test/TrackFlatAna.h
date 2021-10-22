//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 21 11:10:48 2021 by ROOT version 6.20/08
// from TTree CleanTrackFlat/CleanTrackFlat
// found on file: CleanTrack.root
//////////////////////////////////////////////////////////

#ifndef TrackFlatAna_h
#define TrackFlatAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class TrackFlatAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        trackPComb;
   Double_t        trackPtComb;
   Double_t        tglComb;
   Double_t        phi_ROC0;
   Double_t        phi_ROC1;
   Double_t        phi_ROC2;
   Double_t        fTPCsignal;
   Double_t        esdTrack_fITSsignal;
   Double_t        esdTrack_fTRDsignal;
   Double_t        esdTrack_fTPCsignalN;
   Double_t        trackP;
   Double_t        trackPt;
   Double_t        tgl;
   Double_t        trackEta;
   Double_t        trackTPCchi2;
   Double_t        trackITSchi2;
   Double_t        trackTPC_cls;
   Double_t        trackITS_cls;
   Double_t        ntracks;
   Double_t        trackstatus;
   Double_t        track_tpcNsigma_el;
   Double_t        track_tpcNsigma_pi;
   Double_t        track_tpcNsigma_ka;
   Double_t        track_tpcNsigma_pro;
   Double_t        track_tpcNsigma_deut;
   Double_t        track_tpcNsigma_tri;
   Double_t        track_tpcNsigma_He3;
   Double_t        tracktofNsigmaElectron;
   Double_t        tracktofNsigmaPion;
   Double_t        tracktofNsigmaKaon;
   Double_t        tracktofNsigmaProton;
   Double_t        tracktofNsigmaDeuteron;
   Double_t        tracktofNsigmaTrition;
   Double_t        tracktofNsigmaHelium3;
   Double_t        track_ExpectedTPCSignalV0_el;
   Double_t        track_ExpectedTPCSignalV0_pi;
   Double_t        track_ExpectedTPCSignalV0_ka;
   Double_t        track_ExpectedTPCSignalV0_pro;
   Double_t        track_ExpectedTPCSignalV0_deut;
   Double_t        track_ExpectedTPCSignalV0_tri;
   Double_t        track_ExpectedTPCSignalV0_He3;
   Double_t        track_CorrectedTPCSignalV0_el;
   Double_t        track_CorrectedTPCSignalV0_pi;
   Double_t        track_CorrectedTPCSignalV0_ka;
   Double_t        track_CorrectedTPCSignalV0_pro;
   Double_t        track_CorrectedTPCSignalV0_deut;
   Double_t        track_CorrectedTPCSignalV0_tri;
   Double_t        track_CorrectedTPCSignalV0_He3;
   Double_t        track_GetPileupValue;
   Double_t        trackDCAxy;
   Double_t        trackDCAz;
   Double_t        nCrossRows;
   Double_t        trackDCAxy_TPC;
   Double_t        trackDCAz_TPC;
   Double_t        fAlpha;
   Double_t        fITSClusterMap;
   Double_t        fSigned1Pt;
   Double_t        run;
   Double_t        intrate;
   Double_t        timeStampS;
   Double_t        timestamp;
   Double_t        bField;
   Double_t        triggerMask;
   Double_t        gid;
   Double_t        shiftA;
   Double_t        shiftC;
   Double_t        shiftM;
   Double_t        nPileUpPrim;
   Double_t        nPileUpSum;
   Double_t        primMult;
   Double_t        tpcClusterMult;
   Double_t        pileUpOK;
   Double_t        multSSD;
   Double_t        multSDD;
   Double_t        multSPD;
   Double_t        multV0;
   Double_t        multT0;
   Double_t        spdvz;
   Double_t        itsTracklets;
   Double_t        tpcMult;
   Double_t        tpcTrackBeforeClean;
   Double_t        centV0;
   Double_t        centITS0;
   Double_t        centITS1;
   Double_t        TPCRefit;
   Double_t        ITSRefit;
   Double_t        Nucleitrigger_OFF;
   Double_t        track_hasTOF;
   Double_t        isMinBias;
   Double_t        isPileUp;
   Double_t        selectionPtMask;
   Double_t        selectionPIDMask;
   Double_t        logSignalTot0;
   Double_t        logSignalTot1;
   Double_t        logSignalTot2;
   Double_t        logSignalTot3;
   Double_t        logSignalMax0;
   Double_t        logSignalMax1;
   Double_t        logSignalMax2;
   Double_t        logSignalMax3;
   Double_t        signalNcl0;
   Double_t        signalNcl1;
   Double_t        signalNcl2;
   Double_t        signalNcl3;
   Double_t        signalNcr0;
   Double_t        signalNcr1;
   Double_t        signalNcr2;
   Double_t        signalNcr3;
   Double_t        dSector0;
   Double_t        dSector1;
   Double_t        dSector2;
   Double_t        ptTPCMC;
   Double_t        ptMC;
   Double_t        tglMC;
   Double_t        fPdgCode;
   Double_t        T;

   // List of branches
   TBranch        *b_trackPComb;   //!
   TBranch        *b_trackPtComb;   //!
   TBranch        *b_tglComb;   //!
   TBranch        *b_phi_ROC0;   //!
   TBranch        *b_phi_ROC1;   //!
   TBranch        *b_phi_ROC2;   //!
   TBranch        *b_fTPCsignal;   //!
   TBranch        *b_esdTrack_fITSsignal;   //!
   TBranch        *b_esdTrack_fTRDsignal;   //!
   TBranch        *b_esdTrack_fTPCsignalN;   //!
   TBranch        *b_trackP;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_tgl;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackTPCchi2;   //!
   TBranch        *b_trackITSchi2;   //!
   TBranch        *b_trackTPC_cls;   //!
   TBranch        *b_trackITS_cls;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_trackstatus;   //!
   TBranch        *b_track_tpcNsigma_el;   //!
   TBranch        *b_track_tpcNsigma_pi;   //!
   TBranch        *b_track_tpcNsigma_ka;   //!
   TBranch        *b_track_tpcNsigma_pro;   //!
   TBranch        *b_track_tpcNsigma_deut;   //!
   TBranch        *b_track_tpcNsigma_tri;   //!
   TBranch        *b_track_tpcNsigma_He3;   //!
   TBranch        *b_tracktofNsigmaElectron;   //!
   TBranch        *b_tracktofNsigmaPion;   //!
   TBranch        *b_tracktofNsigmaKaon;   //!
   TBranch        *b_tracktofNsigmaProton;   //!
   TBranch        *b_tracktofNsigmaDeuteron;   //!
   TBranch        *b_tracktofNsigmaTrition;   //!
   TBranch        *b_tracktofNsigmaHelium3;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_el;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_pi;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_ka;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_pro;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_deut;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_tri;   //!
   TBranch        *b_track_ExpectedTPCSignalV0_He3;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_el;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_pi;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_ka;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_pro;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_deut;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_tri;   //!
   TBranch        *b_track_CorrectedTPCSignalV0_He3;   //!
   TBranch        *b_track_GetPileupValue;   //!
   TBranch        *b_trackDCAxy;   //!
   TBranch        *b_trackDCAz;   //!
   TBranch        *b_nCrossRows;   //!
   TBranch        *b_trackDCAxy_TPC;   //!
   TBranch        *b_trackDCAz_TPC;   //!
   TBranch        *b_fAlpha;   //!
   TBranch        *b_fITSClusterMap;   //!
   TBranch        *b_fSigned1Pt;   //!
   TBranch        *b_run;   //!
   TBranch        *b_intrate;   //!
   TBranch        *b_timeStampS;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_bField;   //!
   TBranch        *b_triggerMask;   //!
   TBranch        *b_gid;   //!
   TBranch        *b_shiftA;   //!
   TBranch        *b_shiftC;   //!
   TBranch        *b_shiftM;   //!
   TBranch        *b_nPileUpPrim;   //!
   TBranch        *b_nPileUpSum;   //!
   TBranch        *b_primMult;   //!
   TBranch        *b_tpcClusterMult;   //!
   TBranch        *b_pileUpOK;   //!
   TBranch        *b_multSSD;   //!
   TBranch        *b_multSDD;   //!
   TBranch        *b_multSPD;   //!
   TBranch        *b_multV0;   //!
   TBranch        *b_multT0;   //!
   TBranch        *b_spdvz;   //!
   TBranch        *b_itsTracklets;   //!
   TBranch        *b_tpcMult;   //!
   TBranch        *b_tpcTrackBeforeClean;   //!
   TBranch        *b_centV0;   //!
   TBranch        *b_centITS0;   //!
   TBranch        *b_centITS1;   //!
   TBranch        *b_TPCRefit;   //!
   TBranch        *b_ITSRefit;   //!
   TBranch        *b_Nucleitrigger_OFF;   //!
   TBranch        *b_track_hasTOF;   //!
   TBranch        *b_isMinBias;   //!
   TBranch        *b_isPileUp;   //!
   TBranch        *b_selectionPtMask;   //!
   TBranch        *b_selectionPIDMask;   //!
   TBranch        *b_logSignalTot0;   //!
   TBranch        *b_logSignalTot1;   //!
   TBranch        *b_logSignalTot2;   //!
   TBranch        *b_logSignalTot3;   //!
   TBranch        *b_logSignalMax0;   //!
   TBranch        *b_logSignalMax1;   //!
   TBranch        *b_logSignalMax2;   //!
   TBranch        *b_logSignalMax3;   //!
   TBranch        *b_signalNcl0;   //!
   TBranch        *b_signalNcl1;   //!
   TBranch        *b_signalNcl2;   //!
   TBranch        *b_signalNcl3;   //!
   TBranch        *b_signalNcr0;   //!
   TBranch        *b_signalNcr1;   //!
   TBranch        *b_signalNcr2;   //!
   TBranch        *b_signalNcr3;   //!
   TBranch        *b_dSector0;   //!
   TBranch        *b_dSector1;   //!
   TBranch        *b_dSector2;   //!
   TBranch        *b_ptTPCMC;   //!
   TBranch        *b_ptMC;   //!
   TBranch        *b_tglMC;   //!
   TBranch        *b_fPdgCode;   //!
   TBranch        *b_T;   //!

   TrackFlatAna(TTree *tree=0);
   virtual ~TrackFlatAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TrackFlatAna_cxx
TrackFlatAna::TrackFlatAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CleanTrack.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("CleanTrack.root");
      }
      f->GetObject("CleanTrackFlat",tree);

   }
   Init(tree);
}

TrackFlatAna::~TrackFlatAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TrackFlatAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TrackFlatAna::LoadTree(Long64_t entry)
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

void TrackFlatAna::Init(TTree *tree)
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

   fChain->SetBranchAddress("trackPComb", &trackPComb, &b_trackPComb);
   fChain->SetBranchAddress("trackPtComb", &trackPtComb, &b_trackPtComb);
   fChain->SetBranchAddress("tglComb", &tglComb, &b_tglComb);
   fChain->SetBranchAddress("phi_ROC0", &phi_ROC0, &b_phi_ROC0);
   fChain->SetBranchAddress("phi_ROC1", &phi_ROC1, &b_phi_ROC1);
   fChain->SetBranchAddress("phi_ROC2", &phi_ROC2, &b_phi_ROC2);
   fChain->SetBranchAddress("fTPCsignal", &fTPCsignal, &b_fTPCsignal);
   fChain->SetBranchAddress("esdTrack.fITSsignal", &esdTrack_fITSsignal, &b_esdTrack_fITSsignal);
   fChain->SetBranchAddress("esdTrack.fTRDsignal", &esdTrack_fTRDsignal, &b_esdTrack_fTRDsignal);
   fChain->SetBranchAddress("esdTrack.fTPCsignalN", &esdTrack_fTPCsignalN, &b_esdTrack_fTPCsignalN);
   fChain->SetBranchAddress("trackP", &trackP, &b_trackP);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("tgl", &tgl, &b_tgl);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackTPCchi2", &trackTPCchi2, &b_trackTPCchi2);
   fChain->SetBranchAddress("trackITSchi2", &trackITSchi2, &b_trackITSchi2);
   fChain->SetBranchAddress("trackTPC_cls", &trackTPC_cls, &b_trackTPC_cls);
   fChain->SetBranchAddress("trackITS_cls", &trackITS_cls, &b_trackITS_cls);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("trackstatus", &trackstatus, &b_trackstatus);
   fChain->SetBranchAddress("track_tpcNsigma_el", &track_tpcNsigma_el, &b_track_tpcNsigma_el);
   fChain->SetBranchAddress("track_tpcNsigma_pi", &track_tpcNsigma_pi, &b_track_tpcNsigma_pi);
   fChain->SetBranchAddress("track_tpcNsigma_ka", &track_tpcNsigma_ka, &b_track_tpcNsigma_ka);
   fChain->SetBranchAddress("track_tpcNsigma_pro", &track_tpcNsigma_pro, &b_track_tpcNsigma_pro);
   fChain->SetBranchAddress("track_tpcNsigma_deut", &track_tpcNsigma_deut, &b_track_tpcNsigma_deut);
   fChain->SetBranchAddress("track_tpcNsigma_tri", &track_tpcNsigma_tri, &b_track_tpcNsigma_tri);
   fChain->SetBranchAddress("track_tpcNsigma_He3", &track_tpcNsigma_He3, &b_track_tpcNsigma_He3);
   fChain->SetBranchAddress("tracktofNsigmaElectron", &tracktofNsigmaElectron, &b_tracktofNsigmaElectron);
   fChain->SetBranchAddress("tracktofNsigmaPion", &tracktofNsigmaPion, &b_tracktofNsigmaPion);
   fChain->SetBranchAddress("tracktofNsigmaKaon", &tracktofNsigmaKaon, &b_tracktofNsigmaKaon);
   fChain->SetBranchAddress("tracktofNsigmaProton", &tracktofNsigmaProton, &b_tracktofNsigmaProton);
   fChain->SetBranchAddress("tracktofNsigmaDeuteron", &tracktofNsigmaDeuteron, &b_tracktofNsigmaDeuteron);
   fChain->SetBranchAddress("tracktofNsigmaTrition", &tracktofNsigmaTrition, &b_tracktofNsigmaTrition);
   fChain->SetBranchAddress("tracktofNsigmaHelium3", &tracktofNsigmaHelium3, &b_tracktofNsigmaHelium3);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_el", &track_ExpectedTPCSignalV0_el, &b_track_ExpectedTPCSignalV0_el);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_pi", &track_ExpectedTPCSignalV0_pi, &b_track_ExpectedTPCSignalV0_pi);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_ka", &track_ExpectedTPCSignalV0_ka, &b_track_ExpectedTPCSignalV0_ka);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_pro", &track_ExpectedTPCSignalV0_pro, &b_track_ExpectedTPCSignalV0_pro);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_deut", &track_ExpectedTPCSignalV0_deut, &b_track_ExpectedTPCSignalV0_deut);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_tri", &track_ExpectedTPCSignalV0_tri, &b_track_ExpectedTPCSignalV0_tri);
   fChain->SetBranchAddress("track_ExpectedTPCSignalV0_He3", &track_ExpectedTPCSignalV0_He3, &b_track_ExpectedTPCSignalV0_He3);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_el", &track_CorrectedTPCSignalV0_el, &b_track_CorrectedTPCSignalV0_el);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_pi", &track_CorrectedTPCSignalV0_pi, &b_track_CorrectedTPCSignalV0_pi);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_ka", &track_CorrectedTPCSignalV0_ka, &b_track_CorrectedTPCSignalV0_ka);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_pro", &track_CorrectedTPCSignalV0_pro, &b_track_CorrectedTPCSignalV0_pro);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_deut", &track_CorrectedTPCSignalV0_deut, &b_track_CorrectedTPCSignalV0_deut);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_tri", &track_CorrectedTPCSignalV0_tri, &b_track_CorrectedTPCSignalV0_tri);
   fChain->SetBranchAddress("track_CorrectedTPCSignalV0_He3", &track_CorrectedTPCSignalV0_He3, &b_track_CorrectedTPCSignalV0_He3);
   fChain->SetBranchAddress("track_GetPileupValue", &track_GetPileupValue, &b_track_GetPileupValue);
   fChain->SetBranchAddress("trackDCAxy", &trackDCAxy, &b_trackDCAxy);
   fChain->SetBranchAddress("trackDCAz", &trackDCAz, &b_trackDCAz);
   fChain->SetBranchAddress("nCrossRows", &nCrossRows, &b_nCrossRows);
   fChain->SetBranchAddress("trackDCAxy_TPC", &trackDCAxy_TPC, &b_trackDCAxy_TPC);
   fChain->SetBranchAddress("trackDCAz_TPC", &trackDCAz_TPC, &b_trackDCAz_TPC);
   fChain->SetBranchAddress("fAlpha", &fAlpha, &b_fAlpha);
   fChain->SetBranchAddress("fITSClusterMap", &fITSClusterMap, &b_fITSClusterMap);
   fChain->SetBranchAddress("fSigned1Pt", &fSigned1Pt, &b_fSigned1Pt);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("intrate", &intrate, &b_intrate);
   fChain->SetBranchAddress("timeStampS", &timeStampS, &b_timeStampS);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("bField", &bField, &b_bField);
   fChain->SetBranchAddress("triggerMask", &triggerMask, &b_triggerMask);
   fChain->SetBranchAddress("gid", &gid, &b_gid);
   fChain->SetBranchAddress("shiftA", &shiftA, &b_shiftA);
   fChain->SetBranchAddress("shiftC", &shiftC, &b_shiftC);
   fChain->SetBranchAddress("shiftM", &shiftM, &b_shiftM);
   fChain->SetBranchAddress("nPileUpPrim", &nPileUpPrim, &b_nPileUpPrim);
   fChain->SetBranchAddress("nPileUpSum", &nPileUpSum, &b_nPileUpSum);
   fChain->SetBranchAddress("primMult", &primMult, &b_primMult);
   fChain->SetBranchAddress("tpcClusterMult", &tpcClusterMult, &b_tpcClusterMult);
   fChain->SetBranchAddress("pileUpOK", &pileUpOK, &b_pileUpOK);
   fChain->SetBranchAddress("multSSD", &multSSD, &b_multSSD);
   fChain->SetBranchAddress("multSDD", &multSDD, &b_multSDD);
   fChain->SetBranchAddress("multSPD", &multSPD, &b_multSPD);
   fChain->SetBranchAddress("multV0", &multV0, &b_multV0);
   fChain->SetBranchAddress("multT0", &multT0, &b_multT0);
   fChain->SetBranchAddress("spdvz", &spdvz, &b_spdvz);
   fChain->SetBranchAddress("itsTracklets", &itsTracklets, &b_itsTracklets);
   fChain->SetBranchAddress("tpcMult", &tpcMult, &b_tpcMult);
   fChain->SetBranchAddress("tpcTrackBeforeClean", &tpcTrackBeforeClean, &b_tpcTrackBeforeClean);
   fChain->SetBranchAddress("centV0", &centV0, &b_centV0);
   fChain->SetBranchAddress("centITS0", &centITS0, &b_centITS0);
   fChain->SetBranchAddress("centITS1", &centITS1, &b_centITS1);
   fChain->SetBranchAddress("TPCRefit", &TPCRefit, &b_TPCRefit);
   fChain->SetBranchAddress("ITSRefit", &ITSRefit, &b_ITSRefit);
   fChain->SetBranchAddress("Nucleitrigger_OFF", &Nucleitrigger_OFF, &b_Nucleitrigger_OFF);
   fChain->SetBranchAddress("track_hasTOF", &track_hasTOF, &b_track_hasTOF);
   fChain->SetBranchAddress("isMinBias", &isMinBias, &b_isMinBias);
   fChain->SetBranchAddress("isPileUp", &isPileUp, &b_isPileUp);
   fChain->SetBranchAddress("selectionPtMask", &selectionPtMask, &b_selectionPtMask);
   fChain->SetBranchAddress("selectionPIDMask", &selectionPIDMask, &b_selectionPIDMask);
   fChain->SetBranchAddress("logSignalTot0", &logSignalTot0, &b_logSignalTot0);
   fChain->SetBranchAddress("logSignalTot1", &logSignalTot1, &b_logSignalTot1);
   fChain->SetBranchAddress("logSignalTot2", &logSignalTot2, &b_logSignalTot2);
   fChain->SetBranchAddress("logSignalTot3", &logSignalTot3, &b_logSignalTot3);
   fChain->SetBranchAddress("logSignalMax0", &logSignalMax0, &b_logSignalMax0);
   fChain->SetBranchAddress("logSignalMax1", &logSignalMax1, &b_logSignalMax1);
   fChain->SetBranchAddress("logSignalMax2", &logSignalMax2, &b_logSignalMax2);
   fChain->SetBranchAddress("logSignalMax3", &logSignalMax3, &b_logSignalMax3);
   fChain->SetBranchAddress("signalNcl0", &signalNcl0, &b_signalNcl0);
   fChain->SetBranchAddress("signalNcl1", &signalNcl1, &b_signalNcl1);
   fChain->SetBranchAddress("signalNcl2", &signalNcl2, &b_signalNcl2);
   fChain->SetBranchAddress("signalNcl3", &signalNcl3, &b_signalNcl3);
   fChain->SetBranchAddress("signalNcr0", &signalNcr0, &b_signalNcr0);
   fChain->SetBranchAddress("signalNcr1", &signalNcr1, &b_signalNcr1);
   fChain->SetBranchAddress("signalNcr2", &signalNcr2, &b_signalNcr2);
   fChain->SetBranchAddress("signalNcr3", &signalNcr3, &b_signalNcr3);
   fChain->SetBranchAddress("dSector0", &dSector0, &b_dSector0);
   fChain->SetBranchAddress("dSector1", &dSector1, &b_dSector1);
   fChain->SetBranchAddress("dSector2", &dSector2, &b_dSector2);
   fChain->SetBranchAddress("ptTPCMC", &ptTPCMC, &b_ptTPCMC);
   fChain->SetBranchAddress("ptMC", &ptMC, &b_ptMC);
   fChain->SetBranchAddress("tglMC", &tglMC, &b_tglMC);
   fChain->SetBranchAddress("fPdgCode", &fPdgCode, &b_fPdgCode);
   fChain->SetBranchAddress("T", &T, &b_T);
   Notify();
}

Bool_t TrackFlatAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TrackFlatAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TrackFlatAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TrackFlatAna_cxx
