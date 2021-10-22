
#include "TObject.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSpline.h>
#include <TVectorD.h>
#include <vector>
#include "V0FlatAna.h"
#include "TrackFlatAna.h"
//#include "AliESDtrack.h"
class TH2F;
class TH1F;
class TF1;

class TArrayF;
template <class X>
class THnSparseT;
typedef class THnSparseT<TArrayF> THnSparseF;




class AliSkimmedDataAnalysisMaker :public TObject {

  
 public:
  static void doTest();
  AliSkimmedDataAnalysisMaker();
  //const char*
  virtual ~AliSkimmedDataAnalysisMaker();

  void bookHistogram();
  void WriteHistogram();
  void execute();
  void BinLogAxis(THnSparseF *h, Int_t axisNumber);
  void SetAxisNamesFromTitle(const THnSparseF *h);
  void TreeV0_BBFitAnalysis(TString filename_v0);
  void TreePrimary_BBFitAnalysis(TString filename_track);

  protected:

  static Double_t Lund(Double_t* xx, Double_t* par);
  static Double_t SaturatedLund(Double_t* xx, Double_t* par);
  static Double_t Aleph(Double_t* xx, Double_t* par);
  

 private:

  TFile * mOutputFile_BBFitAnalysis;      //! dummy

 TTree* fTree2;                           //! dummy
 //TVectorD tpctof(5);

//tree varialbe BB fit check
  const Int_t g=5;
  Double_t p;
  Double_t oneoverpt;
  Double_t rawtpcsignal;
  Double_t rawtpcsignal_leg;
  Double_t itssignal;
  Double_t pileupcor_tpcsignal;
  Double_t pileupcor_tpcsignal_leg;
  TVectorD tofNsig= TVectorD(9);
  TVectorD tpcNsig= TVectorD(9);
  //TVectorD tpcnsigma(nSpecies);
  Double_t BG;
  Double_t eta;
  Double_t tgl;
  Double_t multiplicty;
  Double_t centrality;
  Double_t PDGcode;
  Double_t isV0;
  Double_t V0Pull;
  Double_t V0PullEff;
  Double_t isPileUp;
  Double_t dEdxExpected_SatLund_woDeut;
  Double_t dEdxExpected_SatLund_wDeut;
  Double_t dEdxExpected_ALEPH_woDeut;
  Double_t dEdxExpected_ALEPH_wDeut;
  Double_t dEdxExpected_SatLund_woDeut_pileUpCut;
  Double_t tpcnsigma;
  Double_t tofnsigma;
  Double_t shiftM;
  Double_t nPileUpPrim;
  Double_t nPileUpSum;
  Double_t primMult;
  Double_t tpcClusterMult;
  Double_t pileUpOK;
  Double_t multSSD;
  Double_t multSDD;
  Double_t SignalTot0;
  Double_t SignalTot1;
  Double_t SignalTot2;
  Double_t SignalTot3;
  Double_t SignalMax0;
  Double_t SignalMax1;
  Double_t SignalMax2;
  Double_t SignalMax3;
  Double_t SignalNcl0;
  Double_t SignalNcl1;
  Double_t SignalNcl2;
  Double_t SignalNcl3;
  Double_t SignalNcr0;
  Double_t SignalNcr1;
  Double_t SignalNcr2;
  Double_t SignalNcr3;
  Double_t tpc_ncls;
  Double_t its_ncls;
  Double_t tpc_chi2;
  Double_t its_chi2;
  Double_t dca_r;
  Double_t dca_tpc_r;
  Double_t dca_z;
  Double_t dca_tpc_z;
  Double_t fITSClusterMap;
  Double_t fITSClusterMap_leg;
  Double_t fAlpha;
  Double_t K0PullEff;
  Double_t EPullEff;
  Double_t LPullEff;
  Double_t ALPullEff;
  Double_t fSigned1Pt;
  Double_t selectionPtMask;
  Double_t selectionPIDMask;
  // comined information
   Double_t        trackPComb;
   Double_t        trackPtComb;
   Double_t        tglComb;
   // local phi information
   Double_t        phi_ROC0;
   Double_t        phi_ROC1;
   Double_t        phi_ROC2;
   // sector information if exist
   Double_t        dSector0;
   Double_t        dSector1;
   Double_t        dSector2;
   //MC true information if exist
   Double_t        ptTPCMC;
   Double_t        ptMC;
   Double_t        tglMC;
   Double_t        fPdgCodeMC;
   //
   // 2 track distance from the V0 finder
   Double_t        dZ_Row0;
   Double_t        dRPhi_Row0;
   Double_t        dZ_Row62;
   Double_t        dRPhi_Row62;
   Double_t        dZ_Row125;
   Double_t        dRPhi_Row125;
   Double_t        dZ_Row158;
   Double_t        dRPhi_Row158;
   Double_t        v0_fRr;



  ClassDef(AliSkimmedDataAnalysisMaker,1);
};
