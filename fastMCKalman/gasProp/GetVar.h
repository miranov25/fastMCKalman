#ifndef _GETVAR_
#define _GETVAR_

#include "stdio.h"
//#include "TAxis.h"
#include "TFile.h"
//#include "TCanvas.h"
//#include "TColor.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TGraphAsymmErrors.h"
//#include "TLegend.h"
//#include "TLine.h"
//#include "TList.h"
#include "TH1I.h"
#include "TSystem.h"
#include "TTree.h"
#include "GasUtils.h"
#include "gasElement.h"
#include "gasMixer.h"
//#include "GetVar.h"

//using namespace GasUtils;

namespace GetVar
{

bool CreateOrganicCompound(gasMixer & gQuencher, int nC, int nH, int nF, gasMixer gC, gasMixer gH, gasMixer gF);


float GetSigmaT(int fquencher=10, int p=1000, int nC=1, int nH=4, int Zgas=18, int Agas=40, int E=165, int B=0);


float GetSigmaTErr(int fquencher=10, int p=1000, int nC=1, int nH=4, int Zgas=18, int Agas=40, int E=165, int B=0);


float GetX0(int fquencher=10, int p=1000 , int nC=1, int nH=4, int Zgas=18, int Agas=40); //fquencher in percentage, p in mBar2

float GetX0F(int fquencher, int p , int nC, int nH, int nF, int Zgas, int Agas, int mC, int mH, int mF);

}

#endif