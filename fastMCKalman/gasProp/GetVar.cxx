//#include "GasUtils.h"
#include "gasElement.h"
#include "gasMixer.h"
#include "GetVar.h"

namespace GetVar{
  TTree * tree=0;
  std::map<std::string,AliNDLocalRegression*> mapRegressionDtStr;
  std::map<int,TString> mapId;
  std::map<int,TVectorF> mapParam;
  std::map<int,AliNDLocalRegression*> mapRegressionDtInt;
}

bool GetVar::Init(){
  static bool isInitialized=kFALSE;
  if (isInitialized==false){
    gSystem->Load("libgasElement.so");
    gSystem->Load("libgasMixer.so");
    gSystem->Load("libGasUtils.so");
    isInitialized=kTRUE;
  }
  return isInitialized;
}

bool GetVar::CreateOrganicCompound(gasMixer & gQuencher, int nC, int nH, int nF, gasMixer gC, gasMixer gH, gasMixer gF)
{
    if(nC==1) gQuencher=gQuencher+gC;
    else if (nC>1) gQuencher=gQuencher+gC*nC;
    if(nH==1) gQuencher=gQuencher+gH;
    else if (nH>1) gQuencher=gQuencher+gH*nH;
    if(nF==1) gQuencher=gQuencher+gF;
    else if (nF>1) gQuencher=gQuencher+gF*nF;
    if(nH==0&&nC==0&&nF==0) 
    {
        std::cout<<"Can't have a null quencher!"<<std::endl;
        return 0;
    }
    return 1;
}

float GetVar::GetSigmaT(int fquencher, int p, int nC, int nH, int Zgas, int Agas, int E, int B)
{ 
    //Setup input file and tree
    TFile * fin=new TFile("all_merged.root");
    TTree * gastree=dynamic_cast<TTree*>(fin->Get("gastree"));

    //Define the cut to get the desired value
    TString cut;
    cut = Form("fquencher==%d && p==%d && nC==%d && nH==%d && Zgas==%d && Agas==%d && E==%d && B==%d",fquencher,p,nC,nH,Zgas,Agas,E,B);

    //obtain the ttree entry correspondent to our parameters
    gastree->Draw("Entry$>>hist(Entries$,0,Entries$)",cut,"goff");
    TH1I *hist = (TH1I*)gDirectory->Get("hist");
    Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
    //std::cout<<"nEntries: "<<hist->GetEntries()<<"; Entry: "<<iEntry<<std::endl;

    if (hist->GetEntries()==0)
    {
        std::cout<<"No entry matches your criteria"<<std::endl;
        return 0;
    }
    if (hist->GetEntries()==0)
    {
        std::cout<<"More than one entry matches your criteria"<<std::endl;
        return 0;
    }
    
    
    //Setup the sigma diffusion branch
    float dt;
    gastree->SetBranchAddress("dt",&dt);
    gastree->GetEntry(iEntry);

    return dt;
}

float GetVar::GetSigmaTErr(int fquencher, int p, int nC, int nH, int Zgas, int Agas, int E, int B)
{ 
    //Setup input file and tree
    TFile * fin=new TFile("all_merged.root");
    TTree * gastree=dynamic_cast<TTree*>(fin->Get("gastree"));

    //Define the cut to get the desired value
    TString cut;
    cut = Form("fquencher==%d && p==%d && nC==%d && nH==%d && Zgas==%d && Agas==%d && E==%d && B==%d",fquencher,p,nC,nH,Zgas,Agas,E,B);

    //obtain the ttree entry correspondent to our parameters
    gastree->Draw("Entry$>>hist(Entries$,0,Entries$)",cut,"goff");
    TH1I *hist = (TH1I*)gDirectory->Get("hist");
    Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
    //std::cout<<"nEntries: "<<hist->GetEntries()<<"; Entry: "<<iEntry<<std::endl;

    if (hist->GetEntries()==0)
    {
        std::cout<<"No entry matches your criteria"<<std::endl;
        return 0;
    }
    if (hist->GetEntries()==0)
    {
        std::cout<<"More than one entry matches your criteria"<<std::endl;
        return 0;
    }
    
    
    //Setup the sigma diffusion branch
    float dterr;
    gastree->SetBranchAddress("dterr",&dterr);
    gastree->GetEntry(iEntry);

    return dterr;
}


float GetVar::GetX0Rel(int fquencher, int p , int nC, int nH, int Zgas, int Agas) //fquencher in percentage, p in mBar2
{
    Init();

    //Define building blocks
    const gasMixer gC(6, 12), gH(1, 1), gF(9,19);
    
    //Define main gas
    gasMixer gMain(Zgas,Agas);
    
    //Define quencher
    gasMixer gQuencher;
    CreateOrganicCompound(gQuencher,nC,nH,0,gC,gH,gF);
    
    //Define mixture
    const gasMixer gTotal(gMain, gQuencher, fquencher/100., p/1000.);

    //gMain.Print();
    //gQuencher.Print();
    //gTotal.Print();
    
    return gTotal.getX0();    
}

float GetVar::GetX0F(int fquencher, int p , int nC, int nH, int nF, int Zgas, int Agas, int mC, int mH, int mF) //fquencher in percentage, p in mBar2
{
  Init();

    //Define building blocks
    const gasMixer gC(6, 12), gH(1, 1), gF(9,19);
    
    //Define main gas
    gasMixer gMain(Zgas,Agas);
    if (mC!=0 || mH!=0 || mF!=0) CreateOrganicCompound(gMain,mC,mH,mF,gC,gH,gF);
    
    //Define quencher
    gasMixer gQuencher;
    CreateOrganicCompound(gQuencher,nC,nH,nF,gC,gH,gF);
    
    //Define mixture
    const gasMixer gTotal(gMain, gQuencher, fquencher/100., p/1000.);

    //gMain.Print();
    //gQuencher.Print();
    //gTotal.Print();
    
    return gTotal.getX0();    
}


void GetVar::makeGasFunction(int minEntries, const char * output){
  if (output) gSystem->RedirectOutput(output);
  int nH[]={1,4,6,8,10};
  int nC[]={1,2,3,4};
  int Agas[]={4,40};
  int Zgas[]={2,18};
  float B[]={0,0.5};
  //
  tree->SetAlias("normdT","sqrt(p)*dt");
  Int_t id=0;
  for (int iH=0; iH<5; iH++)
    for (int iC=0; iC<4; iC++)
      for (int iAgas=0; iAgas<2;iAgas++)
        for (int iZgas=0; iZgas<2;iZgas++)
          for (int iB=0; iB<2; iB++){
            TString selection=TString::Format("sim==1&&nH==%d&&nC==%d&&Agas==%d&&Zgas==%d&&B==%0.2f",nH[iH],nC[iC], Agas[iAgas],Zgas[iZgas],B[iB]);
;
            Int_t entries=tree->GetEntries(selection);
            //gSystem->RedirectOutput(0,0);
            if (entries<minEntries) continue;
            entries=tree->Draw("fquencher:E/p",selection,"goffpara");
            double range[4];
            range[0]=GetVar::tree->GetVal(0)[TMath::LocMin(entries, GetVar::tree->GetVal(0))];
            range[1]=GetVar::tree->GetVal(0)[TMath::LocMax(entries, GetVar::tree->GetVal(0))];
            range[2]=GetVar::tree->GetVal(1)[TMath::LocMin(entries, GetVar::tree->GetVal(1))];
            range[3]=GetVar::tree->GetVal(1)[TMath::LocMax(entries, GetVar::tree->GetVal(1))];
            //
            //TString range0="(20,0,100,40,0.005,1)";
            TString range0=TString::Format("(20,%3f,%3f,40,%3f,%3F)",range[0],range[1],TMath::Max(range[2],0.005),TMath::Min(range[3],1.));
            TString kernel0="10:0.01+(E/p)*0.05";
            TString valueString=TString::Format("normdT:0.1");
            TString fitString=TString::Format("fquencher:E/p");
            printf("%s\t%d\t%s\n",selection.Data(),entries,range0.Data());
            //
            AliNDLocalRegression *fit = AliNDLocalRegression::MakeRegression(tree,Form("normdTFit_%d",id),range0, valueString.Data(), fitString,selection+"&&E/p<1&&E/p>0.005", kernel0,0.001);
            mapId[id]=selection.Data();
            TVectorF &test = mapParam[id];
            mapParam[id].ResizeTo(5);
            mapParam[id][0]=nC[iC];
            mapParam[id][1]=nH[iH];
            mapParam[id][2]=Zgas[iZgas];
            mapParam[id][3]=Agas[iAgas];
            mapParam[id][4]=B[iB];
            mapRegressionDtStr[selection.Data()]=fit;
            mapRegressionDtInt[id]=fit;
            id++;
          }
    if (output)  gSystem->RedirectOutput(0,0);
}


double GetVar::getDt(int id,  float fquencher, float ep){
  AliNDLocalRegression * regresion=mapRegressionDtInt[id];
  if (regresion==0) return 0;
  Double_t inData[2]={fquencher,ep};
  return regresion->Eval(inData);
}


char * GetVar::getIdName(int id){
    return (char*) mapId[id].Data();
}


void GetVar::initTree(){ //TODO - load from the source github

  //using GetVar::tree;
  TFile *f = TFile::Open("all_merged.root");
  if (f== nullptr) f = TFile::Open(gSystem->ExpandPathName("$fastMCKalman/fastMCKalman/gasProp/all_merged.root"));
  tree = (TTree*) f->Get("gastree");
}

float GetVar::GetDensity(int fquencher, int p , int nC, int nH, int Zgas, int Agas) //fquencher in percentage, p in mBar2
{
    Init();

    //Define building blocks
    const gasMixer gC(6, 12), gH(1, 1), gF(9,19);
    
    //Define main gas
    gasMixer gMain(Zgas,Agas);
    
    //Define quencher
    gasMixer gQuencher;
    CreateOrganicCompound(gQuencher,nC,nH,0,gC,gH,gF);
       

    // returns an estimated density in g/cm3 for gas mix at given pressure p in mbar
    
    // constants in g / cm3 - from wolframalpha at STP := 1000mbar, 298K
    double noble_stp;
    double quencher_stp;

    switch ((int)gMain.getZ())
    {
    case 18:  // Ar
      noble_stp = 16.13e-4;
      break;
    
    case 2:  // He
      noble_stp = 1.615e-4;
    
    default:
      std::cout << "Unknown noble gas requested in gasMixer::getGasDensity - aborting" << std::endl;
      printf("nC\t%d\tnH\t%d\tZ\t%d\tA\t%d\n",nC,nH,Zgas,Agas);
      exit(1);
      break;
    }

    switch ((int)gQuencher.getH())
    {
    case 4:  // CH4
      quencher_stp = 6.486e-4;
      break;
    case 6:  // C2H6
      quencher_stp = 12.23e-4;
      break;
    case 8:  // C3H8
      quencher_stp = 18.09e-4;  // could be replaced by a lambda that scales linearly between stp and 10bar
      break;

    default:
      std::cout << "Unknown quenching gas requested in gasMixer::getGasDensity - aborting" << std::endl;
      exit(1);
      break;
    }

    // scale using ideal gas law
    double density = (1-fquencher/100.) * noble_stp + fquencher/100. * quencher_stp;
    density *= p * 1e-3;

    return density;
}

float GetVar::GetX0(int fquencher, int p, int nC, int nH, int Zgas, int Agas)
{ 
  float X0;
  X0=(GetVar::GetX0Rel(fquencher,p,nC,nH,Zgas,Agas)/GetVar::GetDensity(fquencher,p,nC,nH,Zgas,Agas));
  return X0; 
}

float GetVar::GetNumberOfX0(int fquencher, int p , int nC, int nH, int Zgas, int Agas, float length){ //length must be in cm
  double X0;
  X0=length/(GetVar::GetX0Rel(fquencher,p,nC,nH,Zgas,Agas)/GetVar::GetDensity(fquencher,p,nC,nH,Zgas,Agas));
  return X0; 
}





