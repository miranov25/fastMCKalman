# include fastMCKalman to the path
# import sys,os; sys.path.insert(1, os.environ[f"fastMCKalman"]+'/fastMCKalman/MC/');  from test_fastSimulation import *
import sys,os;
sys.path.insert(1, os.environ[f"fastMCKalman"]+'/fastMCKalman/MC/');
import ROOT
from ROOT import TFile, gSystem
from fastSimulation import *
from RootInteractive.Tools.aliTreePlayer import *
from RootInteractive.MLpipeline.MIForestErrPDF import *
from RootInteractive.Tools.RDataFrameTools import filterRDFColumns
import awkward._v2 as ak

# load tree
def loadTree(inputData = "fastParticle.list"):
    ROOT.initTreeFast(inputData)
    tree=ROOT.treeFast
    return tree

# define aliases
def makeAliases(tree):
    tree.SetAlias("X0","geom.fLayerX0[0]")
    tree.SetAlias("sigmaRPhi","geom.fLayerResolRPhi[0]")
    tree.SetAlias("sigmaZ","geom.fLayerResolZ[0]")
    for pType in ["In","Out","Refit","MC"]:
        tree.SetAlias(f"statusMaskFull{pType}",f"partFull.fStatusMask{pType}")
        tree.SetAlias(f"NPointsFull{pType}",f"partFull.fNPoints{pType}")
        tree.SetAlias(f"dEdx{pType}",f"AliExternalTrackParam::BetheBlochSolid(partFull.fParam{pType}[].GetP()/partFull.fMassMC)")
        for iPar in [0,1,2,3,4]:
            tree.SetAlias(f"paramFull{pType}{iPar}",f"partFull.fParam{pType}[].fP[{iPar}]")
            iCov=ROOT.AliExternalTrackParam.GetIndex(iPar, iPar)
            tree.SetAlias(f"pullFull{pType}{iPar}",
                          f"(partFull.fParam{pType}[].fP[{iPar}]-partFull.fParamMC[].fP[{iPar}])/sqrt(partFull.fParam{pType}[].fC[{iCov}])")
            tree.SetAlias(f"deltaFull{pType}{iPar}",
                          f"(partFull.fParam{pType}[].fP[{iPar}]-partFull.fParamMC[].fP[{iPar}])")
            tree.SetAlias(f"sqrtCovFull{pType}{iPar}",
                          f"partFull.fParam{pType}[].fC[{iCov}]")


def loadPanda(tree,entries):
    variables=[".*pullFull.*",".*statusMaskFull.*",".*paramFull.*",".*NPointsFull.*",".*dEdx.*",".*deltaFull.*", ".*CovFull.*",
               ".*gx.*",".*gy.*",".*gz.*", "X0", "sigmaRPhi","sigmaZ",
               "ptMC","tglMC","fPdgCodeMC","fMassMC","pidCode","charge"
               # to add  lever arm at given layer
               ]
    exclude=[".*pullFullMC.*",".*deltaFullMC.*",".*CovFullMC.*",".*dEdxExp*"]
    df=tree2Panda(tree,variables, "partFull.fLengthIn>5",columnMask=[["_fElements",""]],exclude=exclude,nEntries=entries)
    df["statusMaskFullRefit"]=df["statusMaskFullRefit"].astype("int")
    return df


def loadData(inputData = "fastParticle.list"):
    tree= loadTree(inputData)
    makeAliases(tree)
    df=loadPanda(tree)
    return df

def makeRegression(df):
    varIn=["fPdgCodeMC","paramFullMC2","paramFullMC3","paramFullMC4","NPointsFullRefit"]
    target="pullFullRefit4"
    n_estimators=200
    n_jobs=100
    nPoints=500000
    max_depthBase=14
    max_samples=0.05

    regressor = RandomForestRegressor(n_estimators =n_estimators,n_jobs=n_jobs,max_depth=max_depthBase,max_samples=max_samples)
    #dfFit=df.query("(statusMaskFullRefit&0x2000)>0")
    dfFit0=df[((df["statusMaskFullRefit"]&0x1000)>0)].query(f"abs({target})<10")
    dfFit1=dfFit0.sample(frac=0.05).sort_index()
    #
    #
    nAll=dfFit1.shape[0]
    regressor.fit(dfFit1[varIn][:nAll//2], dfFit1[target][:nAll//2])
    dfFit1[f"{target}Pred0"]=regressor.predict(dfFit1[varIn])
    #
    ((dfFit1[f"{target}"]-dfFit1[f"{target}Pred0"])[nAll//2:]).std()
    ((dfFit1[f"{target}"]-dfFit1[f"{target}Pred0"])[:nAll//2]).std()



def loadCode():
    ROOT.gInterpreter.ProcessLine("""gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");""")
    ROOT.gInterpreter.ProcessLine(""".L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g""")
    ROOT.gInterpreter.ProcessLine(""".L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g""")
    ROOT.gInterpreter.ProcessLine(""".L $fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g""")


def loadRDF(input="fastParticle.list",verbosity=0, doTest=True):
    """

    :param input:
    :param verbosity:
    :return:
    example:
    input="fastParticle.list"; verbosity=0; doTest=True
    """

    from RootInteractive.Tools.RDataFrameTools import filterRDFColumns
    tree=loadTree(input)
    tree.SetCacheSize(40000000)
    if doTest == False: ROOT.EnableImplicitMT(10)
    rdf1 =ROOT. makeDataFrame(tree)
    if verbosity>0: print(rdf1.GetColumnNames())
    # define alias namess for some columns with dots not supported for the
    varListAlias=   filterRDFColumns(rdf1, ["partFull.*Status.*","partFull.*NPoint.*"],[],[".*"],[],verbosity)
    for var in varListAlias:
        aliasName=var.replace("partFull.f","")
        #rdf1=rdf1.Alias(aliasName,var)
        rdf1=rdf1.Define(aliasName,var)

    deltaListAlias=   filterRDFColumns(rdf1, ["delta.*"],[],[".*"],[],verbosity)

    for var in deltaListAlias:
        pullName=var.replace("delta","pull")
        covarName=var.replace("delta","covar")
        rdf1=rdf1.Define(pullName,var+"/"+covarName)
    #
    #
    varList = filterRDFColumns(rdf1, ["param.*","covar.*","delta.*",".*pid.*","charge",".*Status.*",".*NPoi.*",".*dEdx.*","pull.*",".*X0.*",".*sigma.*"],
                               ["part.*Para.*","geom.*","part.*",".*InRot.*" ],[".*"],[".*AliExternal.*","Long64.*","Long.*"], verbose=verbosity)


    if doTest:
        rdfTest=rdf1.Range(0,10)
        rdfTest.Snapshot("testVarRDF","testVarRDF.root", varList)
        array = ak.from_rdataframe(rdfTest, columns=varList)
        df=ak.to_dataframe(array)
        print(df.head(5),df.shape)
    #

    import awkward._v2 as ak
    array = ak.from_rdataframe(rdf1, columns=varList)
    df=ak.to_dataframe(array)
    return df