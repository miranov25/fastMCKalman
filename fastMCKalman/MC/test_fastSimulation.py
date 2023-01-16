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
import awkward as ak
# memory/cpu monitoring
import os
import psutil
import time
import gc
pid = os.getpid()
python_process = psutil.Process(pid)
# load tree
def loadTree(inputData = "fastParticle.list"):
    ROOT.initTreeFast(inputData)
    #tree.SetCacheSize()
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


def loadRDF(input="fastParticle.list",verbosity=0, doTest=True, nThreads=0):
    """

    :param input:
    :param verbosity:
    :return:
    example:
    input="fastParticle.list"; verbosity=0; doTest=True; nThreads=0
    """

    from RootInteractive.Tools.RDataFrameTools import filterRDFColumns
    tree=loadTree(input)
    tree.SetCacheSize(40000)
    if nThreads >0:
        ROOT.EnableImplicitMT(nThreads)
    else:
        ROOT.DisableImplicitMT()
    rdf1 =ROOT. makeDataFrame(tree)
    if verbosity>0: print(rdf1.GetColumnNames())
    # define alias namess for some columns with dots not supported for the
    varListAlias=   filterRDFColumns(rdf1, ["partFull.*Status.*","partFull.*NPoint.*"],[],[".*"],[],verbosity)
    for var in varListAlias:
        aliasName=var.replace("partFull.f","")
        #rdf1=rdf1.Alias(aliasName,var)
        rdf1=rdf1.Define(aliasName,var)

    #deltaListAlias=   filterRDFColumns(rdf1, ["delta.*"],[],[".*"],[],verbosity)

    #for var in deltaListAlias:
    #    pullName=var.replace("delta","pull")
    #    covarName=var.replace("delta","covar")
    #    rdf1=rdf1.Define(pullName,var+"/"+covarName)
    #
    #
    varList = filterRDFColumns(rdf1,
                               ["param.*","covar.*","delta.*",".*pid.*","charge",".*Status.*",".*NPoi.*",".*dEdx.*","pull.*",".*X0.*",".*sigma.*",
                                "densScaling","isSecondary","hasDecay","LArm.*"],
                               ["part.*Para.*","geom.*","part.*",".*InRot.*" ],[".*"],[".*AliExternal.*","Long64.*","Long.*"], verbose=verbosity)


    if doTest:
        rdfTest=rdf1.Range(0,1000)
        rdfTest.Snapshot("testVarRDF","testVarRDF.root", varList)
        array = ak.from_rdataframe(rdfTest, columns=varList)
        df=ak.to_dataframe(array)
        print(df.head(5),df.shape)
        print(ROOT.testPullsSnapshot("testVarRDF","In"))
        print(ROOT.testPullsSnapshot("testVarRDF","Out"))
        print(ROOT.testPullsSnapshot("testVarRDF","Refit"))
    #
    else:
        #import awkward._v2 as ak
        t0=time.perf_counter()
        m0=python_process.memory_info()[0]/(2.**30)
        gc.collect()
        array = ak.from_rdataframe(rdf1, columns=varList)
        gc.collect()
        memoryUse = python_process.memory_info()[0]/(2.**30)  # memory use in GB...I think
        print('Resource rdf->awkward', m0, memoryUse-m0,time.perf_counter()-t0)
        df=ak.to_dataframe(array)
        gc.collect()
        print('Resource awkward->df', m0, memoryUse-m0,time.perf_counter()-t0)
    del tree
    gc.collect()
    return df


def testAK(varList,rdf1):
    import os
    import psutil
    import time
    import gc
    pid = os.getpid()
    python_process = psutil.Process(pid)
    #
    gc.collect()
    t0=time.perf_counter()
    m0=python_process.memory_info()[0]/2.**30
    memoryUse={}
    memoryUse["begin"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    print('memory use:', memoryUse["begin"]-m0,time.perf_counter()-t0)

    array = ak.from_rdataframe(rdf1, columns=varList)
    memoryUse["readAll"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    print('memory use read all:', memoryUse["readAll"]-m0,time.perf_counter()-t0)
    gc.collect()
    memoryUse["readAllColl"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    print('memory use read all:', memoryUse["readAllColl"]-m0,time.perf_counter()-t0)

    if 0:
    #
        #tree.SetCacheSize(100000)
        gc.collect()
        m0=python_process.memory_info()[0]/2.**30
        for i in range(0, 10):
            rdfTest = rdf1.Range(0, 4000)
            array = ak.from_rdataframe(rdfTest, columns=varList)
            gc.collect()
            memoryUse[f"array{i}"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
            print(f'memory use with gc - {i}', memoryUse[f"array{i}"]-m0, time.perf_counter()-t0)
        gc.collect()
        #
        gc.collect()
        m0=python_process.memory_info()[0]/2.**30
        for i in range(0, 10):
            rdfTest = rdf1.Range(0, 4000)
            array = ak.from_rdataframe(rdfTest, columns=varList)
            memoryUse[f"array{i}"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
            print(f'memory use without gc - {i}', memoryUse[f"array{i}"]-m0, time.perf_counter()-t0)
        gc.collect()
        #
        for i in [1, 3, 6, 9]:
            rdfTest = rdf1.Range(0, 4000*i)
            array = ak.from_rdataframe(rdfTest, columns=varList)
            gc.collect()
            memoryUse[f"array{i}"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
            print(f'memory iterative - {i}', memoryUse[f"array{i}"]-m0, time.perf_counter()-t0)
        gc.collect()

    ak.to_dataframe(array[:100])

    # t0=time.perf_counter()
    # #
    dfInner=ak.to_dataframe(array, how="inner")
    memoryUse["afterInner"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    print('memory use afterInner:', memoryUse["afterInner"],dfInner.shape,time.perf_counter()-t0)
    # dfOuter=ak.to_dataframe(array,how="outer")
    # memoryUse["afterOuter"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    # print('memory use afterOuter:', memoryUse["afterOuter"],dfOuter.shape,time.perf_counter()-t0)
    # dfInner2=ak.to_dataframe(array,how="inner")
    # memoryUse["afterInner2"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    # print('memory use afterInner2:', memoryUse["afterInner2"],dfInner2.shape,time.perf_counter()-t0)
    # dfOuter2=ak.to_dataframe(array,how="outer")
    # memoryUse["afterOuter2"] = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    # print('memory use afterOuter2:', memoryUse["afterOuter2"],dfOuter2.shape,time.perf_counter()-t0)
    # #

    return array
