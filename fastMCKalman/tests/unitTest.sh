
# source $fastMCKalman/fastMCKalman/tests/unitTest.sh


alias helpCat=cat
[[ -x "$(command -v pygmentize)" ]] && alias helpCat="pygmentize -O style=borland,linenos=1 -l bash"
init(){
  cat <<HELP_USAGE | helpCat
  makeData
  makePullTest
  makePullTestSeed
  analyzeLogs
HELP_USAGE
}

makeData(){
     export nPoints=${1:-40000}

    cat <<EOF >  makeData.sh
#!/bin/bash
    root.exe -n -b -l <<\EOF 2>&1 | tee makeData.log
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
   .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
    AliPDG::AddParticlesToPdgDataBase();
    AliLog::SetPrintRepetitions(0);
    testTPC(${nPoints},kTRUE);            //setup for the looper development
    .q
EOF
   chmod a+x makeData.sh

   ./makeData.sh
}

makePullTest(){
       cat <<EOF >  makePullTest.sh
#!/bin/bash
    echo $(pwd)/fastParticle.root >fastParticle.list
    root.exe -n -b -l <<\EOF 2>&1 | tee makePullTest.log
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
     .L $fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g
     initTreeFast()
     testPulls("","In","&&Iteration$==0")
     testPulls("Full","In","&&Iteration$==0")
     testPulls("Full","Out","&&Iteration$==30")
     testPulls("Full","Refit","&&Iteration$==30")
    .q
EOF
   chmod a+x makePullTest.sh
   ./makePullTest.sh
}

makePullTestSeed(){
       cat <<EOF >  makePullTestSeed.sh
#!/bin/bash
    echo $(pwd)/fastParticle.root >fastParticle.list
    root.exe -n -b -l <<\EOF 2>&1 | tee makePullTestSeed.log
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
     .L $fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g
     initTreeFast()
     testPullsSeed()
    .q
EOF
   chmod a+x makePullTestSeed.sh
   ./makePullTestSeed.sh
}

analyzeLogs(){
   errors=( "short track" "Too few consecutive points" "Rotation failed" "Propagation failed" "Update failed" "Too big chi2"
   "Correct for material failed" "PropagateToMirrorX failed" )
   errorSources=( "fastParticle::reconstructParticleFull:" "fastParticle::reconstructParticle:" "fastParticle::reconstructParticleFullOut:")
   for errorSource in  "${errorSources[@]}"; do
      nErrorstot=$(cat makeData.log | grep -c ${errorSource})
      echo  ${errorSource} ${nErrorstot}
      for error in "${errors[@]}"; do
        nErrors=$(cat makeData.log | grep ${errorSource} | grep -c "${error}")
        echo ${errorSource} ${error} ${nErrors}
      done;
   done
}

logDiff(){
  # TODO
  diff --side-by-side -W 200  makePullTestSeed_v1.log makePullTestSeed.log |sed  s/"testFastTracker seed pull"//g

}

init