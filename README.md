# fastMCKalman
Fast simulation and performance parameterization - primary for ALICE3 and DUNE  


# Assumptions

Code will be installed using containers
* to be pulled from repository or to build from def file
  * to check web server to store container
```
singularity pull https://rootinteractive.web.cern.ch/RootInteractive/Singularity/alidockSingularity6.sif
singularity shell ....
```


Assume we have environment variable with path to the software  e.g.: 
```
export fastMCKalman=$HOME/github/fastMCKalman
```

