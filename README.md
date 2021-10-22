# fastMCKalman
Fast simulation and performance parameterization - primary for ALICE3 and DUNE  


# Visualization using RootInteractive
* RootInteractive activelly developped
  * compression on the client->server
  * interactive histogram
  * interactive histograms with weights
* master can be used for the moment
  * Release to be announced  

# Containers usage
Code will be installed using containers
* to be pulled from repository or to build from def file
  * using wget faster than pull - orders of magnitude
```
#singularity pull https://rootinteractive.web.cern.ch/RootInteractive/Singularity/alidockSingularity6.sif
wget   https://rootinteractive.web.cern.ch/RootInteractive/Singularity/alidockSingularity6.sif .
singularity shell alidockSingularity6.sif
#
wget   https://rootinteractive.web.cern.ch/RootInteractive/Singularity/tensorflow_latest-gpuCUDA.sif
singularity shell tensorflow_latest-gpuCUDA.sif
```

#Environment variables:
Assume we have environment variable with path to the software  e.g.: 
```
export fastMCKalman=$HOME/github/fastMCKalman
export fastMCKalmanData=$HOME/github/fastMCKalman
```

