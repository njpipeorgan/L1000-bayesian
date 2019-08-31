# L1000 peak deconvolution based on Bayesian analysis

## Overview


## Catalogs

### Downloads

LINCS L1000 Phase II (GSE 70138) catalogs generated by our pipeline are currectly available. 

| Description                               | Download                                      |
| ----------------------------------------- | --------------------------------------------- |
| Marginal distributions of peak locations  | [Bayesian_GSE70138_Level2_DPEAK.zip](http://callisto.astro.columbia.edu/files/L1000/Bayesian_GSE70138_Level2_DPEAK.zip)|
| Plate control *z*-scores                  | [Bayesian_GSE70138_Level4_ZSPC_n335465x978.h5](http://callisto.astro.columbia.edu/files/L1000/Bayesian_GSE70138_Level4_ZSPC_n335465x978.h5)|
| Combined *z*-scores by bio-replicates     | [Bayesian_GSE70138_Level5_COMPZ_n116218x978.h5](http://callisto.astro.columbia.edu/files/L1000/Bayesian_GSE70138_Level5_COMPZ_n116218x978.h5)|
| Checksum                                  | [Bayesian_GSE70138_sha512sum.txt](http://callisto.astro.columbia.edu/files/L1000/Bayesian_GSE70138_sha512sum.txt)|

The meta data are available from the publication by L1000 group: [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138). They include perturbagen and cell line information associated with signature and well IDs in the catalogs. 

### Data stuctures

The *z*-score catalogs (as HDF5) are compatible with the original catalogs by L1000 group. Each of them contains three datasets as follows:
* `/colid` are the signature IDs (Level 5) or well IDs (Level 4); 
* `/rowid` are the names of landmark genes; 
* `/data` are the *z*-scores as a matrix. 

Each marginal distribution catalog contain the information of peak locations on one plate. It contains four datasets as follows:
* `/colid` are the well IDs; 
* `/rowid` are the names of landmark genes; 
* `/peakloc` are the loations of the peaks for calculating likelihood function; 
* `/data` are encoded log-likelihoods as a rank-3 array of 16-bit unsigned integers. To retrive the log-likelihoods, the values should be multiplied by a factor of `-0.001`.

## Installation

The Bayesian peak deconvolution algorithm has implementations in two languages &mdash; C++ and CUDA. They can be directly call from other C++ code or compiled and linked to external code. Note that the implementation in CUDA is typically more efficient, so it is preferred if you want to test the algorithm on a large scale. 

Here we give an example of compiling them to a shared library that can be called from *Wolfram Mathematica*. The environment variable `$InstallationDirectory` below is the directory where you installed *Mathematica*. 

* Using the implementation in C++

      g++ -I$InstallationDirectory/SystemFiles/IncludeFiles/C -fPIC --shared -o dpeak.so librarylink.cpp dpeak.cpp
    
* Using the implementation in CUDA

      nvcc -I$InstallationDirectory/SystemFiles/IncludeFiles/C --shared -Xcompiler "-fPIC" -o dpeak.so librarylink.cpp dpeak.cu

See [example/dpeak.nb](example/dpeak.nb) for how to load a *LibraryLink* function from the `dpeak.so`. 


## Citation