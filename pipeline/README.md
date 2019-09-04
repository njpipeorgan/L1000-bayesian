# L1000 Bayesian pipeline

## Prerequisites

To run the pipeline, you should have [*Wolfram Mathematica*](http://www.wolfram.com/mathematica/) or the free [*Wolfram Engine*](https://www.wolfram.com/engine/) installed. Run the following command to validate the installation: 

    wolframscript -c 2+2

If you want to use the CUDA implementation for peak deconvolution, you should install [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) as well. 

## Build shared library

A shared library should build from the peak deconvolution source code to run the pipeline. To build the library, follow the instruction in the [example for *Mathematica*](https://github.com/njpipeorgan/L1000-bayesian/tree/master/example#linked-to-mathematica). The pipeline assumes the library to be put in the current directory and named `dpeak.so`. 

## Run the pipeline

We prepared a small portion of LXB files from `REP.A028_MCF7` in `sample_lxb` directory. To run the pipeline on the sample data, execute

    wolframscript -file pipeline.wl sample_lxb

If it succeed, there will be three files in the `output` directory, storing the LISS parameters, marginal distributions, and *z*-scores. Note that the *z*-scores are different from those in the database, because the sample data is only cover a part of the plate. 
