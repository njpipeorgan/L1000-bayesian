# L1000 Bayesian pipeline

## Pipeline overview

The pipeline takes raw fluorescence intensity data from  LINCS L1000 Phase II datasets as input and gives a combined *z*-score profile for each experiment as its signature. The pipeline is composed the following steps:

1. LISS and quality control. In this step, we perform a two step linear scaling to calibrate the fluorescent intensities. The sameple data we give are of good quality so there is no code to remove bad wells, but you can add different types of quality control on other data based on calibration parameters or other properties when necessary.

2. Peak deconvolution. For beads coupled with two different transcript barcodes, a deconvolution step is involved to infer the peak position for each gene. Two probability distributions will be given to the transcripts as the estimations of their expression levels.

3. Quantile normalization. The shape of expression profile is standardized across all samples on the same plate so that different samples on the same plate are comparable to each other.

4. *z*-score inference. *z*-scores are inferred from the probability distribution for each gene to represent relative gene expression.

5. Combining replicates. The sample data only include the measurements on one of the three bio-replicates, so there is no combination process in the code. If there are replicates in your dataset, you can combine the *z*-scores you get from this pipeline.

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
