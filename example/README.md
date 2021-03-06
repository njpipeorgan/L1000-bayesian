
# Peak deconvolution examples
Here we give two examples for using the Bayesian peak deconvolution algorithm. 

## Background
L1000 uses the Luminex FlexMap 3D platform, which can identify 500 different bead color as tags for different genes. To measure all 978 landmark genes within one scan, L1000 separately coupled two gene barcodes to aliquots of the same bead color and mixed them with a ratio of 2:1. In consequence, two peaks in the distribution of fluorescent intensity are expected, and a deconvolution step is involved here to access the expression of a certain gene. 

## Implementation

The algorithm has [implementations](/src) in two languages &mdash; C++ and CUDA. The implementation in CUDA is typically more efficient, so it is preferred if you want to test the algorithm on a large scale.

The examples will illustrate how to get a *z*-score from raw fluorescent intensity values, which includes

* Generate log2-FI values from a pre-defined distribution that contains two peaks;

* Setup parameters for likelihood calculation;

* Calculate a two-dimensional log-likelihood function for the locations of two peaks;

* Calculate the marginal distributions of the locations;

* Adopt an idealized reference gene expression distribution to get a *z*-score.

Scaling and normalization of the marginal distributions depend on the reads from other beads/wells, so they are not included in the examples.

## Usage

The source code can be compiled into a shared library and linked to various languages, as well as directly called and compiled with other C++ code.

### Called in C++

Enter the directory `C++` and compile the code by

    g++ -o dpeak example_dpeak.cpp ../../src/dpeak.cpp

The program will print the *z*-score of the high-abundance peak in the mock data. See the comments in `example_dpeak.cpp` for details.

### Linked to *Mathematica*

First, find the directory where you installed *Mathematica*. You can find it by evaluate `$InstallationDirectory` in *Mathematica*. We refer to it as `$InstallationDirectory` below.

Next, enter the directory `Mathematica` and compile the code into a shared library by

     g++ -I$InstallationDirectory/SystemFiles/IncludeFiles/C -fPIC --shared -o dpeak.so librarylink.cpp ../../src/dpeak.cpp

Alternatively, you can use the CUDA implementation by

    nvcc -I$InstallationDirectory/SystemFiles/IncludeFiles/C --shared -Xcompiler "-fPIC" -o dpeak.so librarylink.cpp ../../src/dpeak.cu

Finally, open `dpeak.nb` and evaluate the cells in order. `dpeak.pdf` shows an example of the results from a run.
