Module for Ab Initio Structure Evolution (MAISE) features


Hyperdimensional optimization version of MAISE has been developed by

<br /> Alexey Kolmogorov <kolmogorov@binghamton.edu>
<br /> Maxwell Meyers <mtm6098@psu.edu>
<br /> Daviti Gochitashvili <dgochit1@binghamton.edu>

---
## General info

This version of MAISE works on Linux platforms and was designed to optimize structures in higher (3D+) dimensionsions.
Behler–Parrinello symmetry functions map atomic environments into machine learning input of constant length. The functions can be naturally extended into higher-dimensional spaces with N coordinates by calculating distances and angles.

## Download

The source code for MAISE can be obtained from the commandline by running:

```
git clone --branch hyper.maise --single-branch https://github.com/maise-guide/maise.git hyper.maise
cd hyper.maise
```

## Installation

1 Use '**make --jobs**' for full compilation. For recompilation, use 'make clean' to remove
object files or 'make clean-all' to remove object files and external libraries.

2 During MAISE compilation, 'make --jobs' checks if two required
external libraries, [GSL library](https://www.gnu.org/software/gsl/)
and [SPGLIB v1.11.2.1, Feb 2019](https://atztogo.github.io/spglib),
are present. If not, they will be automatically downloaded to
./ext-dep and installed in ./lib on most systems.

3 If the GSL or SPGLIB installation is not completed automatically
please compile them manually and copy (i) libgsl.a, libgslcblas.a and
libsymspg.a into the './lib' subdirectory; (ii ) the 'spglib.h' header
into './lib/include' subdirectory; and (iii) all gsl headers into the
'./lib/include/gsl' subdirectory.

4 The 'TEST' directory contains all necessery files for hyperdimensional optimzation.
User can test compiled code using the following command:
cp maise TEST; cd TEST; maise



The code has been extensively tested on Linux platforms. We will
appreciate users' feedback on the installation and performance of the
package on different platforms.

---
## Input

Main input files that define a simulation are 'setup' with job
settings, 'model' with NN parameters, and 'POSCAR' with atomic
structure parameters in the VASP format. 
For the general input parameters please check the README file of 
the master branch.



## Setup input tag description
---
| TAG | DESCRIPTION |
|:--|:---------|
| <a name="hmap"></a>HMAP | Mapping of normal coordinates for hyperspatial optimization: random displacements of atoms into each extra dimension (0) normal coordinates are uniformly contracted to compensate for the increase of the interatomic distances in the hyper space (1) roughly uniform atomic distribution via a correlated adjustment of the normal and hyper coordinates (2)|
| <a name="step"></a>STEP | Number of equal groups into which the total maximum number of minimization steps is divided|
| <a name="mnot"></a>MNOT | Initial value of harmonic penalty µ|
| <a name="beta"></a>BETA | Multiplicative factor by which a control parameter is increased after each structural optimization cycle: mu_{k+1} = mu_k * beta|
| <a name="delt"></a>DELT | Strength of random displacements in extra coordinate(s)|
