Module for Ab Initio Structure Evolution (MAISE) features 
<br /> * neural network-based description of interatomic interactions
<br /> * evolutionary optimization
<br /> * structure analysis
<br />
<br /> 1. [General info](#general-info)
<br /> 2. [Download and Installation](#download)
<br /> 3. [Input](#input)
<br /> 4. [Examples](#examples)
<br /> 5. [Setup input tag description](#setup-input-tag-description)

MAISE has been developed by

<bra /> Alexey Kolmogorov <kolmogorov@binghamton.edu>  
<bra /> Samad Hajinazar <hajinazar@binghamton.edu>  
<bra /> Ernesto Sandoval <esandov1@binghamton.edu>  

---
## General info

Current version 2.7 works on Linux platforms and combines 3 modules for modeling, optimizing, and analyzing atomic structures.

1 The neural network (NN) module builds, tests, and uses NN models to
describe interatomic interactions with near-ab initio accuracy at a
low computational cost compared to density functional theory
calculations.

With the primary goal of using NN models to accelerate structure
search, the main function of the module is to relax given
structures. To simplify the NN application and comparison, we closely
matched the input and output file formats with those used in the VASP
software. Previously parameterized NN models available in the
['models/'](https://github.com/maise-guide/maise/tree/master/models)
directory have been generated and extensively tested for crystalline
and/or nanostructured materials. First practical applications of NNs
include the prediction of new synthesizable Mg-Ca alloys [1] and
identification of more stable Cu-Pd-Ag nanoparticles [2].

Users can create their own NN models with MAISE which are typically
trained on density functional theory (DFT) total energy and atomic
force data for relatively small structures. The generation of relevant
and diverse configurations is done separately with an 'evolutionary
sampling' protocol detailed in our published work [3]. The code
introduces a unique feature, 'stratified training', of how to build
robust NNs for chemical systems with several elements [3]. NN models
are developed in a hierarchical fashion, first for elements, then for
binaries, and so on, which enables generation of reusable libraries
for extended blocks in the periodic table.

2 The implemented evolutionary algorithm (EA) enables an efficient
identification of ground state configurations at a given chemical
composition. Our studies have shown that the EA is particularly
advantageous in dealing with large structures when no experimental
structural input is available [3,4].

The searches can be performed for 3D bulk crystals, 2D films, and 0D
nanoparticles. Population of structures can be generated either
randomly or predefined based on prior information. Essential
operations are 'crossover', when a new configuration is created based
on two parent structures in the previous generation, and 'mutation',
when a parent structure is randomly distorted. For 0D nanoparticles
we have introduced a multitribe evolutionary algorithm that allows an
efficient simultaneous optimization of clusters in a specified size
range [2].

3 The analysis functions include the comparison of structures based on
the radial distribution function (RDF), the determination of the space
group and the Wyckoff positions with an external SPGLIB package,
etc. In particular, the RDF-based structure dot product is essential
for eliminating duplicate structures in EA searches and selecting
different configurations in the pool of found low-energy structures.
<br />
<br /> [1] https://pubs.rsc.org/en/content/articlelanding/2018/cp/c8cp05314f#!divAbstract
<br /> [2] https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp00837c#!divAbstract
<br /> [3] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.014114
<br /> [4] https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.075501
<br /> [5] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.085131

---
## Download

The source code for MAISE can be obtained from the commandline by running:

```
git clone git://github.com/maise-guide/maise.git
```

or 

```
git clone https://github.com/maise-guide/maise.git
```

or

```
wget -O master.zip https://github.com/maise-guide/maise/archive/master.zip
unzip master.zip
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

4 A 'check' script is available in the './test/' directory which
can be run after compiling the maise executable to ensure the proper
functionality of the code. This script automatically checks for the
performance of the code in parsing the data, training the neural
network, and evaluating a crystal structure. If the compilation is
fine the 'check' script will output so; otherwise error logs
will be provided with further information about the issue.

The code has been extensively tested on Linux platforms. We will
appreciate users' feedback on the installation and performance of the
package on different platforms.

---
## Input

Main input files that define a simulation are 'setup' with job
settings, 'model' with NN parameters, and 'POSCAR' with atomic
structure parameters in the VASP format. Conversion of atomic
environments into NN inputs during the parsing stage of NN development
requires a 'basis' file that specifies Behler-Parrinello symmetry
functions.

<table>
  <tr>
    <td></td>
    <td align="center" colspan="2">EVOS</td>
    <td align="center" colspan="4">NNET</td>
    <td align="center"> CELL</td>
  </tr>
  <tr>
    <td align="center"></td>      <td align="center">SEARCH</td> <td align="center">EXAM</td> <td align="center">PARSE</td> <td align="center">TRAIN</td> <td align="center">TEST</td> <td align="center">SIMUL</td> <td align="center">EXAM</td>
  </tr>
  <tr>
    <td
    align="center"><a href="https://github.com/maise-guide/maise/blob/master/bin/setup">setup</td> <td align="center">+</td> <td align="center">+</td> <td align="center">+</td> <td align="center">+</td> <td align="center">+</td> <td align="center">+</td> <td align="center"> </td> 
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/maise-guide/maise/blob/master/bin/model">model</td> <td align="center"> </td> <td align="center"> </td> <td align="center"> </td> <td align="center">+*</td> <td align="center">+#</td> <td align="center">+#</td> <td align="center"> </td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/maise-guide/maise/blob/master/bin/basis">basis</td> <td align="center"> </td> <td align="center"> </td> <td align="center">+</td> <td align="center">$</td> <td align="center"> </td> <td align="center"> </td> <td align="center"> </td>
  </tr>
  <tr>
    <td colspan="8" class="divider"><hr /></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/atztogo/spglib">SPG</td>   <td align="center"> </td> <td align="center">+</td> <td align="center"> </td> <td align="center"> </td> <td align="center"> </td> <td align="center"> </td> <td align="center">+</td>
  </tr>
  <tr>
    <td align="center"><a href="https://www.gnu.org/software/gsl/">GSL</td>   <td align="center"> </td> <td align="center"> </td> <td align="center"> </td> <td align="center">+</td> <td align="center"> </td> <td align="center">+</td> <td align="center"> </td>
  </tr>
  <tr>
    <td align="left" colspan="8">      *   for stratified training one needs to provide individual models
                                 <br />$  'basis' stored in the parsed directory is appended to 'model' at the end of the training
                                 <br />#  'model' has 'basis' pasted at the end once training is finished</td>
  </tr>
</table>

--- 
The structure examination and manipulation functions are run by calling maise with a flag:

```
maise -flag
```

|Flag| Flag Description|
|:---:|:-|
| man|    output the list of available flags                              |
| rdf|    compute and plot the RDF for POSCAR                             |
| cxc|    compute dot product for POSCAR0 and POSCAR1 using RDF           |
| cmp|    compare RDF, space group, and volume of POSCAR0 and POSCAR1     |
| spg|    convert POSCAR into str.cif, CONV, PRIM                         |
| cif|    convert str.cif into CONV and PRIM                              |
| rot|    rotate  a nanoparticle along eigenvectors of moments of inertia |
| dim|    find    whether POSCAR is periodic (3) or non-periodic (0)      |
| box|    reset   the box size for nanoparticles                          |
| sup|    make    a supercell specified by na x nb x nc                   |
| vol|    compute volume per atom for crystal or nano structures          |

## Examples

Directory 'examples/' has samples of maise jobs for unit cell analysis
and manipulation, parsing data, training neural networks, simulating
structures with neural network models, evolutionary search for ground
state wit neural network model, molecular dynamics run, and phonon
calculation. Eash example has a README file, a setup file with only
relevant tags for the particular job, and reference output files for
comparison.

## Input tags by type
---
[Main job type selector](#main-job-type-selector)\
[Structure-enviroment](#structure-enviroment)\
[Main EVOS](#main-evos)\
[EVOS operations](#evos-operations)\
[EVOS crossover/mutation](#evos-crossover/mutation)\
[Molecular dynamics](#molecular-dynamics)\
[Species related](#species-related)\
[I/O](#i/o)\
[General model](#general-model)\
[Neural Network model](#neural-network-model)\
[Neural Network training](#neural-network-training)\
[Parsing](#parsing)\
[Cell relaxation](#cell-relaxation)

### Main job type selector
[JOBT](#jobt)
### Structure enviroment
[NMAX](#nmax) [MMAX](#mmax)
### Main EVOS         
[CODE](#code) [DENE](#dene) [KMSH](#kmsh) [LBOX](#lbox) [NDIM](#ndim) [NITR](#nitr) [NNJB](#nnjb) [NPOP](#npop) [RAND](#rand) [RUNT](#runt) [SEED](#seed) [SITR](#sitr) [TINI](#tini)  
### EVOS operations      
[BLOB](#blob) [CHOP](#chop) [INVS](#invs) [MATE](#mate) [MUTE](#mute) [PACK](#pack) [PLNT](#plnt) [REFL](#refl) [RUBE](#rube) [SWAP](#swap) [TETR](#tetr)   
### EVOS crossover/mutation 
[ACRS](#acrs) [ADST](#adst) [ELPS](#elps) [LCRS](#lcrs) [LDST](#ldst) [MCRS](#mcrs) [SCRS](#scrs) [SDST](#sdst)   
### Molecular dynamics             
[CPLT](#cplt) [CPLP](#cplp) [ICMP](#icmp) [DELT](#delt) [MOVI](#movi) [NSTP](#nstp) [MDTP](#mdtp) [TMAX](#tmax) [TMIN](#tmin) [TSTP](#tstp)     
### Species related 
[ASPC](#aspc) [NSPC](#nspc) [TSPC](#tspc) 
### I/O 
[COUT](#cout) [DATA](#data) [DEPO](#depo) [EVAL](#eval) [OTPT](#otpt) [WDIR](#wdir)
### Neural Network model
[NCMP](#ncmp) [NNGT](#nngt) [NNNN](#nnnn) [NNNU](#nnnu) [NSYM](#nsym) 
### Neural Network training    
[FMRK](#fmrk) [LREG](#lreg) [NTRN](#ntrn) [NTST](#ntst) [TEFS](#tefs) [NPAR](#npar)
## Parsing
[EMAX](#emax) [FMAX](#fmax) [FMIN](#fmin) [VMAX](#vmax) [VMIN](#vmin) [MMAX](#mmax) 
### Cell Relaxation
[ETOL](#etol) [MINT](#mint) [MITR](#mitr) [PGPA](#pgpa) [RLXT](#rlxt) [TIME](#time)

---
## Setup input tag description
---
| TAG | DESCRIPTION |
|:--|:---------|
| <a name="jobt"></a>JOBT | structure analysis (00) use analysis tools specified by flags, evolutionary search (10) run (11) soft exit (12) hard exit (13) analysis, cell simulation (20) relaxation (21) molecular dynamics (22) phonon calculations, data parsing (30) prepare inputs for NN training , NN training (40) full training (41) stratified training|
| <a name="code"></a>CODE | Type of the code in use. (0) MAISE-INT (1) VASP-EXT (2) MAISE-EXT|
| <a name="npar"></a>NPAR | Number of cores for parallel NN training or cell simulation|
| <a name="mint"></a>MINT | The optimizer algorithm for the neural network training and the cell optimization. (gsl minimizer type (0) BFGS2 (1) CG-FR (2) CG-PR (3) steepest descent |
| <a name="mitr"></a>MITR | Maximum number of the optimization steps; if the desired accuracy is not reached for NN training or cell optimization steps|
| <a name="rlxt"></a>RLXT | Cell optimization type (2) force only (3) full cell (7) volume (ISIF in VASP)|
| <a name="etol"></a>ETOL | Error tolerance for training or cell optimization convergence|
| <a name="tefs"></a>TEFS | Training target value (0) E (1) EF (2) ES (3) EFS (4) TOY|
| <a name="fmrk"></a>FMRK | Fraction of atoms that will be parsed to use for EF or EFS trainings|
| <a name="cout"></a>COUT | Output type in the OUTCAR file in cell evaluation and optimization|
| <a name="nmax"></a>NMAX | Maximum number of atoms in the unit cell|
| <a name="mmax"></a>MMAX | Maximum number of neighbors within the cutoff radius|
| <a name="nspc"></a>NSPC | Number of element types for evolutionary search, parsing the data and neural network training.|
| <a name="tspc"></a>TSPC | Atomic number of the elements specified with NSPC tag|
| <a name="aspc"></a>ASPC | Number of atoms of each element for the evolutionary search|
| <a name="nsym"></a>NSYM | Number of the Behler-Parrinello symmetry functions for parsing data using the "basis" file|
| <a name="ncmp"></a>NCMP | The length of the input vector of the neural network|
| <a name="ntrn"></a>NTRN | Number of structures used for neural network trainin (negative number means percentage)|
| <a name="ntst"></a>NTST | Number of structures used for neural network testing (negative number means percentage)|
| <a name="nnnn"></a>NNNN | Number of hidden layers in the neural network (does not include input vector and output neuron)|
| <a name="nnnu"></a>NNNU | Number of neurons in hidden layers|
| <a name="nngt"></a>NNGT | Activation function type for the hidden layers' neurons (0) linear (1) tanh|
| <a name="emax"></a>EMAX | Parse only this fraction of lowest-energy structures. From 0 to 1|
| <a name="fmax"></a>FMAX | Will not parse data with forces larger than this value|
| <a name="vmin"></a>VMIN | Will not parse data with volume/atom smaller than this value|
| <a name="vmax"></a>VMAX | Will not parse data with volume/atom  larger than this value|
| <a name="ndim"></a>NDIM | Dimensionality of the unit cell in evolutionary search and cell optimization (3) crystal (2) film (0) particle|
| <a name="lbox"></a>LBOX | Box dimension for generating particles in evolutionary search in Angs (ignored for crystals)|
| <a name="npop"></a>NPOP | Population size in the evolutionary search|
| <a name="sitr"></a>SITR | Starting iteration in the evolutionary search (0) start from random or specified structures|
| <a name="nitr"></a>NITR | Number of iterations in the evolutionary search (should be larger than SITR)|
| <a name="tini"></a>TINI | Type of starting the evolutionary search when SITR=0|
| <a name="time"></a>TIME | Maximum time for cell relaxation in evolutionary search and cell optimization|
| <a name="pgpa"></a>PGPA | Pressure in GPa|
| <a name="dene"></a>DENE | Store distinct structures generated in evolutionary search in POOL/ if within this energy/atom (eV/atom) window from the ground state|
| <a name="kmsh"></a>KMSH | K-mesh density used for VASP-EVOS. Suggested values: 0.30 for s/c, 0.05 for metals|
| <a name="seed"></a>SEED | Starting seed for the random number generator in evolutionary search (0) Uses time as seed (+) The seed value|
| <a name="rand"></a>RAND | Starting seed for the parsing of the dataset. (0) Uses time as seed (+) The seed value (-) No randomization: structures are parsed in listing order|
| <a name="tmin"></a>TMIN | Minimum temperature in MD runs (K)|
| <a name="tmax"></a>TMAX | Maximum temperature in MD runs (K)|
| <a name="tstp"></a>TSTP | Temperature step in MD runs (K) in running form TMIN to TMAX|
| <a name="delt"></a>DELT | Time step in the MD runs|
| <a name="nstp"></a>NSTP | Number of steps per temperature in MD runs|
| <a name="cplt"></a>CPLT | Coupling constant in Nose-Hoover thermostat for MD runs. Suggested: 25.0|
| <a name="cplp"></a>CPLP | Coupling constant in Brendsen barostat for MD runs. Suggested: 100.0|
| <a name="icmp"></a>ICMP | Isothermal compressibility in Brendsen barostat for MD runs (in 1/GPa)|
| <a name="movi"></a>MOVI | Number of steps after which a snapshot of structure will be saved during the MD run|
| <a name="mdtp"></a>MDTP | MD run type (10) NVE (20) NVT: Nose-Hoover (30) NPT: Nose-Hoover and Brendsen (40) Isobaric (11,21,31,41) runs with velocisities read in from POSCAR file|
| <a name="depo"></a>DEPO | Path to the DFT datasets to be parsed|
| <a name="data"></a>DATA | Location of the parsed data to parse or read for training (will be overwritten during parsing)|
| <a name="otpt"></a>OTPT | Directory for storing model parameters in the training process|
| <a name="eval"></a>EVAL | Directory for model testing data|
| <a name="wdir"></a>WDIR | Work directory for evolutionary search, MD runs, etc.|
| <a name="tetr"></a>TETR | Fraction of the structures generated randomly using tetris operation. From 0 to 1|
| <a name="plnt"></a>PLNT | Fraction of the structures generated from seeds. From 0 to 1|
| <a name="pack"></a>PACK | Fraction of the structures generated from closed-pack structures. From 0 to 1|
| <a name="blob"></a>BLOB | Fraction of the structures generated randomly using blob shape. From 0 to 1|
| <a name="mate"></a>MATE | Fraction of the structures generated by crossover using two halves from each parent. From 0 to 1|
| <a name="swap"></a>SWAP | Fraction of the structures generated by crossover using core and shell of parents. From 0 to 1|
| <a name="rube"></a>RUBE | Fraction of the structures generated by Rubik's cube operation. From 0 to 1|
| <a name="refl"></a>REFL | Fraction of the structures generated by symmetrization via reflection. From 0 to 1|
| <a name="invs"></a>INVS | Fraction of the structures generated by symmetrization via inversion. From 0 to 1|
| <a name="chop"></a>CHOP | Fraction of the structures generated by chopping to make facets. From 0 to 1|
| <a name="mute"></a>MUTE | Fraction of the structures generated by random distortions to the structure. From 0 to 1|
| <a name="mcrs"></a>MCRS |  0.50              mutation rate in crossover|
| <a name="scrs"></a>SCRS |  0.00              crossover:  swapping rate|
| <a name="lcrs"></a>LCRS |  0.00              crossover:  mutation strength for lattice vectors|
| <a name="acrs"></a>ACRS |  0.10              crossover:  mutation strength for atomic positions|
| <a name="sdst"></a>SDST |  0.00              distortion: swapping rate|
| <a name="ldst"></a>LDST |  0.00              distortion: mutation strength for lattice vectors|
| <a name="adst"></a>ADST |  0.20              distortion: mutation strength for atomic positions|
| <a name="elps"></a>ELPS |  0.30              random:     nanoparticle ellipticity|
