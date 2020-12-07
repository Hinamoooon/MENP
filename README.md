# MENP
Multipole Expansion for NanoPhotonics (MENP)  

## General
MENP is an open-source MATLAB-based package for multipole expansion based on induced current distributions in nanostructures. It imports electric field distribution obtained by full-field simulation techniques (e.g., FDTD, FEM, etc.) and computes electric and magnetic multipoles (ED, MD, EQ, MQ) based on exact expression of multipole expansion. It is also capable of multipole expansion under long-wavelength approximation to find a contribution of toroidal dipole moments (TD).

Target users are researchers in the field of nanophotonics. In particular, recently emerged sub-wavelength Mie resonators which exhibit rich spectral features due to multipole resonances. The resultant multipolar interferences is opening up new opportunities to realize novel functionarity such as unidirectinal scattering (so-called Kerker condition), non-radiating optical Anapole state, and so on. For structural design and interpretation of physics in such systems, full-field simulation accoumpanied by multipole expansion is essential.

This package is basically designed for the use with Lumerical FDTD solutions, but any software can be used by exporting four-dimensional (x,y,z,f) electric field and refractive index data as a MATLAB .mat file.

![Scattering spectra of a silicon nanosphere (R = 100 nm)](https://github.com/Hinamoooon/MENP/blob/main/SiNP_r100.png?raw=true)

## Reference
We ask you to cite the following paper when you publish results obtained with MENP.  
**Tatsuki Hinamto and Minoru Fujii "MENP: An Open-Source MATLAB Implementation of Multipole Expansion for Applications in Nanophotonics", [arXiv:2011.03661](https://arxiv.org/abs/2011.03661)**

## Supplementary references to implemented formulation
1. Alaee, R.; Rockstuhl, C.; Fernandez-Corbaton, I. An Electromagnetic Multipole Expansion beyond the Long-Wavelength Approximation. [Opt. Commun. 2018, 407, 17–21.](https://www.sciencedirect.com/science/article/pii/S003040181730754X)  
2. Baryshnikova, K. V.; Smirnova, D. A.; Luk’yanchuk, B. S.; Kivshar, Y. S. Optical Anapoles: Concepts and Applications. [Adv. Opt. Mater. 2019, 7, 1801350.](https://onlinelibrary.wiley.com/doi/full/10.1002/adom.201801350)  
3. Hasebe, H.; Sugimoto, H.; Hinamoto, T.; Fujii, M. Coupled Toroidal Dipole Modes in Silicon Nanodisk Metasurface: Polarization Independent Narrow Band Absorption and Directional Emission. [Adv. Opt. Mater. 2020, 2001148.](https://onlinelibrary.wiley.com/doi/full/10.1002/adom.202001148)

## License
MIT

## Author
Tatsuki Hinamoto@Kobe University, Japan

## How to use
Run following demo to understand the usage.  
./demo_sphere (exact and approximated multipole expansion for a silicon nanosphere)  
./demo_disk (approximated multipole expansion including toroidal dipole moment fro a silicon nanodisk).  

For the computation, three dimensional electric field distribution and refractive index data computed around a target nanostructure are required. On Lumerical FDTD Solutions, this exporting process can be done by running a lumerical script "./lumerical_script/EField2MAT.lsf". As an example, Lumerical project files (.fsp) are also included in the demo directories.

## Input file format
MENP requires electric field (E) and refractive index (n) distributions and their coordinates (x,y,z,f). To see the spectral dependence of decomposed scattering cross sections, it should have frequency axis (f) in addition to position vector (x,y,z), that is, four dimensional data of E(x,y,z,f) and n(x,y,z,f) as vectors (i.e., Ex, Ey, Ez, n_x, n_y, n_z). Each array of coordinate x,y,z,f should be an array that have a size of *len* x 1, where *len* indicates the length of each array.

The attributes are passed to MENP's main functions (e.g., `exactME.m`) by `exactME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z)`.

## Directory structure
### ./MENP
#### Main functions
- `exactME.m`: This function imports electric field and refractive index distributions and conducts multipole expansion based on an exact expression. The results are returned as total scattering cross section and partial ones from electric dipole (p), magnetic dipole (m), electric quadrupole (Qe), and magnetic quadrupole (Qm).
- `approxME.`: This function imports electric field and refractive index distributions and conducts multipole expansion under an long-wavelength approximation. The results have the same form as `exactME.m`.
- `toroidalME.m`: This function calculates the approximated multipole expansion similarly to `approxME.m` but expands into multipoles including a toroidal dipole moment (T).

#### Supporting functions
- `PhysicalConst.m`: physical constants
- `E2J.m`: conversion of electric field distributions (E) into current density distributions (J)
- `trapz4Dto1D`: Expansin of MATLAB function `trapz()`, which computes trapezoidal numerical integration, for processing an four dimensional array.
- `toroidalME_phase.m`: Supplementary function of `toroidalME.m`. This function provides relative phases of electric and toroidal dipole moments to look into the anapole condition.

### ./lumerical_script
- `EField2MAT.lsf`: Lumerical script file that exports required data (i.e., electric field and refractive index distribution and x,y,z,f coordinates) from a Lumerical FDTD simulation project file (.fsp) as an all-in-one MATLAB file named `ENxyzf.mat`. As can be seen in the script, by default, it reads the results of 3D monitors named "field" and "index" integrated into an analysis group "multipole".

### ./demo_sphere
Demo for a silicon nanosphere with a radius of 100 nm for computation of exact and approximated multipole expansion.
- `demo_exact.m`: This script shows how to import the MATLAB file which contains required data (`ENxyzf.m`) and compuate the multipole expansion with `exact.m`. By running it, total and partial scattering spectra are shown, and a csv file (`demo_exact.csv`) is exported.
- `demo_approx.m`: Sample code which shows how to use `approx.m`
- `demo_exact.csv`: Sample output data from `demo_exact.m`
- `demo_approx.csv`: Sample output data from `demo_approx.m`
- `MENP_SiNP_r100.fsp`: Sample Lumerical Project for comutation of multipole expansion for a silicon nanophere with 100 nm in radius.
- `ENxyzf.mat`: Sample input data exproted from `MENP_SiNP_r100.fsp` with `EField2MAT.lsf`.

### ./demo_disk
Demo for a silicon nanodisk for computation of multipole expansion (under long-wavelength approx.) into multipoles including a toroidal dipole moment.
- `demo_toroidal.m`: This script shows how to compuate the multipole expansion with `toroidalME.m`. In addition to scattering cross sections, it calculates phases of electric and toroidal dipole moments for applications in anapole states.
- `demo_toroidal.csv`: Sample output data of scattering cross sections from `demo_toroidal.m`
- `demo_toroidal_phase.csv`: Sample output data of phases from `demo_toroidal.m`
- `MENP_SiND_D310H50.fsp`: Sample Lumerical Project for comutation of multipole expansion for a silicon nanodisk with diameter of 310 nm and height of 50 nm.
- `ENxyzf.mat`: Sample input data exproted from `MENP_SiND_D310H50.fsp` with `EField2MAT.lsf`.