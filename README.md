# Saturn-Mag-Model
FORTRAN source code for a Saturnian magnetospheric empirical magnetic field model derived from Cassini magnetometer data 

<img src="docs/3dcurrents.png" width="300">

# About
This is a FORTRAN 77ish implementation of an external magnetic field
model for Saturn's magnetosphere developed by G.K.Stephens.

The inputs to the model are the solar wind dynamic pressure in nanoPascals
(nPa), the Sun's latitude in the KSMAG or IAU-Saturn (AKA Saturn SIII) coordinate
system in radians, and the position in KSM coordinates units of Saturn radii (Rs).
The output is the external magnetic field vector in KSM coordinates in units of
nanoTesla (nT). Subroutines are provided that can be used to look-up/compute the
dynamic pressure and Sun latitude as a function of time. Additionally, the
Dougherty et al. (2018) internal magnetic field model is provided.

We recommend looking at example_saturn.for as a comprehensive illustration of
how to use this API.

The code is written and developed using GNU's Fortran Compiler, or gfortran. It
may compile and run with other Fortran compilers as well. The code is said to be
FORTRAN '77ish' because we do not strictly adhere to the 77 standard and may use
functionality from later versions of Fortran. For instance, we DO NOT adhere to
capital letters and variable and subroutine names of 6 or less characters.

The authors also maintain a Java implementation of this code base and the model
is fit using the Java implementation. Agreement between the Fortran and Java
codes is maintained via a series of test codes that comprehensively loops over
model inputs and compares the model outputs. Agreement is good to about 10
digits.

# Getting Started
The model source code is self-contained, meaning it requires no external dependencies or libraries. The model was developed on a Mac using [GFortran](https://gcc.gnu.org/wiki/GFortran).

Assuming GFortran is [installed](https://gcc.gnu.org/wiki/GFortranBinaries), the model can be compiled using the following: 
`gfortran example_saturn.for saturn_int.f saturn_ext.f cart_field_to_sphere.f dipoleshield.f igrf_geo_08.f read_pdyn.f tailsheet.f tailsheet_sym.f tailsheet_asym.f bowldeform.f cartharmonic.f cartharmonic_alt.f igrf_geo_cart.f rotate_about_y.f bessjj.f utc_to_tdt.f cart_field_to_cyl.f find_left_index.f ksun.f bowldist.f -o example_saturn`

This will compile the example program contained in the file `example_saturn.for` into an executable file named `example_saturn`.

The executable can then be run on Unix systems using `./example_saturn`. The output of this program should produce:

>Date:            2013-04-21 22:53:28.39823\
 TDT:                   419856875.58223003     \
 Dp:                    2.3699125945428222E-003\
 sunLat:               0.31991904766201745     \
 (x,y,z):              -12.381273999999999        3.8901874300000001        2.1981273432987001     \
 inside:              T\
 (bx,by,bz):            6.2702298378116481       -1.2576152407132852        7.0949159196898739     \
 (hx,hy,hz):           0.74732969552778661       -1.1478438385940311       -9.4198143792161204     \
 (r,theta,phi):         13.162874607051583        1.4030157792196860        2.8371604325738131     \
 (br,btheta,bphi):     -7.7000229788156362        1.0538352925036358       0.19132520816048570     \
 Closest CS (x,y,z):   -11.952407165927681        3.9131807515159531        3.7265080120331575     \
 Dist. to CS:           1.5875776903077010     

# FAQs
None as of yet.

# License
[3-Clause BSD](LICENSE)

# Author
Grant K. Stephens

Grant.Stephens@jhuapl.edu

[Webpage](https://civspace.jhuapl.edu/people/grant-stephens)

[Twitter](https://twitter.com/GrantKStephens)

# Acknowledgments
This work was funded by the NASA grant 80NSSC21K0533. The authors would like to thank all the people who helped develop and advance this model and source code. The data used to fit the model is from the Cassini Magnetometer (Dougherty et al. 2004) and was obtained from the [Planetary Data System (PDS)](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/inst-mag.html). The internal magnetic field is an implementation of [Dougherty et al. (2018)](https://doi.org/10.1126/science.aat5434) and the source code is inspired from N.A.Tsyganenko's Geopack library. The model structure, particularly for the equatorial and dipole shielding fields, is derived from the TS07D model [(Tsyganenko & Sitnov, 2007)](https://doi.org/10.1029/2007JA012260) for Earth's magnetic field. The date and time API was inspired from the NASA/NAIF SPICE toolkit along with K.Khurana's KMAG/JMAG source code. The modeled magnetopause is from [Pilkington et al. (2015)](https://doi.org/10.1002/2015JA021290), additionally, the magnetopause crossings list provided by C.Jackman is used to filter out magnetometer data outside the magnetosphere. The bowl deformation derives from [Arridge et al. (2008)](https://doi.org/10.1029/2007JA012963) and the implementation is derived from K.Khurana's [KMAG model](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/files/fp/KMAG2012test.f). The solar wind dynamic pressure is derived from the [Tao et al. (2005)](https://doi.org/10.1029/2004JA010959) model and was obtained from the AMDA archive. The subroutine that computes the Sun's position (KSUN) is taken from K.Khurana's KMAG model. The Bessel function evaluator was supplied by J.M.Albert.
Please see NASA's [Cassini Magnetospheric Science](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/sci-fp.html) webpage. K.Khurana's KMAG model is also located on the [PDS](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/files/fp/KMAG2012test.f).


