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
Blank for now.

# FAQs
None as of yet.

#License
3-Clause BSD

# Author
Grant K. Stephens
             Grant.Stephens@jhuapl.edu
https://twitter.com/GrantKStephens
[webpage](https://civspace.jhuapl.edu/people/grant-stephens)

# Acknowledgments
The authors would like to thank all the people who helped
develop and advance this model and source code. The data used to fit the model
is from the Cassini Magnetometer (Dougherty et al. 2004) and was obtained from
the Planetary Data System (PDS). The internal magnetic field is an
implementation of Dougherty et al. (2018) and the source code is derived from
N.A.Tsyganenko's Geopack library. The structure and source code for the
equatorial and dipole shielding fields is derived from the works and code of
N.A.Tsyganenko and M.I.Sitnov. The time conversion API derives from the
NASA/NAIF SPICE toolkit along with K.Khurana's KMAG/JMAG source code. The
modeled magnetopause is from Pilkington et al. (2015), additionally, the
magnetopause crossings list provided by C.Jackman is used to filter out
magnetometer data outside the magnetosphere. The bowl deformation derives from
Arridge et al. (2008) and the implementation is derived from K.Khurana's KMAG
model. The solar wind dynamic pressure is derived from the Tao et al. (2005)
model and was obtained from the AMDA archive. The subroutine that computes the
Sun's position (KSUN) is taken from K.Khurana's KMAG model. The Bessel function
evaluator was supplied by J.M.Albert. Several subroutines derive from
subroutines included in the NASA/NAIF SPICE toolkit. 



