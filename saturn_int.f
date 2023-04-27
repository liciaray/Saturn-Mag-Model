      SUBROUTINE saturn_int (sunAngle,
     .     x,y,z,
     .     bx,by,bz)
c
c     Computes the internal magnetic field for Saturn in the KSM
c     coordinate system in units of nT. The coefficients are from
c     Dougherty et al. (2018) (https://doi.org/10.1126/science.aat5434).
c     The implementation is derived from N.A. Tsyganenko's
c     GEOPACK-2008 library
c     (http://geo.phys.spbu.ru/~tsyganenko/modeling.html).
c      
c     Inputs:
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 x,y,z: The location to evaluate the value of the internal
c       magnetic field in the KSM coordinate system in units of Saturn
c       radii.
c      
c     Outputs:
c     real*8 bx,by,bz: The internal internal field evaluated at the
c       supplied location in KSM coordinates in units of nT.
c     
      implicit none

c     Inputs
      REAL*8 sunAngle, x,y,z

c     Outputs
      REAL*8 bx,by,bz

c     Internal variables
      REAL*8 xgeo,ygeo,zgeo, bxgeo,bygeo,bzgeo
      REAL*8 GS(105),HS(105)

c     The coefficients from Table 1 in Dougherty et al. (2018). Note,
c     the correction to the table has been incorporated.
      DATA GS/0.0D0, 21140.2D0,   0.0D0, 1581.1D0,   0.0D0,
     *        0.0D0,  2260.1D0,   0.0D0,    0.0D0,   0.0D0, 91.1D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,  12.6D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,  17.2D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,-59.6D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,   -10.5D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0, -12.9D0,  0.0D0, 
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,  15.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,    18.2D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.3D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,  0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0/
      DATA HS/0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0,    0.0D0,   0.0D0,   0.0D0,
     *        0.0D0,     0.0D0,   0.0D0/

c     Rotate the KSM position vector into KSMAG.
      call rotate_about_y(sunAngle, x,y,z, xgeo,ygeo,zgeo)
      
c     Since the field is axially symmetric, the magnetic field in KSMAG
c     is identical to the body fixed frame IAU Saturn coordinates.
      call IGRF_GEO_CART(gs,hs, xgeo,ygeo,zgeo, bxgeo,bygeo,bzgeo)

c     Now, rotate the fiedl to KSM
      call rotate_about_y(-sunAngle, bxgeo,bygeo,bzgeo,
     .     bx,by,bz)

      END
