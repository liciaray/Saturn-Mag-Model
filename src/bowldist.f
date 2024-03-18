      SUBROUTINE bowldist (rH,sunAngle, x,y,z, xCS,yCS,zCS, distCS)
c
c     Computes the minimum distance form the supplied point (x,y,z) to
c     the Arridge et al. (2008) Bowl current sheet.
c
c     Inputs:
c     real*8 rH: The hinging distance in units of Saturn radii.
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 x,y,z: The location to be bowl deformed in the KSMAG
c       coordinate system in units of Saturn radii.
c
c     Outputs:
c     real*8 xCS,yCS,zCS: The point on the Bowl current sheet that is
c       the minimum distance from the supplied point
c     real*8 distCS: The minimum distance between the supplied point and
c       the Bowl current sheet.
c
      IMPLICIT NONE
      
c     Inputs      
      REAL*8 rH,sunAngle, x,y,z

c     Outputs
      REAL*8 xCS,yCS,zCS, distCS
 
C     Local Variables
      REAL*8 rho,phi, rhoCS, distSq
      REAL*8 distSquared,solveRhoCS
      
c     convert to cylindrical coords
      call cart_to_cyl(x,y,z, rho, phi)

c     solve the value of rhoCS that is the minimum distance to the Bowl
      rhoCS = solveRhoCS(rH,sunAngle, rho,z, rho)

c     the height of the current sheet at that value of rhoCS
      zCS = (rhoCS-rH*DTANH(rhoCS/rH))*DTAN(sunAngle)

c     compute the distance analytically
      distSq = distSquared(rH,sunAngle, rho,z, rhoCS)
      distCS = DSQRT(distSq)

c     now convert back to Cartesian
      call cyl_to_cart(rhoCS,phi,zCS, xCS,yCS)
      
      return
      END


      real*8 function distSquared(rH,sunAngle, rho,z, rhoCS)
      implicit none
c     The function returning squared distance between a point (rho, z)
c     and the Bowl current sheet (rhoCS, zCS).
      
c     Inputs
      real*8 rH,sunAngle
      real*8 rho,z
      real*8 rhoCS

c     Internal variables
      real*8 tanS,tanhR

      tanS = DTAN(sunAngle)
      
      tanhR = DTANH(rhoCS / rH)

      distSquared = (rhoCS-rho)**2+(rhoCS*tanS-rH*tanS*tanhR-z)**2

      return
      end


      real*8 function dDistSqDrho(rH,sunAngle, rho,z, rhoCS)
      implicit none
c     The derivative with respect to rhoCS of the function returning
c     squared distance between a point (rho, z) and the Bowl current
c     sheet (rhoCS, zCS).
      
c     Inputs
      real*8 rH,sunAngle
      real*8 rho,z
      real*8 rhoCS

c     Internal variables
      real*8 tanS,tanhR

      tanS = DTAN(sunAngle)
      
      tanhR = TANH(rhoCS / rH)

      dDistSqDrho = 2.0d0*(rhoCS - rho) +
     .     2.0d0*tanS*tanhR**2*(rhoCS*tanS - rH*tanS*tanhR - z)

      return
      end      


      real*8 function ddDistSqDrhoDrho(rH,sunAngle, rho,z, rhoCS)
      implicit none
c     The second derivative with respect to rhoCS of the function
c     returning squared distance between a point (rho, z) and the Bowl
c     current sheet (rhoCS, zCS).
      
c     Inputs
      real*8 rH,sunAngle
      real*8 rho,z
      real*8 rhoCS

c     Internal variables
      real*8 tanS,tanhR

      tanS = DTAN(sunAngle)
      
      tanhR = TANH(rhoCS / rH)

      ddDistSqDrhoDrho = 2.0d0 + (4.0d0*tanS/rH)*tanhR*(1 - tanhR**2)*
     .     (rhoCS*tanS - rH*tanS*tanhR - z)
     .     + 2.0d0*tanS**2*tanhR**4
      
      return
      end      


      real*8 function solveRhoCS(rH,sunAngle, rho,z, rho0)
      implicit none
c     solves RhoCS via the Newton-Raphson root finding method
      
c     Inputs
      real*8 rH,sunAngle
      real*8 rho,z
      real*8 rho0

c     Internal variables
      real*8 dDistSqDrho,ddDistSqDrhoDrho
      real*8 fx,dfxdy,xi
      integer numIter, i
      
      numIter = 5
      
c     the initial guess at the solution
      xi = rho0

c     loop over the number of iterations

      do i = 1, numIter

c     evaluate the function
         fx = dDistSqDrho(rH,sunAngle, rho,z, xi)

c     evaluate the derivative
         dfxdy = ddDistSqDrhoDrho(rH,sunAngle, rho,z, xi)

c     the next step of Newton-Raphson, x1 = x0 - f(x0)/f'(x0)
         xi = xi - (fx / dfxdy)
      enddo
      
      solveRhoCS = xi
      
      return
      end
