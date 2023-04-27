      subroutine cart_field_to_cyl(x,y,z, bx,by,bz,
     .     rho,phi, brho,bphi)
c
c     Converts a Cartesian position and field vector at that position to
c     a cylindrical position and field vector.
c
c     Inputs:
c     real*8 x,y,z: The position vector to be converted to cylindrical.
c     real*8 bx,by,bz: The field vector to be converted to cylindrical.
c
c     Outputs:
c     real*8 rho,phi: The cylindrical position vector.
c     real*8 brho,bphi: The cylindrical field vector.
c     
      implicit none
      real*8 x,y,z, bx,by,bz
      real*8 rho,phi, brho,bphi

      real*8 costheta,sintheta,cosphi,sinphi

c     convert the Cartesian position to cylindrical
      call cart_to_cyl(x,y,z, rho,phi)
c     previously we used NAIF-SPICE's reccyl routine, but this was
c     replaced due to licensing issues
c     call recsph(x,y,z, r,theta,phi)
      
      cosphi = DCOS(phi)
      sinphi = DSIN(phi)

c     conver the Cartesian field to spherical
      brho   =  cosphi*bx + sinphi*by
      bphi   = -sinphi*bx + cosphi*by
      
      end

      
      subroutine cart_to_cyl(x,y,z, rho,phi)
c
c     Converts a Cartesian position vector to a cylindrical position
c     vector.
c
c     Inputs:
c     real*8 x,y,z: The position vector to be converted to cylindrical.
c
c     Outputs:
c     real*8 rho,phi: The cylindrical position vector.
c 
      implicit none
      real*8 x,y,z
      real*8 rho,phi

      rho   = DSQRT(x**2+y**2)
      phi   = DATAN2(y, x)

c     If we are along the Z-axis, phi is undefined. Gfortran DATAN2's
c     implementation sets this to zero, however, other compilers might
c     not, so let's add this check in.
      if (x .eq. 0.0d0 .and. y .eq. 0.0d0) then
         phi = 0.0d0
      endif
      
      end


      subroutine cyl_to_cart(rho,phi,z, x,y)
c
c     Converts a Cartesian position vector to a cylindrical position
c     vector.
c
c     Inputs:
c     real*8 x,y,z: The position vector to be converted to cylindrical.
c
c     Outputs:
c     real*8 rho,phi: The cylindrical position vector.
c 
      implicit none
      real*8 x,y,z
      real*8 rho,phi

      x = rho * DCOS(phi)
      y = rho * DSIN(phi)
      
      end
      
