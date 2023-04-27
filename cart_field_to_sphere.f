      subroutine cart_field_to_sphere(x,y,z, bx,by,bz,
     .     r,theta,phi, br,btheta,bphi)
c
c     Converts a Cartesian position and field vector at that position to
c     a spherical position and field vector.
c
c     Inputs:
c     real*8 x,y,z: The position vector to be converted to spherical.
c     real*8 bx,by,bz: The field vector to be converted to spherical.
c
c     Outputs:
c     real*8 r,theta,phi: The spherical position vector.
c     real*8 br,btheta,bphi: The spherical field vector.
c     
      implicit none
      real*8 x,y,z, bx,by,bz
      real*8 r,theta,phi, br,btheta,bphi

      real*8 costheta,sintheta,cosphi,sinphi

c     convert the Cartesian position to spherical
      call cart_to_sphere(x,y,z, r,theta,phi)
c     previously we used NAIF-SPICE's recsph routine, but this was
c     replaced due to licensing issues
c      call recsph(x,y,z, r,theta,phi)

      costheta = DCOS(theta)
      sintheta = DSIN(theta)
      
      cosphi = DCOS(phi)
      sinphi = DSIN(phi)

c     conver the Cartesian field to spherical
      br     = sintheta*cosphi*bx + sintheta*sinphi*by + costheta*bz
      btheta = costheta*cosphi*bx + costheta*sinphi*by - sintheta*bz
      bphi   =         -sinphi*bx +          cosphi*by
      
      end

      
      subroutine cart_to_sphere(x,y,z, r,theta,phi)
c
c     Converts a Cartesian position vector to a spherical position
c     vector.
c
c     Inputs:
c     real*8 x,y,z: The position vector to be converted to spherical.
c
c     Outputs:
c     real*8 r,theta,phi: The spherical position vector.
c 
      implicit none
      real*8 x,y,z
      real*8 r,theta,phi
      
      r     = DSQRT(x**2+y**2+z**2)
      theta = DACOS(z / r)
      phi   = DATAN2(y, x)

c     If we are along the Z-axis, phi is undefined. Gfortran DATAN2's
c     implementation sets this to zero, however, other compilers might
c     not, so let's add this check in.
      if (x .eq. 0.0d0 .and. y .eq. 0.0d0) then
         phi = 0.0d0
      endif
      
c     If the radius is zero, the conversion is undefined, let's just
c     return all zero's for now.
      if (r .eq. 0.0d0) then
         theta = 0.0d0
         phi   = 0.0d0
      endif
      
      end
