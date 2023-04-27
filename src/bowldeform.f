      SUBROUTINE bowldeform (rH,sunAngle, x,y,z, xMap,yMap,zMap)
c
c     Peforms the Bowl deformation described by Arridge et al. (2008).
c     The implementation is based on K.Khurana's KMAG model.
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
c     real*8 xMap,yMap,zMap: The bowl deformed location.
c     
      IMPLICIT NONE
      
c     Inputs      
      REAL*8 rH,sunAngle, x,y,z

c     Outputs
      REAL*8 xMap,yMap,zMap
 
C     Local Variables
      REAL*8 rho,phi, zCS,rCS, sinThCS,cosThCS, rhoMap

c     convert to cylindrical coords
      call cart_to_cyl(x,y,z, rho, phi)
c     previously we used NAIF-SPICE's reccyl routine, but this was
c     replaced due to licensing issues
c     call reccyl(x,y,z, rho, phi)

c     the height of the current sheet surface
      zCS = (rho-rH*DTANH(rho/rH))*DTAN(sunAngle)

c     the radius to the current sheet from the orign
      rCS = DSQRT(rho*rho + zCS*zCS)

c     thetaCS is the angle of the location of the current sheet
      sinThCS = zCS/rCS
      cosThCS = rho/rCS

c     if the current sheet radius is very small, set these to zero
      if (rCS .lt. 1.0E-12) then
         sinThCS = 0.0d0
         cosThCS = 1.0d0
      endif

c     map the rho and z coordinate by rotating by the current sheet
c     angle, this is similar to Arridge's transformation, but maintains
c     the length of the deformed vector
      rhoMap = rho*cosThcs - z*sinThcs;
      zMap   = rho*sinThcs + z*cosThcs;

c     convert back to Cartesian coords
      call cyl_to_cart(rhoMap,phi,zMap, xMap,yMap)

c     previously we used NAIF-SPICE's cylrec routine, but this was
c     replaced due to licensing issues      
c     call cylrec(rhoMap,phi,zMap, xMap,yMap)
      
      return
      END

      
      SUBROUTINE bowldeform_field (rH,sunAngle, xMap,yMap,zMap,
     .     bx,by,bz, bxMap,byMap,bzMap)
c
c     Peforms the Bowl deformation described by Arridge et al. (2008).
c     The implementation is based on K.Khurana's KMAG model.
c
c     Inputs:
c     real*8 rH: The hinging distance in units of Saturn radii.
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 xMap,yMap,zMap: The bowl deformed location.
c     real*8 bx,by,bz: The undeformed field evaluated at the bowl
c       deformed location.
c
c     Outputs:
c     real*8 bxMap,byMap,bzMap: The bowl deformed field.
c   
      IMPLICIT NONE
      
c     Inputs      
      REAL*8 rH,sunAngle, xMap,yMap,zMap, bx,by,bz

c     Outputs
      REAL*8 bxMap,byMap,bzMap
 
C     Local Variables
      REAL*8 xpp,ypp,zpp, xpm,ypm,zpm
      REAL*8 dxpdx,dypdx,dzpdx, dxpdy,dypdy,dzpdy, dxpdz,dypdz,dzpdz
      REAL*8 Txx,Txy,Txz, Tyx,Tyy,Tyz, Tzx,Tzy,Tzz
      REAL*8 dr
      PARAMETER (dr=0.01d0)
      
      call bowldeform(rH,sunAngle, xMap+dr,yMap,zMap, xpp,ypp,zpp)
      call bowldeform(rH,sunAngle, xMap-dr,yMap,zMap, xpm,ypm,zpm)

      dxpdx = (xpp-xpm)/(2.0d0*dr)
      dypdx = (ypp-ypm)/(2.0d0*dr)
      dzpdx = (zpp-zpm)/(2.0d0*dr)
      
      call bowldeform(rH,sunAngle, xMap,yMap+dr,zMap, xpp,ypp,zpp)
      call bowldeform(rH,sunAngle, xMap,yMap-dr,zMap, xpm,ypm,zpm)

      dxpdy = (xpp-xpm)/(2.0d0*dr)
      dypdy = (ypp-ypm)/(2.0d0*dr)
      dzpdy = (zpp-zpm)/(2.0d0*dr)
      
      call bowldeform(rH,sunAngle, xMap,yMap,zMap+dr, xpp,ypp,zpp)
      call bowldeform(rH,sunAngle, xMap,yMap,zMap-dr, xpm,ypm,zpm)
      
      dxpdz = (xpp-xpm)/(2.0d0*dr)
      dypdz = (ypp-ypm)/(2.0d0*dr)
      dzpdz = (zpp-zpm)/(2.0d0*dr)
      
C     Now calculate the T matrix
c     Calculate the mapped location
      Txx = dypdy*dzpdz-dypdz*dzpdy
      Txy = dxpdz*dzpdy-dxpdy*dzpdz
      Txz = dxpdy*dypdz-dxpdz*dypdy
      Tyx = dypdz*dzpdx-dypdx*dzpdz
      Tyy = dxpdx*dzpdz-dxpdz*dzpdx
      Tyz = dxpdz*dypdx-dxpdx*dypdz
      Tzx = dypdx*dzpdy-dypdy*dzpdx
      Tzy = dxpdy*dzpdx-dxpdx*dzpdy
      Tzz = dxpdx*dypdy-dxpdy*dypdx
      
C     Now calculate the field at the mapped location 
      Bxmap = Txx*Bx+Txy*By+Txz*Bz	
      Bymap = Tyx*Bx+Tyy*By+Tyz*Bz	
      Bzmap = Tzx*Bx+Tzy*By+Tzz*Bz
      
      return
      END
      
