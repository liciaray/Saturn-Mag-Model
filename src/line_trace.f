c     Jan 2025 Licia Ray
      SUBROUTINE line_trace (sunLat, pDyn,
     .     inside,
     .     x,y,z,s,
     .     bx,by,bz,
     .     hx,hy,hz,
     .     bTotx,bToty,bTotz,bMag,
     .     r,theta,phi,
     .     direction, geometry,
     .     geo_file)
c
c     This is workhorse to trace Saturn's magnetic field to the planet's
c     surface. The planetary coordinates of the footprint and magnetude
c     of the magnetic field are provided. It uses the JPL Saturn-Mag_Field
c     package. Units are in Saturn radii (Rs) and nanoTesla
c     (nT).
c
c     The coordinates and magnitude of the field can be written into
c     an output file for visual verification
c
c     The driver is find_footprint.for
c     
      implicit none

      integer io
      integer direction
      
      real*8 pDyn

      real*8 sunLat

      logical inside,inside_magpause

      logical geometry
      
      real*8 x,y,z,s,r_oblate
      
      real*8 bx,by,bz, hx,hy,hz, bTotx,bToty,bTotz
      real*8 r,theta,phi
      real*8 br,btheta,bphi
      real*8 bMag, dx,dy,dz, ds
      character(len=60) geo_file

      ds = 1000.0

           
      if (geometry) then

         io = 15

         open(unit=io, file=geo_file, status="new", action="write")
         write(io,'(A,A,A,A,A,A)')
     .        "x               ",
     .        "y               ",
     .        "z               ",
     .        "r_evaluation    ",
     .        "line_length     ",
     .        "Bmag            "
         write(io,'(6F12.5)') x,y,z,r,s,Bmag
      endif

c      Write(*,*) "Coordinates in line trace:   ", x, y, z
      

      r_oblate = DSQRT(1./(1.+(1./(1-1./11.1)**2 - 1)*
     .     (DCOS(theta)**2)))


c     Trace the field to an altitude of 10000 above the surface or
c     until the trace exceeds 90 RS in length as it's likely that 
c     the code has gone off piste by then, often due to proximity
c     to the magnetopause boundary
      do while ((r > (0.166 + r_oblate)).and.(s < 90.0))

c     determine the unit vectors for the field to use in
c     tracing the field line towards the planet
         dx = bTotx/bMag
         dy = bToty/bMag
         dz = bTotz/bMag
      
         x = x + direction/(ds/(r))*dx
         y = y + direction/(ds/(r))*dy
         z = z + direction/(ds/(r))*dz

         s = s + (r/ds)
      
c     Evaluates the external magnetic field in KSM coordinates in units
c     of nT. As stated above, this subroutine will still return vectors
c     even if the position vector is outside the modeled magnetopause.
         call saturn_ext(sunLat, pDyn,
     .        x,y,z,
     .        bx,by,bz)

c     Evaluates the internal magnetic field (Dougherty et al. 2018) in
c     KSM coordinates in units of nT.
         call saturn_int(sunLat, x,y,z, hx,hy,hz)
      
c     The total magnetic field is the sum of the external and internal
c     fields.
         bTotx = bx + hx
         bToty = by + hy
         bTotz = bz + hz
         bMag  = DSQRT(bTotx**2+bToty**2+bTotz**2)

      
c     Converts the Cartesian coordinates to spherical coordinates.
         call cart_to_sphere(x,y,z,r,theta,phi)

         if (geometry) then
c     Write coordinates and field magnitude to file
            write(io,'(6F12.5)') x,y,z,r,s,Bmag
         endif

         r_oblate = DSQRT(1./(1.+(1./(1-1./11.1)**2 - 1)*
     .        (DCOS(theta)**2)))
         
      end do

      if (geometry) then
         close(io)
      endif
      

      end
      
