c     June 2019, G.K.Stephens: Initial prototype.
      SUBROUTINE line_trace (sunLat, pDyn,
     .     inside,
     .     x,y,z,
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
      
      real*8 x,y,z
      
      real*8 bx,by,bz, hx,hy,hz, bTotx,bToty,bTotz
      real*8 r,theta,phi
      real*8 br,btheta,bphi
      real*8 bMag, dx,dy,dz, ds
      character(len=37) outfilename
      character(len=27) geo_file

      ds = 600.0

      if (geometry) then

         io = 15

         outfilename = "../output/"//geo_file

         open(unit=io, file=outfilename, status="new", action="write")
         write(io,'(A,A,A,A,A)')
     .        "x               ",
     .        "y               ",
     .        "z               ",
     .        "r_evaluation    ",
     .        "Bmag            "
         write(io,'(5F12.5)') x,y,z,r,Bmag
      endif

      
      do while (r > (1.024 - 1./11.1*DABS(DCOS(theta))**2))

c     determine the unit vectors for the field to use in
c     tracing the field line towards the planet
         dx = bTotx/bMag
         dy = bToty/bMag
         dz = bTotz/bMag
      
         x = x + direction/ds*dx
         y = y + direction/ds*dy
         z = z + direction/ds*dz
      
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
            write(io,'(5F12.5)') x,y,z,r,Bmag
         endif
                 
      end do

      if (geometry) then
         close(io)
      endif
      

      end
      
