c     June 2019, G.K.Stephens: Initial prototype.
      program find_footprint
c
c     This is an example program that demonstrates how to evaluate both
c     the external and internal magnetic field (Dougherty et al. 2018)
c     for Saturn. The coordinate system is KSM (Kronocentric Solar
c     Magnetospheric) and units are in Saturn radii (Rs) and nanoTesla
c     (nT).
c
c     The program can be compiled using the following using gfortran:
c     
c     gfortran example_saturn.for -o example_saturn
c      
c     The expected output of this program is the following:
c     
c     Date:    2013-04-21 22:53:28.39823
c     TDT:        419856875.58223003     
c     sunLat:    0.31991904766201745     
c     (x,y,z):
c       -12.381273999999999    3.8901874300000001    2.1981273432987001
c     inside:      T
c     (bx,by,bz):
c        6.2702298378118524   -1.2576152407132852    7.0949159196904112
c     (hx,hy,hz):
c       0.74732969552778661   -1.1478438385940311   -9.4198143792161204
c     (r,theta,phi):
c        13.162874607051583    1.4030157792196860    2.8371604325738131
c     (br,btheta,bphi):
c       -7.7000229788157384    1.0538352925030736   0.19132520816042442
c     
      implicit none

      character*(512) pdynfile
      character*(512) eventsfile

      integer noEvents
      parameter (noEvents = 4)
      integer final_count

      integer i
      
      integer event_year(noEvents)
      integer event_month(noEvents)
      integer event_dom(noEvents)
      integer event_hour(noEvents)
      integer event_min(noEvents)
      real*8 event_x(noEvents)
      real*8 event_y(noEvents)
      real*8 event_z(noEvents)
      
      integer year,month,dom,hour,min
      integer direction
      real*8 sec

      real*8 pDyn

      real*8 utc_to_tdt,tdt
      
      real*8 sunLat,sunLong, zPhi,zTheta

      logical inside,inside_magpause
      logical geometry
      
      real*8 x,y,z
      
      real*8 bx,by,bz, hx,hy,hz, bTotx,bToty,bTotz
      real*8 r,theta,phi, br,btheta,bphi, bMag
      real*8 xCS,yCS,zCS, distCS

      character(len=60) geo_file

c     Set whether you want the field line coordinates saved for each field line
      geometry = .FALSE.
      
c     The file containing the dynamic pressure of the solar wind at
c     Saturn. This file is derived from the Tao et al. (2005) model and
c     was obtained from the AMDA archive.
      pdynfile = "Dp.txt"
      pdynfile = trim(pdynfile)

      eventsfile = "../input/test2.csv"
      eventsfile = trim(eventsfile)

c     Reads and stores the results of the Dynamic pressure file.
c     This call must be made once prior to evaluating the dynamic
c     pressure. Once called, this does not need called again.
      call store_pdyn(pdynfile)

      call read_events(eventsfile, event_year,
     .     event_month, event_dom, event_hour, event_min,
     .     event_x,event_y,event_z,final_count,noEvents)

      write(*,*) "Year =",event_year
      write(*,*) "Month =",event_month
      write(*,*) "Day =",event_dom
      write(*,*) "Hour =",event_hour
      write(*,*) "Min =",event_min
      write(*,*) "x =",event_x
      write(*,*) "y =",event_y
      write(*,*) "x =",event_z
      write(*,*) "Number events = ", noEvents
      
      do i = 1, noEvents
      
c     The UTC date and time to evaluate the model. The format follows
c     the standard year-month-dayofmonth format, where these fields
c     start from 1, e.g. month=1 corresonds to January, and the time
c     fields start at 0. Year, month, dom, hour, and min are integer,
c     while sec is a double.
      year = event_year(i)
      month = event_month(i)
      dom = event_dom(i)
      hour = event_hour(i)
      min = event_min(i)
      sec = 0d0

c     Create a filename for the field line coordinates if
c     they are to be saved
      if (geometry) then
         call create_filename(year,month,
     .        dom,hour,min,geo_file)
      endif
      
c     Evaluates the dynamic pressure at Saturn for the given time.
c     Pressures are linearly interpolated.
      call eval_pdyn(year,month,dom, hour,min,sec,
     .     pdyn)

c     Evalute the Terrestrial Dynamical Time (TDT or Terrestrial Time or
c     TT) for the given UTC date and time.
      tdt = utc_to_tdt(year,month,dom, hour,min,sec)

c     Evaluates the location of the Sun at the supplied TDT time in
c     radians.
      call ksun(tdt, sunLat,sunLong, zPhi,zTheta)

c     The position vector in KSM coordinates in units of Saturn radii.
      x = event_x(i)
      y = event_y(i)
      z = event_z(i)

c     Evalutes if the position vector is inside or outside the modeled
c     magnetopause (true if inside, false if on or outside). Note, if
c     the position is outside the magnetopause, the subroutines
c     saturn_ext and saturn_int will still return vectors, but these
c     vectors will not have any physical signficance. We leave it up to
c     the user to decide how to handle vectors outside the
c     magnetosphere.     
      inside = inside_magpause(sunLat, pDyn,
     .     x,y,z)

      if (.not. inside) then
         write(*,*) "Outside"
         stop
      endif
      
c     Evaluates the external magnetic field in KSM coordinates in units
c     of nT. As stated above, this subroutine will still return vectors
c     even if the position vector is outside the modeled magnetopause.
      call saturn_ext(sunLat, pDyn,
     .     x,y,z,
     .     bx,by,bz)

c     Evaluates the internal magnetic field (Dougherty et al. 2018) in
c     KSM coordinates in units of nT.
      call saturn_int(sunLat, x,y,z, hx,hy,hz)
      
c     The total magnetic field is the sum of the external and internal
c     fields.
      bTotx = bx + hx
      bToty = by + hy
      bTotz = bz + hz
      bMag  = DSQRT(bTotx**2+bToty**2+bTotz**2)

c     Converts the Cartesian total field value to spherical coordinates.
      call cart_field_to_sphere(x,y,z,
     .     bTotx,bToty,bTotz,
     .     r,theta,phi, br,btheta,bphi)

c     Determine if field is above or below current sheet using the
c     radial component of the field
      if (br < 0) then
         direction = 1
         write(*,*) "Will go to southern hemisphere"
      else
         direction = -1
         write(*,*) "Will go to northern hemisphere"
      end if

c     Evaluates the distance and the closest location from the supplied
c     location to the Bowl shaped current sheet.
c      call distance_to_sheet(sunLat, pDyn, x,y,z, xCS,yCS,zCS, distCS)
      
c     Prints out the results.
      write(*,'(A,I4,A,I0.2,A,I0.2,A,I0.2,A,I0.2,A,F8.5)')
     .     "Date:            ",year,"-",month,"-",dom,
     .     " ",hour,":",min,":",sec
      write(*,*) "TDT:                ",TDT
      write(*,*) "Dp:                 ",pdyn
      write(*,*) "sunLat:             ",sunLat
      write(*,*) "(x,y,z):            ",x,y,z
      write(*,*) "inside:             ",inside
      write(*,*) "(bx,by,bz):         ",bx,by,bz
      write(*,*) "(hx,hy,hz):         ",hx,hy,hz
      write(*,*) "(r,theta,phi):      ",r,theta*180/3.14,phi*180/3.14

c     Evaluates the total magnetic field at the provided location
      call line_trace(sunLat, pDyn,
     .     inside,
     .     x,y,z,
     .     bx,by,bz,
     .     hx,hy,hz,
     .     bTotx,bToty,bTotz,bMag,
     .     r,theta,phi,
     .     direction,geometry,
     .     geo_file)

c     Prints out the results.
      write(*,*) "Field footprint details"
      write(*,*) "(x,y,z):            ",x,y,z
      write(*,*) "(bx,by,bz):         ",bx,by,bz
      write(*,*) "(hx,hy,hz):         ",hx,hy,hz
      write(*,*) "(r,theta,phi):      ",r,theta*180./3.14,phi*180/3.14
      write(*,*) "(bMag):             ",bMag

      end do
      
      end program
      
