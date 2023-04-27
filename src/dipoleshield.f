C     May 2019, G.K.Stephens, this subroutine was adabpted to the Saturn
C     magnetic field model.
C      
C     THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
C     REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
C     to the z=0 plane (see NB#4, p.74-74)
C     
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c     harmonics (A(1)-A(36).
c     The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
C     entering the arguments of exponents, sines, and cosines in each of the
C     18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
C     (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE  dipoleshield(sunAngle, x,y,z, bx,by,bz)
c     
c     Computes the dipole shielding magnetic field for Saturn in the KSM
c     coordinate system in units of nT.
c
c     Inputs:
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 x,y,z: The location to evaluate the value of the dipole
c       sheilding magnetic field in the KSM coordinate system in units
c       of Saturn radii.
c     
c     Outputs:
c     real*8 bx,by,bz: The dipole shielding magnetc field evaluated at
c       the supplied location in KSM coordinates in untis of nT.
c      
      IMPLICIT NONE
      
c     Inputs
      REAL*8 sunAngle
      REAL*8 x,y,z
      
c     Outputs
      REAL*8 bx,by,bz
      
c     Internal variables
      REAL*8 pis(1:5)
      REAL*8 pisI(1:5)
      DATA   pis /51.361793505102675d0,27.451107992675237d0,
     .     36.93371761990351d0, 38.56965680962369d0,
     .     31.628812615079276d0/
      REAL*8 rks(1:5)
      REAL*8 rksI(1:5)
      DATA   rks /149.22817690462293d0, 179.84124219936155d0,
     .     92.34168588076675d0, 235.88621885954436d0,
     .     60.82606721603004d0/
      REAL*8 aPerp(1:5,1:5)
      REAL*8 aik(1:5,1:5)
      DATA   aik /
     .     -53.70374727845774d0, 312.2566262198816d0,
     .     69.69840316544105d0,  70.27085743445696d0,
     .     93.54067986428709d0,  -164.23204211289703d0,
     .     387.66685685911216d0, 116.86278373609093d0,
     .     102.35198517604658d0, 176.508253300266d0,
     .     307.38282640310354d0, -395.6673868510843d0,
     .     -372.8749339742717d0, -300.55001861609344d0,
     .     -537.132276940567d0,  -235.87848137380206d0,
     .     394.66152467047505d0, 131.0801800266836d0,
     .     107.20536964395887d0, 210.25214021584543d0,
     .     18.939385790228698d0, -31.778921060908942d0,
     .     127.03672231573d0,    112.54092955558008d0,
     .     97.260143858839d0/
      REAL*8 bik(1:5,1:5)
      DATA   bik /
     .     27.105441699870426d0, 150.23119173362284d0,
     .     93.9432666093453d0,   110.96256516927315d0,
     .     33.915890370230045d0, -69.65251712528698d0,
     .     240.29760859489033d0, 141.94178459575778d0,
     .     143.79174850250274d0, 122.64369638420158d0,
     .     254.7539318088966d0,  -529.7266885004938d0,
     .     -357.36017097142394d0,-279.2478065603482d0,
     .     -581.1350663421326d0, -141.72436288789322d0,
     .     273.16590913384425d0, 154.4368670396916d0,
     .     145.1265456670062d0,  166.34183727949312d0,
     .     -212.44018034313922d0,129.7694371946136d0,
     .     185.10747602555784d0, 146.3140099168013d0,
     .     212.32970976208162d0/
     
      REAL*8 qis(1:5)
      REAL*8 qisI(1:5)
      DATA   qis /186.33341088948478d0, 93.40115095973002d0,
     .     132.81448110260635d0, 135.31562338507575d0,
     .     227.74480939593587d0/
      REAL*8 sks(1:5)
      REAL*8 sksI(1:5)
      DATA   sks /34.01116096397626d0, 106.83007664071285d0,
     .     65.4777551388745d0, 98.94999816905575d0,
     .     31.807773086282975d0/
      REAL*8 aParr(1:5,1:5)
      REAL*8 cik(1:5,1:5)
      DATA   cik /
     .     89.92041440636967d0,  -575.4243264084216d0,
     .     97.58516581146978d0,  104.46133244899102d0,
     .     37.77903186756157d0,  165.73640203270043d0,
     .     151.34367880201899d0, 51.007805158020346d0,
     .     50.903186128271045d0, 290.5710359282093d0,
     .     -29.121352129528532d0,-57.310464271955425d0,
     .     153.56074310033d0,    144.89370246237377d0,
     .     -116.87466278026113d0,-186.34612167128944d0,
     .     -45.057424933620496d0,-195.67160155097372d0,
     .     -201.1550737285579d0, -117.88846556644512d0,
     .     -155.45344871256384d0,305.3987023318109d0,
     .     208.59246082924074d0, 188.64275516575435d0,
     .     -327.44286062114406d0/
      REAL*8 dik(1:5,1:5)
      DATA   dik /
     .     203.18810868973378d0, -1030.3886395706795d0,
     .     217.9443888757378d0,  230.67325245414395d0,
     .     106.1785151171207d0,  285.15162330379826d0,
     .     258.60346450685756d0, 74.71416954344022d0,
     .     74.48663448976004d0,  514.5125080734724d0,
     .     -34.66126675459964d0, -46.7600919904653d0,
     .     318.8601920182118d0,  301.76788999920245d0,
     .     -203.06449031824013d0,-367.1076049407129d0,
     .     -97.57527641340857d0, -379.70399937813636d0,
     .     -390.09193960169796d0,-242.98053284783964d0,
     .     -323.2289471644326d0, 512.7522877848241d0,
     .     346.1253883028403d0,  309.57017752557294d0,
     .     -640.3028149751481d0/

      integer i,k
      REAL*8 T1,T2, tilt,cosTilt,sinTilt,sin2Tilt, anglePerp,angleParr
      REAL*8 xr1,yr1,zr1, xr2,yr2,zr2
      REAL*8 bx1,by1,bz1, bxr1,byr1,bzr1
      REAL*8 bx2,by2,bz2, bxr2,byr2,bzr2

c     the dipole tilt angle is the same as the sun latitude angle
      tilt = sunAngle
      
c     invert the coefficients, the coefficients are relatively small,
c     so the inverse is fit instead, so before we evaluate them, we must
c     invert them back
      call invert_a(pis, 1,5, pisI)
      call invert_a(rks, 1,5, rksI)
      call invert_a(qis, 1,5, qisI)
      call invert_a(sks, 1,5, sksI)
      
      T1  = 0.6479062130963051d0
      T2  = 0.7817092718726241d0
      
      cosTilt = DCOS(tilt)
      sinTilt = DSIN(tilt)
      sin2Tilt = 2.0d0*sinTilt*cosTilt

c     Compute the coefficients that are a function of the dipole tilt
c     angle
      do i=1,5
         do k=1,5
            aPerp(i,k) = aik(i,k) + bik(i,k)*cosTilt
            aParr(i,k) = cik(i,k)*sinTilt + dik(i,k)*sin2Tilt
         enddo
      enddo
      
      anglePerp = tilt*T1
      angleParr = tilt*T2

c     MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
      bx1 = 0.0d0
      by1 = 0.0d0
      bz1 = 0.0d0
      call rotate_about_y(anglePerp, x,y,z, xr1,yr1,zr1)
      call cartharmonic(.true., .false., 5,5,
     .     pisI,rksI,aPerp,
     .     xr1,yr1,zr1,
     .     bx1,by1,bz1)
      call rotate_about_y(-anglePerp, bx1,by1,bz1, bxr1,byr1,bzr1)
      
      
c     MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
      bx2 = 0.0d0
      by2 = 0.0d0
      bz2 = 0.0d0
      call rotate_about_y(angleParr, x,y,z, xr2,yr2,zr2)
      call cartharmonic(.true., .true., 5,5,
     .     qisI,sksI,aParr,
     .     xr2,yr2,zr2,
     .     bx2,by2,bz2)
      call rotate_about_y(-angleParr, bx2,by2,bz2, bxr2,byr2,bzr2)

c     add the perpendicular and parallel fields together
      bx = bxr1+bxr2
      by = byr1+byr2
      bz = bzr1+bzr2      
       
      RETURN
      END
      
      
      SUBROUTINE invert_a(array, firstindex,lastindex, arrayI)
c     
c     Inverts the elements of an array of doubles.
c
c     Inputs:
c     real*8 array: The array whose elements will be inverted and
c       returned.
c     integer firstindex: The first index of the array.
c     integer lastindex: The last index of the array.
c     
c     Outputs:
c     real*8 arrayI: An array whose elements will be the inverse of the
c       elements of the supplied array.
c      
      
c     Inputs      
      integer i,firstindex,lastindex
      real*8 array(*)

c     Outputs
      real*8 arrayI(*)
      
      do i=firstindex,lastindex
         arrayI(i) = 1.0d0/array(i)
      enddo
      
      return
      END
