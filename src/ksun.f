      Subroutine KSun(time,stheta,sphi,Ztheta,Zphi)
c     INPUT:  J2000 time of the data point 
C     OUTPUTS: stheta, sphi, latitude and longitude  (in radians) of the Sun in system III (RH).  
C     OUTPUTS: Ztheta, Zphi, latitude and longitude  (in radians) of the Zkso in system III (RH).  
c     The equations are written in J2000 convention followed by PDS. 
c     
      
c     Last updated August 26, 2004.
c     
c     theta=a1*cos(omegay*t)+a2*sin(omegay*t)+
c     a3*cos(2.*omegay*t)+a4*sin(2.*omegay*t)+
c     a5*cos(3.*omegay*t)+a6*sin(3.*omegay*t)+
c     a7*cos(4.*omegay*t)+a8*sin(4.*omegay*t)+
c     a9*t**2+a10*t+a11
c     fphi=b1*cos(omegay*t)+b2*sin(omegay*t)+
c     b3*cos(2.*omegay*t)+b4*sin(2.*omegay*t)+
c     b5*cos(3.*omegay*t)+b6*sin(3.*omegay*t)+
c     b7*cos(4.*omegay*t)+b8*sin(4.*omegay*t)+
c     b9*t**2+b10*t+b11	(I assume I have despun Saturn first.)
      
c     Then we rotate into the System III coordinates
      
c     time is a double precision variable and is assumed to be in J2000. theta and phi are single precision 
c     variables.
c     
      
c     
c     omega is Saturn's rotation rate
c     omegaz is saturn's rotation rate for use of Zkso lat and lon
c     omegay is Saturn's yearly orbital rate
      
c     Initialize variables
      Implicit Real*8(A-H,O-Z)
      Real*8 aa(11),bb(11), cc(11),x(11)
      Real*8 time,jtime1,fphi,yrsat,omega,D360,omegay,t,year,omegaz
c     Real*4 stheta,sphi,ztheta,Zphi
      Parameter (PI=3.141592654d0)
      PARAMETER (twopi=2.*PI,radian=PI/180.,degree=180./PI)
      Parameter (yrsat=10759.22d0*86400.0d0)
      Parameter (omega=360.D0/(10.65622222D0*3600.D0))
      PARAMETER (omegay=2.*PI/yrsat)
      Parameter (omegaz=360.D0/(10.65622222D0*3600.D0)+3.880688963D-7)
      Parameter (year=86400D0*.36525D3,D360=360.D0)
      Parameter (jtime1=-.825767955816D+09)
      Parameter (Zthetd=63.271d0)
c     Parameter (Zthetd=-26.729)
      Data aa  /-.26237934D+02,.30525104D+01,-.13686733D+01,
     +     .10182726D+00,-.30897805D+00,.89277949D-01,-.39638771D-01,
     +     .97499653D-02,.40974159D-04,.14075684D-02,.13505775D+01/
c     
      Data bb /-.32108124D+00,.58569866D+01,.71266272D+00,.33244539D+01,
     +     .47529267D-01,.32770362D+00,.46935622D-01,.10720536D+00,
     +     -.50594764D-03,.29439612D-01,.26423581D+03/
c
      Data cc/-.13200649D-02,-.18358429D-02,.64658927D-03,.12800487D-02,
     +     .17618936D-04,-.58790898D-03,.49804081D-04,.42372137D-03,
     +     .14744891D-04,.25369744D-01,.77943328D+02/

c     First calculate the latitude and longitude in non-rotating Saturn coordinates.
      
c     Calculate the best fit theta and fphi
      t=time-jtime1
      x(1)=cos(omegay*t)
      x(2)=sin(omegay*t)
      x(3)=cos(2.*omegay*t)
      x(4)=sin(2.*omegay*t)
      x(5)=cos(3.*omegay*t)
      x(6)=sin(3.*omegay*t)
      x(7)=cos(4.*omegay*t)
      x(8)=sin(4.*omegay*t)
      x(9)=(t/year)**2
      x(10)=t/year
      x(11)=1.0
      stheta=0.
c     fphi is phi in Saturn fixed (non-rotating) coordinate
      fphi=0.0
      Zfphi=0.0
      Do j=1,11
         fphi=fphi+bb(j)*x(j)
         Zfphi=Zfphi+cc(j)*x(j)
         stheta=stheta+aa(j)*x(j)
      enddo
C     Now rotate the longitude to Saturn System III
c     First Add the rotation of Saturn around the Sun.
c     fphi is the phi of the Sun as unspinning Saturn goes around the sun
      fphi=DMod(fphi+t/yrsat*360.d0, D360)
c     Next add the rotation of Saturn around its axis.
      sphi=DMod(fphi-t*omega, D360)
      if (sphi .lt. 0.0) sphi = sphi+360.0
      sphi=sphi*radian
      stheta=stheta*radian
c     
C     Next rotate the longitude of Zfphi to Saturn System III
c     First Add the rotation of Saturn around the Sun.
c     Zfphi is the phi of Zkso as unspinning Saturn goes around the sun
      Zfphi=DMod(Zfphi+t/yrsat*360.d0+180.D0, D360)
c     Next add the rotation of Saturn around its axis.
      Zphi=DMod(Zfphi-t*omegaz, D360)
      if (Zphi .lt. 0.0) Zphi = Zphi+360.0
      Zphi=Zphi*radian
      Ztheta=Zthetd*radian
      Return
      end
      
      
