      SUBROUTINE saturn_ext (sunAngle, pDyn,
     .     x,y,z,
     .     bx,by,bz)
c      
c     Computes the external magnetic field for Saturn in the KSM
c     coordinate system in units of nT.
c
c     Inputs:
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 pDyn: The dynamic pressure of the solar wind measured in
c       nPa.
c     real*8 x,y,z: The location to evaluate the value of the external
c       magnetic field in the KSM coordinate system in units of Saturn
c       radii.
c      
c     Outputs:
c     real*8 bx,by,bz: The external magnetc field evaluated at the
c       supplied location in KSM coordinates in units of nT.
c      
      IMPLICIT NONE

c     Inputs
      REAL*8 sunAngle,pDyn, x,y,z
      
c     Outputs
      REAL*8 bx,by,bz
      
c     Internal variables
      REAL*8 scale, xScale,yScale,zScale, xKSMAG,yKSMAG,zKSMAG
      REAL*8 xMap,yMap,zMap
      REAL*8 bdx,bdy,bdz, btx,bty,btz, btMapx,btMapy,btMapz
      REAL*8 btKSMx,btKSMy,btKSMz

c     Constants
      REAL*8     a2
      PARAMETER (a2 = 1.0d0/5.7d0)
      REAL*8     pDyn0
      PARAMETER (pDyn0 = 0.017d0)
      INTEGER    Mend
      PARAMETER (Mend = 6)
      INTEGER    Nend
      PARAMETER (Nend = 10)
      REAL*8     rho0
      PARAMETER (rho0 = 50.0d0)
      REAL*8     D0
      PARAMETER (D0 = 2.0d0)
      REAL*8     rH
      PARAMETER (rH = 23.0d0)   

c     the coefficients for the equatorial current system
      REAL*8     sym(Nend)
      DATA       sym/
     .     481.7854461018803d0, -551.3012456231409d0,
     .     112.14327344552166d0, 295.2815270708938d0,
     .     -29.02961024189838d0, -154.9371611030641d0,
     .     69.54170605988301d0, 173.25497239890817d0,
     .     -284.68457458149567d0, 192.39276672456603d0/
      REAL*8     odd(Mend,Nend)
      DATA       odd/
     .     -10.772735477793702d0, -18.913225357221574d0,
     .     -48.2695591097975d0, -103.9309093075073d0,
     .     -138.40454388326856d0, 150.41637104937172d0,
     .     -6.127128886733822d0, 35.65996621137313d0,
     .     40.45335520559493d0, 33.66986718333211d0,
     .     6.210935743280192d0, -62.747101290784805d0,
     .     21.422231619553973d0, -52.45436354743525d0,
     .     -3.3697224089465796d0, 6.286741481268737d0,
     .     18.96994762578857d0, 41.59521930320715d0,
     .     28.071107342608876d0, 32.842695564806824d0,
     .     -34.4929061808727d0, 0.6634148029631938d0,
     .     -15.31828240036242d0, -30.77775696446307d0,
     .     -69.80398539141038d0, 6.873233466093325d0,
     .     8.751639918416458d0, -27.443158934023693d0,
     .     26.32137605336531d0, 28.91249705924795d0,
     .     43.3459723752866d0, -5.416099666571757d0,
     .     20.137734327120008d0, -11.755111692546343d0,
     .     -60.606189154611954d0, -26.98316137404707d0,
     .     -1.6018728351918827d0, 24.77895770885043d0,
     .     41.46746822396604d0, 82.8822996837279d0,
     .     66.47817916095046d0, 15.940555780177203d0,
     .     -60.38202763747571d0, -98.17995618671314d0,
     .     -99.06855450696239d0, -81.89719204488559d0,
     .     -35.08179572838114d0, -5.122772745558465d0,
     .     98.7867290543529d0, 101.17759258153403d0,
     .     62.19014479466227d0, 31.539420297061064d0,
     .     8.638743111784425d0, 0.5208609078240295d0,
     .     -44.86242371561234d0, -32.563760486805315d0,
     .     -11.70076269351273d0, -3.8597513943662967d0,
     .     -0.9757210479190807d0, 0.3068128786804058d0/
      REAL*8     even(Mend,Nend)
      DATA       even/
     .     -11.532529083847432d0, -29.74105183427352d0,
     .     -61.032521949546606d0, -81.70749507257412d0,
     .     -35.66624193597269d0, 94.16878399142774d0,
     .     55.29720817705375d0, 53.89552262220369d0,
     .     37.58600474038564d0, 6.165332872061285d0,
     .     -19.391240179605223d0, -26.433254887191367d0,
     .     -95.21572008753628d0, -23.124880633880693d0,
     .     20.79816813675172d0, 31.15986472297058d0,
     .     19.05885493809634d0, 10.470451417526323d0,
     .     49.15090595708274d0, -59.25101516590921d0,
     .     -53.14242148476542d0, -15.22998205391208d0,
     .     1.6195721223542883d0, -2.289021327580942d0,
     .     48.25339342611149d0, 76.85649980638112d0,
     .     -5.2016689241359d0, -34.933131053107715d0,
     .     -15.31456754837573d0, 1.0030269547635804d0,
     .     -40.62701458717349d0, 19.771740082533668d0,
     .     81.98025531751843d0, 41.18527983949692d0,
     .     2.0274312297191743d0, -4.090512083970708d0,
     .     -73.4158651368815d0, -88.4410883948333d0,
     .     -62.11954574273583d0, 3.6932040820339638d0,
     .     12.254802971786683d0, 4.050387973087059d0,
     .     103.74036722405172d0, 41.30998575895324d0,
     .     -5.256871567924476d0, -24.69042158551766d0,
     .     -8.118147357920629d0, -1.8329206669403135d0,
     .     -31.656859679351555d0, 13.31740167739348d0,
     .     20.914947178972398d0, 10.347975635856265d0,
     .     1.5339086071185215d0, 0.6549925939176265d0,
     .     -6.836837023655434d0, -10.39250559614622d0,
     .     -5.240615493693548d0, -0.6678768137058495d0,
     .     -0.4291085077697023d0, -0.03538317559905341d0/

c     Perform the self-similar scaling to the position vector.
c     This is the same scaling that is done to the magnetopause.
      scale = (pDyn/pDyn0)**a2
      xScale = scale*x
      yScale = scale*y
      zScale = scale*z      

c     convert the scaled KSM coordinate to KSMAG
      call rotate_about_y(sunAngle,
     .     xScale,yScale,zScale,
     .     xKSMAG,yKSMAG,zKSMAG)

c     evaluate the dipole shielding field
      call dipoleshield(sunAngle,
     .     xScale,yScale,zScale,
     .     bdx,bdy,bdz)
c     scale the field's magnitude
      bdx = bdx*scale*scale*scale
      bdy = bdy*scale*scale*scale
      bdz = bdz*scale*scale*scale

c     compute the bowl deformed position, input is KSMAG coordinates
      call bowldeform(rH,sunAngle,
     .     xKSMAG,yKSMAG,zKSMAG,
     .     xMap,yMap,zMap)

c     evaluate the equatorial currents
      call tailsheet_sum(Mend,Nend,rho0,D0, sym,odd,even,
     .     xMap,yMap,zMap,
     .     btx,bty,btz)

c     now deform the field, this is based on Tsy.'s general deformation
c     technique
      call bowldeform_field(rH,sunAngle,
     .     xMap,yMap,zMap, btx,bty,btz, btMapx,btMapy,btMapz)

c     now, rotate tail field back to KSM coordinates
      call rotate_about_y(-sunAngle,
     .     btMapx,btMapy,btMapz,
     .     btKSMx,btKSMy,btKSMz)

c     the total field is the sum of the dipole shielding field and the
c     equatorial field
      bx = bdx+btKSMx
      by = bdy+btKSMy
      bz = bdz+btKSMz      
      
      END


      logical function inside_magpause(sunAngle, pDyn,
     .     x,y,z)
c      
c     Determines if the supplied location is inside (returns true) or
c     outside or on (returns false) the modeled magnetopause surface.
c     Note, the magnetopause surface is warped using the bowl
c     deformation and is rescaled using self-similar pressure rescaling
c     in the same fashion as the other magnetospheric current systems.
c
c     Inputs:
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 pDyn: The dynamic pressure of the solar wind measured in
c       nPa.
c     real*8 x,y,z: The location to evaluate the value of the external
c       magnetic field in the KSM coordinate system in units of Saturn
c       radii.
c      
c     Outputs:
c     logical inside_magpause: Returns true if the supplied location is
c     inside the modeled magnetopause and false if it is outside or on
c     the surface.
c      
      IMPLICIT NONE

c     Inputs
      REAL*8 sunAngle,pDyn, x,y,z
      
c     Internal variables
      REAL*8 scale, xScale,yScale,zScale, xKSMAG,yKSMAG,zKSMAG
      REAL*8 xMap,yMap,zMap, rMap
      REAL*8 r0,theta,rMagpause

c     Constants
      REAL*8     a1
      PARAMETER (a1 = 10.5d0)
      REAL*8     a2
      PARAMETER (a2 = 1.0d0/5.7d0)
      REAL*8     a3
      PARAMETER (a3 = 0.67d0)
      REAL*8     a4
      PARAMETER (a4 = 0.17d0)
      REAL*8     pDyn0
      PARAMETER (pDyn0 = 0.017d0)
      REAL*8     rH
      PARAMETER (rH = 23.0d0)   
      REAL*8     plasmaBeta
      PARAMETER (plasmaBeta = 3.0d0)  
      
c     Perform the self-similar scaling to the position vector.
c     This is the same scaling that is done to the magnetopause.
      scale = (pDyn/pDyn0)**a2
      xScale = scale*x
      yScale = scale*y
      zScale = scale*z      

c     convert the scaled KSM coordinate to KSMAG
      call rotate_about_y(sunAngle,
     .     xScale,yScale,zScale,
     .     xKSMAG,yKSMAG,zKSMAG)

c     compute the bowl deformed position, input is KSMAG coordinates
      call bowldeform(rH,sunAngle,
     .     xKSMAG,yKSMAG,zKSMAG,
     .     xMap,yMap,zMap)

c     compute the radius of the bowl deformed position
      rMap = DSQRT(xMap**2+yMap**2+zMap**2)

c     if the radius is zero, return true
      if (rMap .eq. 0.0d0) then
         inside_magpause = .true.
      else

c     Now evaluate the Pilkington et al. (2015) model. Note, we evaluate
c     using the baseline dynamic pressure (pDyn0) and not the
c     instantaneous dynamic pressure (pDyn). This is because the dynamic
c     pressure rescaling has already been performed above.
         theta=DATAN2(DSQRT(yMap**2+zMap**2),xMap)

c     the values of Pilkington using the actual value of Pdyn
c     r0 = a1*(pDyn/(1.0d0+plasmaBeta))**(-a2)
c     rMagpause = r0*(2.0d0/(1.0d0+COS(theta)))**(a3+a4*pdyn)

c     eq. (12) 
         r0 = a1*(pDyn0/(1.0d0+plasmaBeta))**(-a2)
c     eq. (1)
         rMagpause = r0*(2.0d0/(1.0d0+COS(theta)))**(a3+a4*pDyn0)

         inside_magpause = rMap .lt. rMagpause
      endif
      
      end

      
