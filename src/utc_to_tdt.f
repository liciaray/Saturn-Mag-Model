      

      real*8 function utc_to_tdt(year,month,dom, hour,min,sec)
c
c     Converts a UTC date and time into the elapsed seconds since
c     Greenwich noon on January 1st 2000 (J2000 Epoch) in the TDT
c     (Terrestrial Dynamical Time) time system.
c
c     Inputs:
c     integer year: The standard calendar year, which starts from 1.
c     integer month: The standard calendar month, starting from 1, which
c       corresponds to January.
c     integer dom: The day of month, starting from 1, which accounts for
c       leap years.
c     integer hour: The hour of the day, ranging from 0 - 23.
c     integer min: The minute of the hour, ranging from 0 - 59.
c     real*8 sec: The second of the minute, ranging from 0.0 - <60.0 for
c       minutes not containing a leap second
c      
c     Outputs:
c     real*8 utc_to_tdt: Elapsed seconds since Greenwich noon on January
c       1st 2000 (J2000 Epoch) in the TDT (Terrestrial Dynamical Time)
c       time system.
c           
      implicit none
c     Inputs
      integer year,month,dom
      integer hour,min
      real*8 sec

c     Internal variables
      real*8 tai,utc_to_tai

c     compute the TAI time
      tai = utc_to_tai(year,month,dom, hour,min,sec)

c     now add 32.184 seconds
      utc_to_tdt = tai+32.184d0

      return
      end
      

      real*8 function utc_to_tai(year,month,dom, hour,min,sec)
c
c     Converts a UTC date and time into the elapsed seconds since
c     Greenwich noon on January 1st 2000 (J2000 Epoch) in the TAI
c     (International Atomic Time) time system.
c
c     Inputs:
c     integer year: The standard calendar year, which starts from 1.
c     integer month: The standard calendar month, starting from 1, which
c       corresponds to January.
c     integer dom: The day of month, starting from 1, which accounts for
c       leap years.
c     integer hour: The hour of the day, ranging from 0 - 23.
c     integer min: The minute of the hour, ranging from 0 - 59.
c     real*8 sec: The second of the minute, ranging from 0.0 - <60.0 for
c       minutes not containing a leap second
c      
c     Outputs:
c     real*8 utc_to_tdt: Elapsed seconds since Greenwich noon on January
c       1st 2000 (J2000 Epoch) in the TAI (International Atomic Time)
c       time system.
c                 
      implicit none
c     Inputs
      integer year,month,dom
      integer hour,min
      real*8 sec

c     Internal variables
      real*8 secs,secs_since_epoch
      real*8 timeCor
      integer i
      integer    numleaps
      parameter (numleaps=27)
      
      real*8 leapseconds(numleaps)
      data   leapseconds /-867888000.0d0,-851990400.0d0,-820454400.0d0,
     .     -788918400.0d0,-757382400.0d0,-725760000.0d0,-694224000.0d0,
     .     -662688000.0d0,-631152000.0d0,-583891200.0d0,-552355200.0d0,
     .     -520819200.0d0,-457660800.0d0,-378691200.0d0,-315532800.0d0,
     .     -283996800.0d0,-236736000.0d0,-205200000.0d0,-173664000.0d0,
     .     -126230400.0d0,-78969600.0d0,-31536000.0d0,189388800.0d0,
     .     284083200.0d0,394416000.0d0,489024000.0d0,536544000.0d0/
      
      secs = secs_since_epoch(year,month,dom, hour,min,sec)

      timeCor = 10.0d0
      do i=1, numleaps
         if(secs .ge. leapseconds(i)) timeCor = timeCor+1.0d0
      enddo

      utc_to_tai = secs+timeCor-43200.0d0
      
      return
      end

            
      real*8 function secs_since_epoch(year,month,dayofmonth,
     .     hour,min,sec)
c
c     Computes the number of elapsed seconds since some epoch (in this
c     case Greenwich midnight on January 1st 2000) not accounting for
c     leapseconds.
c     
c     Inputs:
c     integer year: The standard calendar year, which starts from 1.
c     integer month: The standard calendar month, starting from 1, which
c       corresponds to January.
c     integer dom: The day of month, starting from 1, which accounts for
c       leap years.
c     integer hour: The hour of the day, ranging from 0 - 23.
c     integer min: The minute of the hour, ranging from 0 - 59.
c     real*8 sec: The second of the minute, ranging from 0.0 - <60.0 for
c       minutes not containing a leap second
c      
c     Outputs:
c     real*8 secs_since_epoch: Elapsed seconds since some epoch (in this
c       case Greenwich midnight on January 1st 2000) not accounting for
c       leapseconds.
c        
      implicit none
c     Inputs
      integer year,month,dayofmonth, hour,min
      real*8 sec

c     Internal variables      
      integer ndays,days_since_epoch
      integer    epochyear
      parameter (epochyear = 2000)
      
c     First calculate the number of days since the epoch
      ndays=days_since_epoch(year,month,dayofmonth)

c     Calculate the seconds for years after the epoch
      if (year .ge. epochyear) then                  
         secs_since_epoch = dble(ndays)*86400D0+dble(hour)*3600D0+
     .        dble(min)*60D0+dble(sec)
         
         return
         
c     Calculate the seconds for years prior to the epoch
      else 
         secs_since_epoch = dble(ndays-1)*86400D0+dble(hour)*3600D0+
     .        dble(min)*60D0+dble(sec)
         
         return
      end if
      end

      
      integer function days_since_epoch(year,month,dayofmonth)
c
c     Computes the number of elapsed days since some epoch (in this case
c     January 1st 2000). So inputing Jaunary 1st 2000 (2000, 1, 1) would
c     return 0. Likewise Jaunary 2nd 2000 (2000, 1, 2) would return 1.
c     
c     Inputs:
c     integer year: The standard calendar year, which starts from 1.
c     integer month: The standard calendar month, starting from 1, which
c       corresponds to January.
c     integer dom: The day of month, starting from 1, which accounts for
c       leap years.
c      
c     Outputs:
c     integer days_since_epoch: Elapsed days since some epoch (in this
c       case January 1st 2000).
c   
      implicit none
c     Inputs
      integer year,month,dayofmonth

c     Internal variables      
      integer ndays,doy, i
      logical isleapyear
      integer    epochyear
      parameter (epochyear = 2000)

c     First calculate the number of days since the epoch
      ndays=0
      if (year .ge. epochyear) then
         
         do i=epochyear,year-1
            ndays = ndays+365
            if (isleapyear(i)) ndays=ndays+1
         enddo
         
c     Now add the number of days of the current year
         ndays = ndays+doy(year,month,dayofmonth)-1

         days_since_epoch = ndays
         
         return
         
c     Calculate the days for the years prior to the epoch
      else 
         
         do i=year, epochyear-1
            ndays = ndays-365
            if (isleapyear(i)) ndays=ndays-1
         enddo
         
c     Now subtract the number of days of the current year	
         ndays = ndays+doy(year,month,dayofmonth)
         
         days_since_epoch = ndays
         
         return
      end if
      end
      
      
      logical function isleapyear(year)
c     
c     Is this year a leap year? True for yes, false for no. (Gregorian
c     Calendar every 4 years is a leap year, except centennial years not
c     divisible by 400).
c     
c     Inputs:
c     integer year: The standard calendar year, which starts from 1.
c     
c     Outputs:
c     logical isleapyear: True if the year is a leap year false if it is
c       not
c   
      implicit none
c     Inputs
      integer year
      
      isleapyear=((mod(year,4) .eq. 0) .and. (mod(year,100) .ne. 0))
     .     .or.  (mod(year,400) .eq. 0)

      return
      end

      
      integer function doy(year,month,dayofmonth)
c     
c     Computes the day of the year, starting from 1, which corresponds
c     to January 1st
c     
c     Inputs:
c     integer year: The standard calendar year, which starts from 1.
c     integer month: The standard calendar month, starting from 1, which
c       corresponds to January.
c     integer dom: The day of month, starting from 1, which accounts for
c       leap years.
c     
c     Outputs:
c     integer doy: The day of the year, starting from 1
c         
      implicit none
c     Inputs      
      integer year,month,dayofmonth

c     Internal variables      
      logical isleapyear
      integer i,daysInMonth(12),daysInMonthReg(12),daysInMonthLeap(12)
      data daysInMonthReg  /31,28,31,30,31,30,31,31,30,31,30,31/
      data daysInMonthLeap /31,29,31,30,31,30,31,31,30,31,30,31/

c     if it is a leap year, set dayInMonth to the one with the extra day
      daysInMonth = daysInMonthReg
      if (isleapyear(year)) daysInMonth = daysInMonthLeap

c     first, add the days in the month, and then loop through all the
c     months up to but not including the current month and increment
c     the days
      doy = dayofmonth
      do i = 1,month-1
         doy = doy + daysInMonth(i)
      enddo
      
      return
      end
