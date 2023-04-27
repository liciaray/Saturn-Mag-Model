
      
      subroutine read_pdyn(pdynfile, year,month,dom, hour,min,sec)
      implicit none
c     Inputs
      character*(*) pdynfile
      integer year,month,dom, hour,min
      real*8 sec
      
c     Outputs
      real*8 pdyn
      
c     Internal variables
      
      real*8 utc_to_tdt         ! function that converts UTC to TDT
!     variables used for interpolating pDyn
      real*8 timeLeft,timeRight,valueLeft,valueRight,slope
      integer n,index,index2,find_left_index
      real*8 time
      
c     Constants
      integer    fileSize
      parameter (fileSize = 736849)

c     Outputs
      real*8 times(fileSize)
      real*8 pDyns(fileSize)
      
      save times
      save pDyns

      
      entry store_pdyn(pdynfile)
      call read_pdyn_file(pdynfile, times, pDyns)
      return

      
      entry eval_pdyn(year,month,dom, hour,min,sec, pdyn)
c
c     Evaluates the dynamic pressure for the supplied UTC date and time.
c     The value is liearnly interpolated from the values in the file.
c     NOTE: The subroutine store_pdyn MUST be called first to load the
c     dynamic pressure values from the file.
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
c       minutes not containing a leap second.
c      
c     Outputs:
c     real*8 pdyn: The interpolated dynamic pressure.
c
      
c     converts a UTC date and time into a double
      time = utc_to_tdt(year,month,dom, hour,min,sec)

c     find the index to the left of the supplied time
      index = find_left_index(time,times,filesize)
      
c     retrieve the bracketing times
      timeLeft = times(index)
      timeRight = times(index+1)

c     retrieve the bracketing dynamic pressures      
      valueLeft = pDyns(index)
      valueRight = pDyns(index+1)

c     evaluate the slope from these two points
      slope = (valueRight-valueLeft)/(timeRight-timeLeft)

c     now interpolate the function by evaluating the line
      pdyn = slope*(time-timeLeft)+valueLeft
      
      return

      end


      subroutine read_pdyn_file(pdynfile, times,pDyns)
      implicit none
c     Inputs
      character*(*) pdynfile

c     Constants
      integer    fileSize
      parameter (fileSize = 736849)
      integer    lun
      parameter (lun = 23)
      
c     Outputs
      real*8 times(fileSize)
      real*8 pDyns(fileSize)

c     Internal variables      
      integer year,month,dom, hour,min
      real*8 pDyn
      real*8 sec,utc_to_tdt
      integer EOF,count
      character*(256) c
      
      open (unit=lun,file=pdynfile)

c     skip past the header
      read(lun, '(A256)') c
      read(lun, '(A256)') c
      read(lun, '(A256)') c

      count = 0
      
      DO WHILE(1.EQ.1)
         count = count+1
         read(lun, FMT='(I4,A1,I2,A1,I2,A1,I2,A1,I2,A7,F11.9)',
     .        iostat=eof)
     .        year,c,month,c,dom, c, hour,c,min,c,pDyn

         sec = 0.0d0
         times(count) = utc_to_tdt(year,month,dom, hour,min,sec)
         pDyns(count) = pDyn
         
         IF (EOF<0) EXIT
      ENDDO      
      
      close(lun)

      return
      END  


