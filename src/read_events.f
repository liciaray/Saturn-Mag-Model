      subroutine read_events(eventsfile,event_year,
     .     event_month, event_dom, event_hour, event_min,
     .     event_x,event_y,event_z,final_count,noEvents)
      implicit none
c     Inputs
      character*(*) eventsfile

c     Constants
      integer    noEvents
c      parameter (noEvents = 500)
      integer    lun
      parameter (lun = 34)
      
c     Outputs
      integer event_year(noEvents)
      integer event_month(noEvents)
      integer event_dom(noEvents)
      integer event_hour(noEvents)
      integer event_min(noEvents)
      real*8 event_x(noEvents)
      real*8 event_y(noEvents)
      real*8 event_z(noEvents)

c     Internal variables      
      integer sc_year,sc_month,sc_dom, sc_hour,sc_min
      real*8 sc_x,sc_y,sc_z
      real*8 sec
      integer EOF,count,final_count
      character*(256) c
      
      open (unit=lun,file=eventsfile)

c     skip past the header
      read(lun, '(A256)') c

      count = 0
      
      DO WHILE(1.EQ.1)
         count = count+1
         read(lun, FMT=
     .        '(I4,A1,I2,A1,I2,A1,I2,A1,I2,A1,F13.8,A1,F13.8,A1,F13.8)',
     .        iostat=eof),
     .        sc_year,c,sc_month,c,sc_dom, c,
     .        sc_hour,c,sc_min,c,sc_x,c,
     .        sc_y,c,sc_z

         event_year(count) = sc_year
         event_month(count) = sc_month
         event_dom(count) = sc_dom
         event_hour(count) = sc_hour
         event_min(count) = sc_min
         event_x(count) = sc_x
         event_y(count) = sc_y
         event_z(count) = sc_z

         IF (EOF<0) EXIT
      ENDDO            
      
      close(lun)

      final_count=count

      return
      END  


