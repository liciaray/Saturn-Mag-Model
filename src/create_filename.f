      SUBROUTINE create_filename(year_txt,month_txt,
     .     dom_txt,hour_txt,min_txt,
     .     geo_file)

      implicit none
      character(len=60) geo_file
      integer year_txt, month_txt
      integer dom_txt
      integer hour_txt, min_txt

      write(geo_file, '(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A)')
     .     '../output/event_', year_txt,'-', 
     .     month_txt, '-', dom_txt, 'T', 
     .     hour_txt, 'h', min_txt,'m.txt'

      end


    
