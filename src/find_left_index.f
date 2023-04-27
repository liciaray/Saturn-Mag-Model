
      
      integer function find_left_index (t, ts, size)
      real*8 t
      real*8 ts(*)
      integer size
      integer index,lstled

      real*8 t0,t1,tmm,tm,tmp
      
      integer l,r,m
c
c     Finds the position of a double (e.g., time) within an array of
c     doubles (e.g., times). Specifically, this returns the index within
c     the supplied array that is left of (less than) or equal to the
c     supplied double. This is particularly useful for making
c     discretized values continuous functions functions, for example, as
c     is done while interpolating. If the supplied value is less than
c     the first value in the array, zero is returned. If the supplied
c     value is greater than or equal to the last value in the array,
c     the array size is returned. Note, this function is similar to
c     NAIF/SPICE's LSTLED function, but has been independently
c     implemented, so it may not always give the same results.
c
c     Inputs:
c     real*8 t: A double value (e.g., time).
c     real*8 ts: An array of monotonically increasing values
c       (e.g., times).
c     integer size: The size (number of elements) of the array.
c      
c     Outputs:
c     integer find_left_index: The index within the supplied array that
c       is left of (less than) or equal to the suppled value. 
c

c     the first and last values in the array
      t0 = ts(1)
      t1 = ts(size)

c     First, let's check if the t is less than the first value in the
c     array or greater than or equal to the last value.
      if (t .lt. t0) then
         find_last_index_lte = 0
      else if (t .ge. t1) then
         find_last_index_lte = size
      else 

c     Now we will start a binary search algorithm. To begin, set the
c     left and the right values as the first and last values in the
c     array.
         l = 1
         r = size

c     Keep looping until we converge on the solution.
         do while (l .le. r)
            
c     the midpoint between the left and the right index
            m = floor((l+r)/2.0d0)
c            m = ceiling((l+r)/2.0d0)

c     determine if our value is left or right of the midpoint
c     whichever is the case, bisect the l or r.
            if (ts(m) .lt. t) then
               l = m+1
            else if (ts(m) .gt. t) then
               r = m-1
c     if neither of the above if statements are true, then we have
c     reached convergence and it is time to exit the loop
            else
               exit
            endif
         enddo

c     The above binary search alogorithm assumes t is an element of ts.
c     The result is that m may be shifted by 1 to the left or the right,
c     so add in the following checks and adjust m as needed.
         tmm = ts(m-1)          ! the value to the left of m
         tm = ts(m)             ! the value at m
         tmp = ts(m+1)          ! the value to the right of m

c     if the value to the left of m is > t, we need to subtract 1 from m
         if (tm .gt. t) then
            m = m-1
         endif
c     if the value to the right of m is <= t, we need to add 1 to m
         if (tmp .le. t) then
            m = m+1
         endif

c     set the return value to m
         find_last_index_lte = m
         
      end if
     
      return
      end
