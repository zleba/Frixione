      subroutine sysdep
c use newver='NEW' for vaxes, 'UNKNOWN' in other machines,
      character * 7 newver
      common/newver/newver
      newver = 'UNKNOWN'
      end

      subroutine delete(fname)
      character * 80 str
      character * (*) fname
      l = len(fname)
      k = 1
      dowhile(fname(k:k).eq.' '.and.k.lt.l)
         k = k+1
      enddo
      dowhile(fname(l:l).eq.' '.and.l.gt.k)
         l = l-1
      enddo
      if(l-k.gt.77) then
         write(*,*) 'delete: filename > 77 chars not allowed'
         stop
      endif
      if(l.eq.k) then
         write(*,*) 'delete: void filename'
         stop
      endif
      str(1:3) = 'rm '
      str(4:) = fname(k:l)
      call system(str)
      end

      subroutine time(ctime)
      character * (*) ctime
      call _time(ctime)
      end

      subroutine idate(iday,imon,iyear)
      call _idate(iday,imon,iyear)
      end
