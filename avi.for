      subroutine sysdep
c use newver='NEW' for vaxes, 'UNKNOWN' in other machines,
      character * 7 newver
      common/newver/newver
      newver = 'UNKNOWN'
      end

      subroutine delete(fname)
      character * (*) fname
      integer unlink
      j=1
      l = len(fname)
      dowhile(fname(j:j).eq.' '.and.j.lt.l)
         j = j+1
      enddo
      k = l
      dowhile(fname(k:k).eq.' '.and.k.gt.0)
         k = k-1
      enddo
      if(k.lt.l.and.k.gt.j) then
      fname(k+1:k+1)=char(0)
      else
         write(*,*)'delete: error with fname'
         stop
      endif
      ir = unlink(fname(1:1))
      if(ir.ne.0) then
         write(*,*)' delete: cannot delete ','"'//fname(j:k)//'"'
      endif
      end
