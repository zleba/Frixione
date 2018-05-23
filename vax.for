      subroutine sysdep
c use newver='NEW' for vaxes, 'UNKNOWN' in other machines,
      character * 7 newver
      common/newver/newver
      newver = 'NEW    '
      end

      subroutine delete(fname)
      character * (*) fname
      ires = lib$delete_file(fname//';')
      if(ires.ne.1) then
          write(*,*) ' delete: cannot delete ',fname
      endif
      end
