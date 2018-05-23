c-
c Integration-histogramming package using vegas.
c Calls:
c      subroutine delete(fname)         ! deletes the file fname
c      character * (*) fname
c
c      subroutine time(ctime)           ! returns a string with the time
c      character * 8 ctime
c
c      subroutine date(cdate)           ! returns a string with the date
c      character * 9 cdate
c

      subroutine strnum(string,num)
c- writes the number num on the string string starting at the blank
c- following the last non-blank character
      character * (*) string
      character * 20 tmp
      l = len(string)
      write(tmp,'(i15)')num
      j=1
      dowhile(tmp(j:j).eq.' ')
        j=j+1
      enddo
      ipos = istrl(string)
      ito = ipos+1+(15-j)
      if(ito.gt.l) then
         write(*,*)'error, string too short'
         write(*,*) string
         stop
      endif
      string(ipos+1:ito)=tmp(j:)
      end

      function istrl(string)
c returns the position of the last non-blank character in string
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end

      subroutine strcat(str1,str2,str)
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end

      subroutine
     # integrate(init,sig,str,n,inew,idim,icalls,avgi,sd,chi2a,isave)
      implicit real * 8 (a-h,o-z)
      external init,sig
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/newver/newver
      character * (*) str
      character * 7 statf,newver
      character * 25 fname
      logical isave
c- if inew = 0 ignore
      if(inew.eq.0) return
c- add str to list of save files
      if(isave) then
         fname = str
         call addfil(fname)
      endif
c- initialize histograms
      call init
c
      if(inew.ge.2) then
         call resume(str,init,'IN')
c also fixes the value of it
         call vegas4(sig,avgi,sd,chi2a)
      endif
c
      if(inew.eq.10) then
c- save the file '.sv2' in ASCII format
         call strcat(str,'.sv2',fname)
c- ndim = idim for upward compatibility with old restart files
         ndim = idim
         call savea(fname,newver)
         return
      endif
c
      if(inew.eq.11) then
c- save the file '.sv1' in standard format
         call strcat(str,'.sv1',fname)
         call save(fname,newver)
         return
      endif
c
c Initialize parameters for Vegas
      do j=1,10
          xl(j) = 0
          xu(j) = 1
      enddo
      acc = -1
      ndim = idim
      ncall = icalls
      itmx = 1
c
c for negative n abs(n) is the number of iterations to reach,
c for positive n, n is the number of iterations to perform
      if(n.lt.0) then
         if(inew.eq.1.or.inew.eq.3) then
            nsteps = abs(n)
         elseif(inew.eq.2) then
            nsteps = abs(n) - it
         endif
      else
         nsteps = n
      endif
c create new versions of save file
      statf = newver
c----------------------------------------------------------
      do j=1,nsteps
         if(j.eq.1) then
            if(inew.eq.1) then
                  call vegas(sig,avgi,sd,chi2a)
            elseif(inew.eq.2) then
                  call vegas2(sig,avgi,sd,chi2a)
            elseif(inew.eq.3) then
                  call init
                  call vegas1(sig,avgi,sd,chi2a)
            else
                  write(*,*) 'irregular inew flag ',inew
                  stop
            endif
         else
            call vegas3(sig,avgi,sd,chi2a)
         endif
         call accum
         if(isave) then
            call strcat(str,'.sv2',fname)
            call save(fname,statf)
            call strcat(str,'.sv1',fname)
            call save(fname,statf)
            call strcat(str,'.sv2',fname)
            call delete(fname)
            statf = 'UNKNOWN'
         endif
      enddo
      return
      entry integrms1
c # of iterations message
      write(*,*)
      write(*,*)'if you enter a number of iterations n>0'
      write(*,*)'the program will execute n iterations'
      write(*,*)'if n<0, it will execute iterations until'
      write(*,*)'total number of iterations performed = abs(n)'
      write(*,*)'(it makes a difference only if in restart mode)'
      write(*,*)
      return
      entry integrms2
c inew values
      write(*,*)
      write(*,*) 'enter 0 to exclude,1 for new run, 2 to restart'
      write(*,*) '3 to restart keeping the grid but reset histo''s '
      write(*,*) '10 to produce ASCII save files from standard ones'
      write(*,*) '11 to produce standard save files from ASCII ones'
      write(*,*) '(10/11 used to transport save files across machines)'
      write(*,*)
      return
      entry nopr(nnn)
      nprn = nnn
      end

      subroutine resume(string,init,fl)
      character * (*) string, fl *2
      character * 25 fname
      if(fl.eq.'IN') call init
      call strcat(string,'.sv1',fname)
      call restart(fname,itmp)
      if(itmp.eq.1.and.fl.eq.'NO') goto 99
      if(itmp.eq.0) return
c- retry with spare save file
      if(fl.eq.'IN') call init
      call strcat(string,'.sv2',fname)
      call restart(fname,itmp)
      if(itmp.eq.1.and.fl.eq.'NO') goto 99
      if(itmp.eq.0) return
c- retry with ASCII spare save file
      if(fl.eq.'IN') call init
      call strcat(string,'.sv2',fname)
      call restara(fname,itmp)
      if(itmp.eq.1.and.fl.eq.'NO') goto 99
      if(itmp.eq.0) return
c- give up.
      write(*,*) 'unusable restart file'
      stop
99    continue
c In this case the accumulated values in the histograms common block
c could have been altered in the attempt to read the save file.
c Abort. (if itmp=2 this is not the case)
      write(*,*) 'corrupted file',fname,'.'
      write(*,'(a,/)')
     # 'The attempt to read the corrupted file may have corrupted',
     # 'the histogram common block. Delete the corrupted file, or',
     # 'rerun one more iteration, or run with nit=11.'
      stop
      end

      subroutine savea(name,statf)
      PARAMETER (NMB=200)
      character * 80 runstr,string
      common/histo/hist(nmb,100),xhis(nmb,100),hdel(nmb),hmin(nmb)
     &,hmax(nmb),uscore(nmb),oscore(nmb)
     &,nbin(nmb),ihis(nmb,100),iuscore(nmb),ioscore(nmb)
     &,ient(nmb),havg(nmb),hint(nmb),hsig(nmb),book(nmb),title(nmb)
     &,nhist
      character title*50,book*3
      character * 3 books(nmb)
      real * 8 xl,xu,acc
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      real * 8 xi,si,si2,swgt,schi
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/seed/num1,num2
      common/runstr/runstr
      character * (*) name, statf
c
c Salva gli istogrammi in uso (marcati 'YES' in book(j))
c e tutte le informazioni per vegas, il seme del numero
c casuale, etc.
c Marchia la fine del file con i numeri acos(-1.),exp(1.),log(10.),
c
      open(unit=97,file=name,status=statf)
101   format(a)
102   format(25a)
103   format(6d13.7)
104   format(6i12)
105   format(4d20.14)
c version of save file
      string = 'SAVEA,VER 8-3-92'
      write(97,101) string
c runstring
      write(97,101)  runstr
c booked histograms
      write(97,104)  nmb
      write(97,102)  (book(j),j=1,nmb)
c Mbook block
      do j=1,nmb
         if(book(j).eq.'YES') then
            write(97,103)
     &      (hist(j,k),k=1,nbin(j)),uscore(j),oscore(j)
            write(97,104)
     &      (ihis(j,k),k=1,nbin(j)),iuscore(j),ioscore(j),ient(j)
         endif
      enddo
c Vegas block
      write(97,104)
     & ndo,it,num1,num2,ndim
      write(97,105)
     & ((xi(j,k),j=1,ndo),k=1,ndim),si,si2,swgt,schi
c
      a=acos(-1.)
      b=exp(1.)
      c=log(10.)
c end mark
      write(97,103)
     & a,b,c
      close(97)
      return
c
c Loads the histograms found in the (ascii) save file into the common block.
c Previous values are overwritten.
c Does not modify histograms not present in the file.
c If all goes well, it riturns iflag= 0.
c If there are errors (e.g. IO errors, the mark at the end of the file
c is not found, etc.) and in the input operation the common blocks
c may have been corrupted, returns iflag=1.
c If there are errors, but no common block has been altered, returns iflag=2.
c
      entry restara(name,iflag)
      open(unit=97,file=name,status='OLD',err=2)
c version of save file
      string = ' '
      read(97,101,err=2,end=2) string
      if(string.ne.'SAVEA,VER 8-3-92') goto 2
c run string
      string = ' '
      read(97,101,err=2,end=2) string
c marked histograms
      read(unit=97,fmt=104,err=2,end=2) nhs
      if(nhs.gt.nmb) then
         write(*,*) 'Error: save files had nmb=',nhs
         write(*,*) 'Running program has nmb=',nmb,'<',nhs
         write(*,*) 'Make nmb >= ',nhs,' in running program'
         stop
      endif
      read(unit=97,fmt=102,err=1,end=1) (books(j),j=1,nhs)
c Mbook block
      do j=1,nhs
         if(books(j).eq.'YES') then
            read(unit=97,fmt=103,err=1,end=1)
     &      (hist(j,k),k=1,nbin(j)),uscore(j),oscore(j)
            read(unit=97,fmt=104,err=1,end=1)
     &      (ihis(j,k),k=1,nbin(j)),iuscore(j),ioscore(j),ient(j)
         endif
      enddo
c vegas block
      read(unit=97,fmt=104,err=1,end=1)
     & ndo,it,num1,num2,ndim
      read(unit=97,fmt=105,err=1,end=1)
     & ((xi(J,k),j=1,ndo),k=1,ndim),si,si2,swgt,schi
c end mark
      read(unit=97,fmt=103,err=1,end=1)
     & a0,b0,c0
c
      tiny = .1e-5
      a = abs(a0-acos(-1.))
      b = abs(b0-exp(1.))
      c = abs(c0-log(10.))
      if( a.gt.tiny .or. b.gt.tiny .or. c.gt.tiny ) goto 1
      close(unit=97,iostat=ios)
      iflag = 0
      if(string.ne.runstr) then
         write(6,*)
     #  'You are using save files belonging to a different run,'
         write(6,*) string
         write(6,*) 'instead of'
         write(6,*) runstr
      endif
      write(6,*) string
      return
 1    continue
      iflag = 1
      goto 3
 2    continue
      iflag = 2
      goto 3
 3    continue
      close(unit=97,iostat=ios)
      end

      subroutine save(name,statf)
      PARAMETER (NMB=200)
      character * 80 runstr,string
      common/histo/hist(nmb,100),xhis(nmb,100),hdel(nmb),hmin(nmb)
     &,hmax(nmb),uscore(nmb),oscore(nmb)
     &,nbin(nmb),ihis(nmb,100),iuscore(nmb),ioscore(nmb)
     &,ient(nmb),havg(nmb),hint(nmb),hsig(nmb),book(nmb),title(nmb)
     &,nhist
      character title*50,book*3
      character * 3 books(nmb)
      real * 8 xl,xu,acc
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      real * 8 xi,si,si2,swgt,schi
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/seed/num1,num2
      common/runstr/runstr
      character * (*) name, statf
c
c Salva gli istogrammi in uso (marcati 'YES' in book(j))
c e tutte le informazioni per vegas, il seme del numero
c casuale, etc.
c Marchia la fine del file con i numeri acos(-1.),exp(1.),log(10.),
c
      open(unit=97,file=name,form='UNFORMATTED',status=statf)
c version of save program
      string = 'SAVE,VER 8-3-92'
      write(97) string
c run string
      write(97)  runstr
c Mbook Block
      write(97) nmb
      write(97)  (book(j),j=1,nmb)
      do j=1,nmb
         if(book(j).eq.'YES') then
           write(97)
     &     (hist(j,k),k=1,nbin(j))
     &     ,uscore(j),oscore(j)
     &     ,(ihis(j,k),k=1,nbin(j))
     &     ,iuscore(j),ioscore(j)
     &     ,ient(j)
         endif
      enddo
c Vegas Block
      write(97)
     & ndo,it,num1,num2,ndim
      write(97)
     & ((xi(J,k),j=1,ndo),k=1,ndim),si,si2,swgt,schi
c
      write(97)   num1,num2
      a=acos(-1.)
      b=exp(1.)
      c=log(10.)
      write(97) a,b,c
      close(97)
      return
c
c Loads the histograms found in the file into the common block.
c Previous values are overwritten.
c Does not modify histograms not present in the file.
c If all goes well, it riturns iflag= 0.
c If there are errors (e.g. IO errors, the mark at the end of the file
c is not found, etc.) and in the input operation the common blocks
c may have been corrupted, returns iflag=1.
c If there are errors, but no common block has been altered, returns iflag=2.
c
      entry restart(name,iflag)
      open(unit=97,file=name,form='unformatted',status='OLD',err=2)
c version of save file
      string = ' '
      read(unit=97,err=2,end=2) string
      if(string.ne.'SAVE,VER 8-3-92') goto 2
c run string
      string = ' '
      read(unit=97,err=2,end=2) string
c Mbook block
      read(97,err=2,end=2) nhs
      if(nhs.gt.nmb) then
         write(*,*) 'Error: save files had nmb=',nhs
         write(*,*) 'Running program has nmb=',nmb,'<',nhs
         write(*,*) 'Make nmb >= ',nhs,' in running program'
         stop
      endif
      read(97,err=1,end=1) (books(j),j=1,nhs)
      do j=1,nhs
         if(books(j).eq.'YES') then
            read(97,err=1,end=1)
     &      (hist(j,k),k=1,nbin(j))
     &      ,uscore(j),oscore(j)
     &      ,(ihis(j,k),k=1,nbin(j))
     &      ,iuscore(j),ioscore(j)
     &      ,ient(j)
         endif
      enddo
c Vegas Block
      read(97,err=1,end=1)
     & ndo,it,num1,num2,ndim
      read(97,err=1,end=1)
     & ((xi(J,k),j=1,ndo),k=1,ndim),si,si2,swgt,schi
c
      read(97,err=1,end=1)   num1,num2
      a0 = 0
      b0 = 0
      c0 = 0
      read(97,err=1,end=1)   a0,b0,c0
      a=acos(-1.)
      b=exp(1.)
      c=log(10.)
      if( a.ne.a0 .or. b.ne.b0 .or. c.ne.c0 ) goto 1
      close(97)
      iflag = 0
      if(string.ne.runstr) then
         write(6,*)
     #   'You are using save files belonging to a different run,'
         write(6,*) string
         write(6,*) 'instead of'
         write(6,*) runstr
      endif
      write(6,*) string
      return
 1    continue
      iflag = 1
      goto 3
 2    continue
      iflag = 2
      goto 3
 3    continue
      close(unit=97,iostat=ios)
      end

      subroutine bookup(n,string,del,xl,xu)
      character * (*) string
c
c Per ogni istogramma da fare, ne sono richiesti quattro
c In n si accumulano i valori in outfun.
c A ogni iterazione di vegas l'istogramma n viene sommato a n+1,
c e il quadrato del suo valore viene sommato a n+2,
c L'istogramma n viene anche usato alla fine per combinare
c i totali dei vari contributi sig0,sig2, etc. mentre
c l'istogramma n+3 viene usato per combinare gli
c errori dei vari contributi sig0,sig2, etc.
c Non si vuole che n e n+3 vengano salvati o riesumati.
c Cambiando il tag in N,N+3, questo non avviene (si guardi
c in mbook e save/restart e anche mclear in questo programma.
c
      call mbook(n,  string,del,xl,xu)
      call mbook(n+1,'tmp ',del,xl,xu)
      call mbook(n+2,'tmp square',del,xl,xu)
      call mbook(n+3,'error ',del,xl,xu)
      call puttag(n,'YST')
      call puttag(n+3,'NST')
      return
      end

      subroutine accum
      PARAMETER (NMB=200)
      implicit real * 4 (a-h,o-z)
      character * 3 tag
c
c     Accumula i valori e i valori al quadrato per l'analisi statistica,
c     e svuota l'istogramma di accumulo.
c
      do j=1,nmb-3
         call gettag(j,tag)
         if(tag.eq.'YST') then
             call mopera(j,'A',j+1,j+2,dum,dum)
         endif
      enddo
      return
      end

      subroutine mclear
      implicit real * 8 (a-h,o-z)
      character * 25 files(20), filn
      data nfil/0/
c MODIFICA
      common/mttini/ini
      ini=0
c
c sum up files
      call sumfil(files,nfil)
      nfil=0
      return
      entry addfil(filn)
c adds file filn to the list of save files.
      nfil = nfil+1
      if(nfil.gt.20) then
         write(*,*) 'mclear: too many files.'
         stop
      endif
      files(nfil) = filn
      end

      subroutine sumfil(files,nfil)
      PARAMETER (NMB=200)
      implicit real * 8 (a-h,o-z)
      external dummyinit
      character * 3 tag
      character * 25 files(*)
c
      do j=1,nfil
         call resume(files(j),dummyinit,'NO')
c
c completa l'analisi statistica
c
         do k=1,nmb-3
            call gettag(k,tag)
            if(tag.eq.'YST') then
               call addup(k)
            endif
         enddo
      enddo
      end

      subroutine dummyinit()
      return
      end

      subroutine addup(j)
c
c accumula j+1 riscalato in j e pone in j+2 la stima dell'errore
c
      call mopera(j+1,'E',j+2,j,dum,dum)
c
c accumula l'errore in quadratura
c
      call mopera(j+2,'Q',j+3,j+3,dum,dum)
      end

       block data vegas0
       implicit  real*8 (a-h,o-z)
       common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
       common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
       data ncall/10000/,itmx/10/,nprn/ 1000/,acc/-1.d0/,
     1    xl/1.d-3, 1.d-3,1.d-3,1.d-3,1.d-3,1.d-3,1.d-3,1.d0,1.d0,1.d0/,
     2   xu/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/
       data xi/500*1.d0/
       end
C
C
C      NCALL IS THE NUMBER OF CALLS TO VEGAS.
C      NPRN >  0 VEGAS PRINTS THE RESULTS OF EACH ITERATION.
C      NPRN 0 VEGAS PRINTS NOTHING.
C      NPRN < 0 VEGAS PRINTS ALL.
C      XL(I) IS LOWER INTEGRATION LIMIT ON I TH AXIS.
C       XU(I) IS UPPER INTEGRATION LIMIT ON I THE AXIS.
c
         subroutine vegas(fxn,avgi,sd,chi2a)
c
c   routine performs n dim monte carlo inte
c written by p lepage
c
         implicit real*8 (a-h,o-z)
       implicit integer*4 (i-n)
         common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
         common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
       COMMON/SEED/NUM,NUM2
         dimension d(50,10),di(50,10),xin(50),r(50),dx(10),dt(10),
     1   x(10),kg(10),ia(10)
       dimension RAND(10)
         data ndmx/50/,alph/1.5d0/,one/1.d0/,mds/1/
c
       NUM=1
c       NUM2 e' irrilevante
         ndo=1
         do 1 j=1,ndim
 1       xi(1,j)=one
c
         entry vegas1(fxn,avgi,sd,chi2a)
c       initialises  cumulative  variables but not grid
         it=0
         si=0.
         si2=si
         swgt=si
         schi=si
c
         entry vegas2(fxn,avgi,sd,chi2a)
c        no initialisation
         nd=ndmx
         ng=1
         if(mds.eq.0)go to 2
         ng=(ncall/2.)**(1./ndim)
         mds=1
         if((2*ng-ndmx).lt.0)go to 2
         mds=-1
         npg=ng/ndmx+1
         nd=ng/npg
         ng=npg*nd
 2       k=ng**ndim
         npg=ncall/k
         if(npg.lt.2)npg=2
         calls=npg*k
         dxg=one/ng
         dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
         xnd=nd
         ndm=nd-1
         dxg=dxg*xnd
         xjac=one/calls
         do 3 j=1,ndim
         dx(j)=xu(j)-xl(j)
 3       xjac=xjac*dx(j)
c
c    rebin preserving bin density
c
         if(nd.eq.ndo)go to 8
         rc=ndo/xnd
         do 7 J=1,ndim
         k=0
         xn=0.
         dr=xn
         i=k
 4       k=k+1
         dr=dr+one
         xo=xn
         xn=xi(k,j)
 5       if(rc.gt.dr)go to 4
         i=i+1
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr
         if(i.lt.ndm)go to 5
         do 6 i=1,ndm
 6       xi(i,j)=xin(i)
 7       xi(nd,j)=one
         ndo=nd
c
 8       if(nprn.ne.0)write(6,200)ndim,calls,it,itmx,acc
     1   ,mds,nd,(xl(j),xu(j),j=1,ndim)
c
         entry vegas3(fxn,avgi,sd,chi2a)
c         main integration loop
 9       it=it+1
          ti=0.
         tsi=ti
         do 10 j=1,ndim
         kg(j)=1
         do 10 i=1,nd
         d(i,j)=ti
 10      di(i,j)=ti
c
 11      fb=0.
         f2b=fb
         k=0
 12      k=k+1
       call randa(ndim,rand)
         wgt=xjac
         do 15 j=1,ndim
         xn=(kg(j)-rand(j))*dxg+one
         ia(j)=xn
         if(ia(j).gt.1)go to 13
         xo=xi(ia(j),j)
         rc=(xn-ia(j))*xo
         go to 14
13       xO=xi(ia(j),j)-xi(ia(j)-1,j)
         rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
 14      x(j)=xl(j)+rc*dx(j)
 15      wgt=wgt*xo*xnd
c
         f=wgt
         f=f*fxn(x,wgt)
         f2=f*f
         fb=fb+f
         f2b=f2b+f2
         do 16 j=1,ndim
         di(ia(j),j)=di(ia(j),j)+f
 16      if(mds.ge.0)d(ia(j),J)=d(ia(j),J)+f2
         if(k.lt.npg) go to 12
c
888    FORMAT(1X,'F',G14.6,'F2',G14.6,'FB',G14.6,'F2B',G14.6)
         f2b= sqrt(f2b*      NPG)
         f2b=(f2b-fb)*(f2b+fb)
1661   FORMAT(1X,'F2B',G14.6,'NPG',  I10)
         ti=ti+fb
         tsi=tsi+f2b
33     FORMAT(1X,'TSI',G14.6,'F2B',G14.6)
         if(mds.ge.0)go to 18
         do 17 j=1,ndim
 17      d(ia(j),j)=d(ia(j),j)+f2b
 18      k=ndim
 19      kg(k)=mod(kg(k),ng)+1
         if(kg(k).ne.1)go to 11
         k=k-1
         if(k.gt.0)go to 19
c
c final results for this iteration
c
        tsi=tsi*dv2g
        ti2=ti*ti
88     format(1x,'tsi',g14.6)
       if(tsi.eq.0.d0)then 
         rpp=1.d15
       else
         rpp=abs(ti2/tsi)
       endif
       if(rpp.lt.1.d14) then
        wgt=ti2/tsi
        si=si+ti*wgt
        si2=si2+ti2
        swgt=swgt+wgt
        schi=schi+ti2*wgt
995    FORMAT(1X,'SWGT',G14.6,'SI2',G14.6)
        avgi=si/swgt
        sd=swgt*it/si2
        chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999)
        sd=dsqrt(one/sd)
       else
        write(*,*) '***VEGAS WARNING***: zero error!'
        write(*,*) ' we guess that integral is exact '
        avgi=ti
        si=ti
        si2=0.d0
        sd=0
        chi2a=-1
        return
       endif
c
        if(nprn.eq.0)go to 21
        tsi=dsqrt(tsi)
        write(6,201)it,ti,tsi,avgi,sd,chi2a
        if(nprn.ge.0)go to 21
        do 20 j=1,ndim
 20     write(6,202) j,(xi(i,j),di(i,j),d(i,j),i=1,nd)
c
c      refine grid
c
 21     do 23 j=1,ndim
        xo=d(1,j)
        xn=d(2,j)
        d(1,j)=(xo+xn)/2.
        dt(j)=d(1,j)
        do 22 i=2,ndm
        d(i,j)=xo+xn
        xo=xn
        xn=d(i+1,j)
        d(i,j)=(d(i,j)+xn)/3.
 22     dt(j)=dt(j)+d(i,j)
        d(nd,j)=(xn+xo)/2.
 23     dt(j)=dt(j)+d(nd,j)
c
        do 28 j=1,ndim
        rc=0.
        do 24 i=1,nd
        r(i)=0.
        if(d(i,j).le.0.)go to 24
        xo=dt(j)/d(i,j)
        r(i)=((xo-one)/xo/dlog(xo))**alph
 24     rc=rc+r(i)
        rc=rc/xnd
        k=0
        xn=0.
        dr=xn
        i=k
 25     k=k+1
        dr=dr+r(k)
        xo=xn
        xn=xi(k,j)
 26     if(rc.gt.dr)go to 25
        i=i+1
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
        if(i.lt.ndm)go to 26
        do 27 i=1,ndm
 27     xi(i,j)=xin(i)
 28     xi(nd,j)=one
c
        if(it.lt.itmx.and.acc*dabs(avgi).lt.sd)go to 9
 200    format(1X,'0input parameters for vegas:  ndim=',i3,
     1  '   ncall=',f8.0/28x,'  it=',i5,'    itmx=',i5/28x,
     2  '  acc=',g9.3/28x,'  mds=',i3,'     nd=',i4/28x,
     3  '  (xl,xu)=',(t40,'( ',g12.6,' , ',g12.6,' )'))
 201    format(///' integration by vegas' / '0iteration no.',i5,
     1  ':  integral=',g14.8/21x,'std dev =',g10.4 /
     2  ' accumulated results:   integral=',g14.8/
     3  24x,'std dev =',g10.4 / 24x,'chi**2 per it''n =',g10.4)
 202    format(1X,'0data for axis',i2,/,' ',6x,'x',7x,'  delt i ',
     1  2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  '
     2  ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/,
     3  (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))
        return
        entry vegas4(fxn,avgi,sd,chi2a)
        if(si2.eq.0.d0)then
          avgi=si
          sd=0
          chi2a=-1
        else
          avgi=si/swgt
          sd=swgt*it/si2
          chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999)
          sd=dsqrt(one/sd)
        endif
        if(nprn.ne.0) write(6,201)it,0.d0,0.d0,avgi,sd,chi2a
        return
        end

c        subroutine save(ndim)
c        implicit real*8 (a-h,o-z)
c       implicit integer*4 (i-n)
c        common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
c
c      stores vegas data   (unit 7) for later initialisation
c
c        write(7,200) ndo,it,si,si2,swgt,schi,
c     1       ((xi(i,j),i=1,ndo),j=1,ndim)
c        return
c        entry restr(ndim)
c
c         enters initialisation data for vegas
c
c        read(7,200) ndo,it,si,si2,swgt,schi,
c     1    ((xi(i,j),i= 1,ndo),j=1,ndim)
c 200    format(2i8,4z16/(5z16))
c        return
c        end

      
      FUNCTION RANDOM(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION MINV,RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      RANDOM = SEED*MINV
      END

      subroutine randa(n,rand)
      implicit double precision (a-h,o-z)
      COMMON/SEED/NUM,NUM2
      common/caso/caso(5)
      dimension rand(10)
      do 1 i=1,n
      rand(i)=random(NUM)
1     continue
      do 2 i=1,5
      caso(i)=random(NUM)
2     continue
      return
      end

C**********************************************************************
C    SIMPLE HISTOGRAMMING PACKAGE --  SIMPLIFIED VERSION OF HBOOK
C    BY Michelangelo Mangano    NOVEMBER 1988
C    LAST REVISED NOVEMBER 9, 1988
C    LAST REVISED JUNE 12, 1989  (ADD SCATTER PLOTS)
C    LAST REVISED oct 1990 (Add multi-plots on one page, routines MULTITOP,
C                               MTFILL,...)
C**********************************************************************
C
C Fills up to 100 histograms with up to 100 bins.
C Gives a data file (to be specified in the calling program by assigning
C a file name to unit 98) and a topdrawer file (to be specified in the
C calling program by assigning a file name to unit 99).
C
C INITIALIZATION:
C Call once INIHIST; this just resets a few counters and logicals
C Call MBOOK(N,'TITLE',DEL,XMIN,XMAX) for each histogram to be booked.
C N (an integer) is the label of the histogram;
C 'TITLE' is the name of the histogram (no more then 100 characters);
C DEL (real*4) is the bin size;
C XMIN (real*4) is the lower limit of the first bin;
C XMAX (real*4)is the upper limit of the last  bin
C Example:
C      call mbook(2,'pt distribution',1.,10,70)
C This call initializes histogram number 2, called 'pt distribution';
C The bin size will be 1. (possibly GeV, if that's what you want), the
C first bin being  10<x<11. and the last one being 69.<x<70
C
C FILLING:
C When it's time, call MFILL(N,X,Y); this will add Y (real*4) to the bin
C in which X (real*4) happens to be, within histogram N.
C
C PLAYING AROUND:
C At the end of the day you may want to sum, divide, cancel, etc.etc.
C various histograms (bin by bin). Then you call MOPERA(I,'O',J,K,X,Y).
C The 1-character string O can take the following values:
C +  : sums       X*(hist I) with Y*(hist J) and puts the result in hist K;
C -  : subtracts  X*(hist I) with Y*(hist J) and puts the result in hist K;
C *  : multiplies X*(hist I) with Y*(hist J) and puts the result in hist K;
C /  : divides    X*(hist I) with Y*(hist J) and puts the result in hist K;
C F  : multiplies hist I by the factor X, and puts the result in hist K;
C R  : takes the square root of  hist  I, and puts the result in hist K;if
C      the value at a given bin is less than or equal to 0, puts 0 in K
C S  : takes the square      of  hist  I, and puts the result in hist K;
C L  : takes the log_10 of  hist  I, and puts the result in hist K; if the
C      value at a given bin is less than or equal to 0, puts 0 in K
C M  : statistical analysis; if I contains the weights (let's say WGT),
C      J contains variable times weight (F*WGT) and K contains the
C      variable squared times the weight (F**2*WGT), then, after using 'M',
C      J will contain the average value of the variable <F> and K will
C      contain the sigma of the average: sigma=sqrt(<F**2>-<F>**2).
C      If WGT=1. for all the entries, then it is enough to put I=J, and
C      it is not necessary to book a hist with the weights.
C V  : estimates errors for vegas evaluation of differential distributions.
C      Fill I with the values of
C      the functions do integrate times the Vegas weight (fun*wgt); fill
C      J with fun**2*wgt; then K will contain an estimate of the error
C      of the integration. Putting X=1/(#of iterations) performs the
C      avegare over the iterations, and gives the right normalization to
C      the differential distribution, I, and to the errors, K. J stays the same.
C
C FINAL ACCOUNTING:
C Now we can finalize our histograms; MFINAL(N) will calculate the integral
C of the histogram N, the mean value of the X variable and its RMS.
C If we now want to renormalize the hist's, we can call MNORM(N,X), which
C will normalize the integral to X  -- CAUTION: do not call MNORM before
C MFINAL, it will blow up.
C
C OUTPUT:
C To get a .dat file containing the values of the histograms, together with
C some information (like integral, mean values, etc.etc.) call MPRINT(N),
C for each hist N that you want in the .dat file. Before the call to MPRINT
C you want to open unit 98 and give it a name:
C     OPEN(UNIT=98,NAME='NAME.DAT',STATUS='NEW')
C If you want a topdrawer file with a plot of the hist values, call
C MTOP(N,M,'X','Y','SCALE'). The points of the plot will be taken from histogram
C N, the error bars from histogram M. 'SCALE', character*(*), determines
C the scale for y, logarithmic or linear (SCALE=LOG,LIN).
C If you do not want error bars, keep
C a histogram of zeros, or just call a hist that had not been booked.
C X will appear as a 'bottom title', and Y will appear as a 'left title'.
C The top title is by default the name of the histogram itself.
C A little box below the plot will contain some information on the plot
C itself. Before calling MTOP,
C     OPEN(UNIT=99,NAME='NAME.TOP',STATUS='NEW')
C--------------------------------------------------------------------------
C
C  COMMON/HISTO/  Histogram N
C
C   BOOK(N),      Three-letter character-string: 'NO' if histogram was not
C                 Booked, 'YES' otherwise
C   TITLE(N),     Title of the histogram
C
C   HMIN(N),      Min value of x range
C   HMAX(N),      Max value of x range
C   HDEL(N),      Bin width
C   NBIN(N),      Total number of bins
C   USCORE(N),    Total integral of underscores with x < HMIN(N)
C   OSCORE(N),    Total integral of onderscores with x > HMAX(N)
C   IUSCORE(N),   Number of entries with x < HMIN(N)
C   IOSCORE(N),   Number of entries with x > HMAX(N)
C   IENT(N),      Total number of entries within x range HMIN(N)<x<HMAX(N)
C   HINT(N),      Integral of the histogram within HMIN(N)<x<HMAX(N)
C   HAVG(N),      Average value of x, weighted over the x range of the histo
C   HSIG(N),      Quadratic dispersion of x around the average
C   HIST(N,L),    Value of bin L-th
C   XHIS(N,L),    Central x value of bin L-th
C   IHIS(N,L),    Number of entries within bin L-th
C   NHIST         Total number of booked histograms
C


      BLOCKDATA INIHIS0
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE*50,BOOK*3
      CHARACTER TITLE2*50,BOOK2*3
      DATA BOOK/NMB*' NO'/
      DATA BOOK2/10*' NO'/
      END

      SUBROUTINE INIHIST
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE*50,BOOK*3
      CHARACTER TITLE2*50,BOOK2*3
      DO 1, I=1,NMB
   1  BOOK(I)=' NO'
      DO 2 I=1,10
   2  BOOK2(I)=' NO'
      END

      SUBROUTINE MBOOK(N,TIT,DEL,XMIN,XMAX)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
      CHARACTER*(*) TIT
      DATA ONEP/1.00001/
      IF(N.GT.NMB) THEN
         CALL MBKWRN('MBOOK')
         WRITE(*,*) 'MBOOK: no more then ',NMB,' histograms'
         STOP
      ENDIF
      NHIST = MAX(N,NHIST)
      IF(BOOK(N).NE.' NO') THEN
         CALL MBKWRN('MBOOK')
         WRITE(*,*) 'Histogram',N,TITLE(N),' already in use. '
         WRITE(*,*) 'superseded by ',TIT
      ENDIF
      BOOK(N) = 'YES'
      TITLE(N) = ' '//TIT
1     HDEL(N) = DEL
      NBIN(N) = INT(ONEP*(XMAX-XMIN)/DEL)
      IF(NBIN(N).GT.100) THEN
        WRITE(*,*) 'TOO MANY BINS (',NBIN(N),') REQUIRED IN HIST ',N
        WRITE(*,*) 'RE-ENTER BIN SIZE DELTA (OLD BIN = ',DEL,' ):'
        READ(*,*) DEL
        GO TO 1
      ENDIF
      HMIN(N) = XMIN
      HMAX(N) = NBIN(N)*DEL+XMIN
      IF(ABS(HMAX(N)-XMAX)/(XMAX-XMIN).GT.0.1E-5) THEN
         CALL MBKWRN('MBOOK')
         WRITE(*,*)
     #'Histogram ', TIT, ' Change of upper limit:',xmax,'-->',HMAX(N)
      ENDIF
      IENT(N) = 0
      IUSCORE(N) = 0
      IOSCORE(N) = 0
      USCORE(N) = 0
      OSCORE(N) = 0
      HAVG(N) = 0
      HINT(N) = 0
      HSIG(N) = 0
      DO I=1,NBIN(N)
         XHIS(N,I)=HMIN(N)+HDEL(N)*(FLOAT(I)-0.5)
         IHIS(N,I)=0
         HIST(N,I)=0
      ENDDO
      END

      SUBROUTINE MFILL(N,X,Y)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
C
      XI=((X-HMIN(N))/HDEL(N))
      I=INT(XI)+1
C
      IF(X.LT.HMIN(N)) I=0
      IF(X.GT.HMAX(N)) I=NBIN(N)+1
C
      IF(I.LT.1) THEN
         USCORE(N) = USCORE(N) + Y
         IUSCORE(N) = IUSCORE(N) + 1
      ELSEIF(I.GT.NBIN(N)) THEN
         OSCORE(N) = OSCORE(N) + Y
         IOSCORE(N) = IOSCORE(N) + 1
      ELSE
         IENT(N)=IENT(N)+1
         IHIS(N,I)=IHIS(N,I)+1
         HIST(N,I)=HIST(N,I)+Y
      ENDIF
      END


      SUBROUTINE MINTEG(NIN,NOUT,IDIR,IPOW)
C If IPOW=1 performs the integral of the distribution contained in histogram
C NIN up to the value specified by the abscissa (if IDIR=1) or from this
C value on (if IDIR=-1). The resulting integral distribution is put into
C NOUT, which is automatically booked if NOUT.ne.NIN .  Choosing IPOW=2
C the routine will return the square root of the integral of the squares,
C as is required, for example, for the propagation of the mean quadratic error
C of a given distribution. Overscores and underscores are included.
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
      CHARACTER*14  C
      DIMENSION C(2)
      DATA C/' INTEG BELOW X',' INTEG ABOVE X'/
      M = NBIN(NIN)
      I = (IDIR + 3)/2
      IF(NOUT.NE.NIN) THEN
        CALL MBOOK(NOUT,TITLE(NIN)//C(I),
     &                HDEL(NIN),HMIN(NIN),HMAX(NIN))
      ENDIF
      IF(IDIR.EQ.1) THEN
         HIST(NOUT,1) = SUMPOW(HIST(NIN,1),USCORE(NIN),IPOW)
         IHIS(NOUT,1) = IHIS(NIN,1) + IUSCORE(NIN)
         XHIS(NOUT,1) = XHIS(NIN,1) + HDEL(NIN)/2
         DO L=2,M
            HIST(NOUT,L) = SUMPOW(HIST(NIN,L),HIST(NOUT,L-1),IPOW)
            IHIS(NOUT,L) = IHIS(NIN,L) + IHIS(NOUT,L-1)
            XHIS(NOUT,L) = XHIS(NIN,L) + HDEL(NIN)/2
         ENDDO
         OSCORE(NOUT) = SUMPOW(OSCORE(NIN),HIST(NIN,M),IPOW)
         IOSCORE(NOUT) = IOSCORE(NIN) + IHIS(NIN,M)
      ELSEIF(IDIR.EQ.-1) THEN
         HIST(NOUT,M) = SUMPOW(HIST(NIN,M),OSCORE(NIN),IPOW)
         IHIS(NOUT,M) = IHIS(NIN,M) + IOSCORE(NIN)
         XHIS(NOUT,M) = XHIS(NIN,M) - HDEL(NIN)/2
         DO L=M-1,1,-1
            HIST(NOUT,L) = SUMPOW(HIST(NIN,L),HIST(NOUT,L+1),IPOW)
            IHIS(NOUT,L) = IHIS(NIN,L) + IHIS(NOUT,L+1)
            XHIS(NOUT,L) = XHIS(NIN,L) - HDEL(NIN)/2
         ENDDO
         USCORE(NOUT) = SUMPOW(USCORE(NIN),HIST(NIN,1),IPOW)
         IUSCORE(NOUT) = IUSCORE(NIN)+IHIS(NIN,1)
      ELSE
         CALL MBKWRN('MINTEG')
         WRITE(*,*) 'Wrong idir in minteg: OPERATION NOT PERFORMED'
         STOP
      ENDIF
      END

      FUNCTION SUMPOW(X,Y,IPOW)
      IF(IPOW.EQ.1) THEN
         SUMPOW = X + Y
      ELSEIF(IPOW.EQ.2) THEN
         SUMPOW = SQRT(X**2+Y**2)
      ELSEIF(IPOW.EQ.0) THEN
         CALL MBKWRN('SUMPOW')
         WRITE(*,*)'Error: IPOW=0 not allowed in SUMPOW'
         STOP
      ELSE
         SUMPOW = (X**IPOW+Y**IPOW)**(1./IPOW)
      ENDIF
      END

      SUBROUTINE MOPERA(I,OPER,J,K,X,Y)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
      CHARACTER OPER*1
      IF(NBIN(I).NE.NBIN(J).AND.(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.
     &    '*'.OR.OPER.EQ.'/'.OR.OPER.EQ.'M'.OR.OPER.EQ.'A')) THEN
          CALL MBKWRN('MOPERA')
          WRITE(20,*) I,J
  20      FORMAT(' ****** INCOMPATIBLE OPERATION HIST ',I2,' &',I2,
     &    '*******'/)
          STOP
      ENDIF
      IF(OPER.EQ.'E') THEN
c If I contains the accumulated weights, J the accumulated squares of the
c weights and IHIS(J,1) the number of accumulated entries, 'E' will add
c the average value of I to K and will put in J the quadratic dispersion.
         IF(IHIS(J,1).NE.0) THEN
            XXX = 1./IHIS(J,1)
         ELSE
            XXX = 0
         ENDIF
         DO L=1,NBIN(I)
            XSUM   = HIST(I,L)
            XSUMSQ = HIST(J,L)
            HIST(K,L)=HIST(K,L) + XXX*XSUM
            IHIS(K,L)=IHIS(K,L) + IHIS(I,L)
            HIST(J,L)=XXX*SQRT(ABS(XSUMSQ-XSUM**2*XXX))
         ENDDO
         IENT(K)=IENT(K)+IENT(I)
         XSUM = USCORE(I)
         XSUMSQ = USCORE(J)
         USCORE(K) = USCORE(K)+XXX*XSUM
         IUSCORE(K) = IUSCORE(K)+IUSCORE(I)
         USCORE(J) = XXX*SQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = OSCORE(I)
         XSUMSQ = OSCORE(J)
         OSCORE(K) = OSCORE(K)+XXX*XSUM
         IOSCORE(K) = IOSCORE(K)+IOSCORE(I)
         OSCORE(J) = XXX*SQRT(ABS(XSUMSQ-XSUM**2*XXX))
      ELSEIF(OPER.EQ.'Q') THEN
         DO L=1,NBIN(I)
            HIST(K,L) = SQRT(HIST(J,L)**2+HIST(I,L)**2)
         ENDDO
         USCORE(K) = SQRT(USCORE(J)**2+USCORE(I)**2)
         OSCORE(K) = SQRT(OSCORE(J)**2+OSCORE(I)**2)
      ELSEIF(OPER.EQ.'A') THEN
         DO L=1,NBIN(I)
            HIST(J,L) = HIST(J,L) + HIST(I,L)
            IHIS(J,L) = IHIS(J,L) + IHIS(I,L)
            HIST(K,L) = HIST(K,L) + HIST(I,L)**2
            IHIS(K,L) = IHIS(K,L) + 1
            IENT(K) = IENT(K)+1
            HIST(I,L) = 0
            IHIS(I,L) = 0
         ENDDO
         IENT(J) = IENT(J)+IENT(I)
         IUSCORE(J) = IUSCORE(J) + IUSCORE(I)
         USCORE(J) = USCORE(J) + USCORE(I)
         IOSCORE(J) = IOSCORE(J) + IOSCORE(I)
         OSCORE(J) = OSCORE(J) + OSCORE(I)
         IUSCORE(K) = IUSCORE(K) + 1
         USCORE(K) = USCORE(K) + USCORE(I)**2
         IOSCORE(K) = IOSCORE(K) + 1
         OSCORE(K) = OSCORE(K) + OSCORE(I)**2
         IENT(I) = 0
         IUSCORE(I) = 0
         IOSCORE(I) = 0
         USCORE(I) = 0
         OSCORE(I) = 0
      ELSE
        DO L=1,NBIN(I)
        IF(OPER.EQ.'+') THEN
                  HIST(K,L)=X*HIST(I,L) + Y*HIST(J,L)
        ELSEIF(OPER.EQ.'-') THEN
          HIST(K,L)=X*HIST(I,L) - Y*HIST(J,L)
        ELSEIF(OPER.EQ.'*') THEN
          HIST(K,L)=X*HIST(I,L) * Y*HIST(J,L)
        ELSEIF(OPER.EQ.'/') THEN
          IF(Y.EQ.0..OR.HIST(J,L).EQ.0.) THEN
            HIST(K,L)=0.
          ELSE
            HIST(K,L)=X*HIST(I,L) / (Y*HIST(J,L))
          ENDIF
                ELSEIF(OPER.EQ.'F') THEN
          HIST(K,L)=X*HIST(I,L)
        ELSEIF(OPER.EQ.'R') THEN
          IF(HIST(I,L).GT.0.) THEN
            HIST(K,L)=X*SQRT(HIST(I,L))
          ELSE
            HIST(K,L)=0.
          ENDIF
        ELSEIF(OPER.EQ.'S') THEN
          HIST(K,L)=X*HIST(I,L)**2
        ELSEIF(OPER.EQ.'L') THEN
          IF(HIST(I,L).EQ.0..OR.J.EQ.0.) THEN
             HIST(K,L)=0.
           ELSE
             HIST(K,L)=X*LOG10(Y*HIST(I,L))
           ENDIF
        ELSEIF(OPER.EQ.'M') THEN
           IF(I.NE.J) XNORM=HIST(I,L)
           IF(I.EQ.J) XNORM=FLOAT(IHIS(J,L))
           IF(XNORM.NE.0.) THEN
             XAVG=HIST(J,L)/XNORM
             HIST(K,L)=
     &       SQRT(ABS(-XAVG**2+HIST(K,L)/XNORM)/FLOAT(IHIS(I,L)))
             HIST(J,L)=XAVG
           ELSE
             HIST(K,L)=0.
             HIST(J,L)=0.
           ENDIF
        ELSEIF(OPER.EQ.'V') THEN
           XAVG=HIST(I,L)*X
           XSQAVG=HIST(J,L)*X
           XNORM=FLOAT(IHIS(I,L))*X
           IF(XNORM.NE.0.) THEN
              HIST(K,L)=SQRT(ABS(XSQAVG-XAVG**2)/XNORM)
              HIST(I,L)=XAVG
           ELSE
              HIST(K,L)=0.
           ENDIF
        ELSE
         CALL MBKWRN('MOPERA')
         WRITE(*,*) OPER
   5     FORMAT(' ****** OPERATION ="',A1,'" UNKNOWN ********'/)
         STOP
        ENDIF
        ENDDO
        IF(OPER.EQ.'+') THEN
                  USCORE(K)=X*USCORE(I) + Y*USCORE(J)
                  OSCORE(K)=X*OSCORE(I) + Y*OSCORE(J)
        ELSEIF(OPER.EQ.'-') THEN
          USCORE(K)=X*USCORE(I) - Y*USCORE(J)
          OSCORE(K)=X*OSCORE(I) - Y*OSCORE(J)
        ELSEIF(OPER.EQ.'*') THEN
          USCORE(K)=X*USCORE(I) * Y*USCORE(J)
          OSCORE(K)=X*OSCORE(I) * Y*OSCORE(J)
        ELSEIF(OPER.EQ.'/') THEN
          IF(Y.EQ.0..OR.USCORE(J).EQ.0.) THEN
            USCORE(K)=0.
          ELSE
            USCORE(K)=X*USCORE(I) / (Y*USCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.OSCORE(J).EQ.0.) THEN
            OSCORE(K)=0.
          ELSE
            OSCORE(K)=X*OSCORE(I) / (Y*OSCORE(J))
          ENDIF
                ELSEIF(OPER.EQ.'F') THEN
          USCORE(K)=X*USCORE(I)
          OSCORE(K)=X*OSCORE(I)
        ELSEIF(OPER.EQ.'R') THEN
          IF(USCORE(I).GT.0.) THEN
            USCORE(K)=X*SQRT(USCORE(I))
          ELSE
            USCORE(K)=0.
          ENDIF
          IF(OSCORE(I).GT.0.) THEN
            OSCORE(K)=X*SQRT(OSCORE(I))
          ELSE
            OSCORE(K)=0.
          ENDIF
        ELSEIF(OPER.EQ.'S') THEN
          USCORE(K)=X*USCORE(I)**2
          OSCORE(K)=X*OSCORE(I)**2
        ELSEIF(OPER.EQ.'L') THEN
          IF(USCORE(I).EQ.0..OR.J.EQ.0.) THEN
             USCORE(K)=0.
           ELSE
             USCORE(K)=X*LOG10(Y*USCORE(I))
           ENDIF
          IF(OSCORE(I).EQ.0..OR.J.EQ.0.) THEN
             OSCORE(K)=0.
           ELSE
             OSCORE(K)=X*LOG10(Y*OSCORE(I))
           ENDIF
        ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE MZERO(N)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
      BOOK(N)='RES'
      IENT(N)=0
      IUSCORE(N)=0
      IOSCORE(N)=0
      HAVG(N)=0.
      HINT(N)=0.
      DO 1 I=1,NBIN(N)
   1  HIST(N,I)=0.
      END

      SUBROUTINE MRESET(N)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
      BOOK(N)='RES'
      END

      SUBROUTINE PUTTAG(J,NAME)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
c Per marcare un istogramma
      CHARACTER * (*) NAME, TAG
      BOOK(J) = NAME
      RETURN
      ENTRY GETTAG(J,TAG)
      TAG = BOOK(J)
      END

      SUBROUTINE MFINAL(N)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3
      AVG=0
      XIN=0
      SIG=0
      IF=0
      DO J=1,NBIN(N)
         X=HIST(N,J)
         AVG=AVG+X*XHIS(N,J)
         XIN=XIN+X
         IF(X.NE.0) IF=1
      ENDDO
      IF(XIN.EQ.0) GO TO 10
      AVG = AVG/XIN
      DO J=1,NBIN(N)
         SIG=HIST(N,J)*(XHIS(N,J)-AVG)**2+SIG
      ENDDO
      SIG=SQRT(ABS(SIG/XIN))
 10   CONTINUE
      HINT(N) = XIN
      HAVG(N) = AVG
      HSIG(N) = SIG
      IF(IF.EQ.0) BOOK(N)='RES'
      END

      SUBROUTINE MNORM(N,X)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3,CTIME*8,CDATE*9
      IF(BOOK(N).NE.'YES')RETURN
      IF(HINT(N).EQ.0.) THEN
        CALL MBKWRN('MNORM')
        WRITE(*,*)' INTEGRAL HIST ',N,' IS ZERO: CANNOT RENORMALIZE'
        RETURN
      ELSE
        Y=X/HINT(N)
      ENDIF
      DO 1, I=1,NBIN(N)
    1 HIST(N,I)=HIST(N,I)*Y
      HINT(N)=X
      OSCORE(N)=OSCORE(N)*Y
      USCORE(N)=USCORE(N)*Y
      END

      SUBROUTINE MPRINT(N)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3,CTIME*8,CDATE*9
      DATA INI/0/
      IF(INI.EQ.0) THEN
c      CALL DATE(CDATE) RADEK
c      CALL TIME(CTIME) RADEK
      INI=1
      ENDIF
      IF(BOOK(N)(1:1).NE.'Y') RETURN
      WRITE(98,7) N,CDATE,CTIME(1:5)
      WRITE(98,*) TITLE(N)
      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      WRITE(98,'(3X,F10.4,2X,E15.4)')
     &                            XHIS(N,J),HIST(N,J)
    1 CONTINUE
      WRITE(98,15) HAVG(N),HSIG(N),HINT(N)
      WRITE(98,20) IENT(N),IUSCORE(N),IOSCORE(N)
    7 FORMAT(4X,'HIST = ',I3,1X,A9,1X,A5/)
   10 FORMAT(4X,2E10.3)
   15 FORMAT(/' AVG =',E10.3,4X,' RMS =',E10.3,' INTEGRAL =',E10.3,/)
   20 FORMAT(' ENTRIES=',I5,4X,'UNDERSCORE=',I5,4x,'OVERSCORE=',I5,//)
      END

      SUBROUTINE MTOP(N,M,BTIT,LTIT,SCALE)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3,CTIME*8,CDATE*9
      CHARACTER*(*) LTIT,BTIT,SCALE
      DATA INI/0/
      IF(INI.EQ.0) THEN
c      CALL DATE(CDATE) RADEK
c      CALL TIME(CTIME) RADEK
      INI=1
      ENDIF
      IF(BOOK(N)(1:1).NE.'Y') RETURN
      WRITE(99,100) TITLE(N),BTIT,LTIT,SCALE,HMIN(N),HMAX(N)
  100 FORMAT( /1x,
     &' SET WINDOW Y 2.5 TO 9.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET FONT DUPLEX '/1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A,'"',/1X,
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE LEFT ','"',A,'"',/1X,
     &' SET SCALE Y ',A,/1X,
     &' (SET TICKS TOP OFF)   '/1x,
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ')
      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      WRITE(99,'(3X,F10.4,2(2X,E15.4))')
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
    1 CONTINUE
      WRITE(99,200)
  200 FORMAT('   PLOT')
      WRITE(99,300) HINT(N),HAVG(N),HSIG(N),IENT(N),IUSCORE(N)
     &   ,IOSCORE(N),USCORE(N),OSCORE(N),CDATE,CTIME(1:5)
  300 FORMAT( /1x,
     &' BOX 6.25 0.9 SIZE 7.5 1.2'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &'(SET FONT SIMPLEX '/1X,
     &' TITLE 2.8 1.2 "INT =',1PE10.3,'   AVG =',1PE10.3,
     &             '   RMS =',1PE10.3,'"',/1X,
     &' TITLE 2.8 0.9 "Entries =',I8,2x,'Undersc =',I6,2X
     &                                 ,'Oversc =',I6,'"',/1X,
     &' TITLE 2.8 0.6 "ufloat=',1PE10.3,'ofloat=',
     &      1PE10.3,'"',/1X,
     &' TITLE 7.5 0.6 "',A9,2X,A5,'"',/1X,
     &' SET TITLE SIZE -2')
      WRITE(99,400)
  400 FORMAT('   NEW PLOT')
      END

      SUBROUTINE MNRTOP(N,M,BTIT,LTIT,SCALE)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      CHARACTER TITLE*50,BOOK*3,CTIME*8,CDATE*9
      CHARACTER*(*) LTIT,BTIT,SCALE
      DATA INI/0/
      IF(INI.EQ.0) THEN
c         CALL DATE(CDATE) RADEK
c         CALL TIME(CTIME) RADEK
         INI=1
      ENDIF
      CALL MFINAL(N)
      WRITE(99,100) TITLE(N)(1:20),BTIT,LTIT,SCALE,HMIN(N),HMAX(N)
  100 FORMAT( /1x,
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &'(SET FONT DUPLEX '/1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A,'"',/1X,
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE LEFT ','"',A,'"',/1X,
     &' SET SCALE Y ',A,/1X,
     &' (SET TICKS TOP OFF)   '/1x,
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ')
      DO 1 J=1,NBIN(N)
      WRITE(99,'(3X,F10.4,2(2X,E10.3))')
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
    1 CONTINUE
      WRITE(99,200)
  200 FORMAT('   PLOT')
      WRITE(99,300) HINT(N),HAVG(N),HSIG(N),IENT(N),IUSCORE(N)
     &   ,IOSCORE(N),USCORE(N),OSCORE(N),CDATE,CTIME(1:5)
  300 FORMAT( /1x,
     &' BOX 6.25 0.9 SIZE 7.5 1.2'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &'(SET FONT SIMPLEX '/1X,
     &' TITLE 2.8 1.2 "INT =',1PE10.3,'   AVG =',1PE10.3,
     &             '   RMS =',1PE10.3,'"',/1X,
     &' TITLE 2.8 0.9 "Entries =',I9,2x,'Usc =',I8,2X
     &                                 ,'Osc =',I8,'"',/1X,
     &' TITLE 2.8 0.6 "ufloat=',1PE10.3,' ofloat=',
     &      1PE10.3,'"',/1X,
     &' TITLE 7.5 0.6 "',A9,2X,A5,'"',/1X,
     &' SET TITLE SIZE -2')
      WRITE(99,400)
  400 FORMAT('   NEW PLOT')
      END


      SUBROUTINE MBOOK2(N,TIT,DELX,XMIN,XMAX,DELY,YMIN,YMAX)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3
      CHARACTER*(*) TIT
      NHIST2=MAX(N,NHIST2)
      IF(BOOK2(N)(1:1).EQ.'Y') THEN
	 CALL MBKWRN('MBOOK2')
         WRITE(*,*) 'Histogram',N,TITLE2(N),' already in use. '
         WRITE(*,*) 'superseded by ',TIT
      ENDIF
      TITLE2(N)='   '//TIT
      BOOK2(N)='YES'
C-- setup x-boundaries
      HDEL2(1,N)=DELX
      HMIN2(1,N)=XMIN
      HMAX2(1,N)=XMAX
      NBIN2(1,N)=INT((XMAX-XMIN)/DELX)
C-- setup y-boundaries
      HDEL2(2,N)=DELY
      HMIN2(2,N)=YMIN
      HMAX2(2,N)=YMAX
      NBIN2(2,N)=INT((YMAX-YMIN)/DELY)
      IENT2(N)=0
      IOSCORE2(N)=0
      HAVG2(1,N)=0.
      HAVG2(2,N)=0.
      HINT2(N)=0.
      DO 1 I=1,NBIN2(1,N)
      XHIS2(N,I)=HMIN2(1,N)+HDEL2(1,N)*(FLOAT(I)-0.5)
   1  CONTINUE
      DO 2 I=1,NBIN2(2,N)
      YHIS2(N,I)=HMIN2(2,N)+HDEL2(2,N)*(FLOAT(I)-0.5)
   2  CONTINUE
      DO 3 I=1,NBIN2(1,N)
      DO 3 J=1,NBIN2(2,N)
      HIST2(N,I,J)=0.
   3  CONTINUE
      END

      SUBROUTINE MFILL2(N,X,Y,WGT)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3
C-- find x-coordinate
      XI=((X-HMIN2(1,N))/HDEL2(1,N))+1.
      IF(XI.GT.100.OR.XI.LT.0.) THEN
      IOSCORE2(N)=IOSCORE2(N)+1
      RETURN
      ENDIF
      I=INT(XI)
C-- find y-coordinate
      YI=((Y-HMIN2(2,N))/HDEL2(2,N))+1.
      IF(YI.GT.100.OR.YI.LT.0.) THEN
      IOSCORE2(N)=IOSCORE2(N)+1
      RETURN
      ENDIF
      J=INT(YI)
      IF(I.GT.0.AND.I.LE.NBIN2(1,N).AND.J.GT.0.AND.J.LE.NBIN2(2,N))
     *                                             THEN
      IENT2(N)=IENT2(N)+1
      IHIS2(N,I,J)=IHIS2(N,I,J)+1
      HIST2(N,I,J)=HIST2(N,I,J)+WGT
      ELSE
      IOSCORE2(N)=IOSCORE2(N)+1
      ENDIF
      END


      SUBROUTINE MOPERA2(I,OPER,J,K,X,Y)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3
      CHARACTER OPER*1
      IF(NBIN2(1,I).NE.NBIN2(1,J).OR.NBIN2(2,I).NE.NBIN2(2,J).
     &AND.(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.
     &'*'.OR.OPER.EQ.'/'.OR.OPER.EQ.'M')) GO TO 10
      DO L1=1,NBIN2(1,I)
      DO L2=1,NBIN2(2,I)
      IF(OPER.EQ.'+') THEN
      HIST2(K,L1,L2)=X*HIST2(I,L1,L2) + Y*HIST2(J,L1,L2)
      ELSEIF(OPER.EQ.'-') THEN
      HIST2(K,L1,L2)=X*HIST2(I,L1,L2) - Y*HIST2(J,L1,L2)
      ELSEIF(OPER.EQ.'*') THEN
      HIST2(K,L1,L2)=X*HIST2(I,L1,L2) * Y*HIST2(J,L1,L2)
      ELSEIF(OPER.EQ.'/') THEN
        IF(Y.EQ.0..OR.HIST2(J,L1,L2).EQ.0.) THEN
          HIST2(K,L1,L2)=0.
          ELSE
          HIST2(K,L1,L2)=X*HIST2(I,L1,L2) / (Y*HIST2(J,L1,L2))
        ENDIF
      ELSEIF(OPER.EQ.'F') THEN
      HIST2(K,L1,L2)=X*HIST2(I,L1,L2)
      ELSEIF(OPER.EQ.'R') THEN
        IF(HIST2(I,L1,L2).GT.0.) THEN
        HIST2(K,L1,L2)=X*SQRT(HIST2(I,L1,L2))
        ELSE
        HIST2(K,L1,L2)=0.
        ENDIF
      ELSEIF(OPER.EQ.'S') THEN
      HIST2(K,L1,L2)=X*HIST2(I,L1,L2)**2
      ELSEIF(OPER.EQ.'L') THEN
        IF(HIST2(I,L1,L2).EQ.0..OR.J.EQ.0.) THEN
             HIST2(K,L1,L2)=0.
             ELSE
             HIST2(K,L1,L2)=X*LOG10(Y*HIST2(I,L1,L2))
        ENDIF
      ELSE
      WRITE(98,5) OPER
   5  FORMAT(' ****** OPERATION ="',A1,'" UNKNOWN ********'/)
      RETURN
      ENDIF
      END DO
      ENDDO
      RETURN
  10  WRITE(98,20) I,J
  20  FORMAT(' ****** INCOMPATIBLE OPERATION HIST2 ',I2,' &',I2,
     &                                                   '*******'/)
      END

      SUBROUTINE MFINAL2(N)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3
      IF(BOOK2(N).NE.'YES') RETURN
      XIN=0.
C-- projection on the x-axis
      DO 2 I=1,NBIN2(1,N)
        DO 1 J=1,NBIN2(2,N)
   1    XPROJ(N,I)=XPROJ(N,I)+HIST2(N,I,J)
   2  XIN=XIN+XPROJ(N,I)
      IF(XIN.EQ.0.) GO TO 10
C-- projection on the y-axis
      DO 3 J=1,NBIN2(2,N)
        DO 3 I=1,NBIN2(1,N)
   3    YPROJ(N,J)=YPROJ(N,J)+HIST2(N,I,J)
      HINT2(N)=XIN
      RETURN
  10  BOOK2(N)=' NO'
      END

      SUBROUTINE PROHIS(N,M,IAX)
C-- projects the scatter plot N onto the IAX axis (x,y->IAX=1,2) and
C   put the contents in histogram M, after automatically booking it.
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE*50,BOOK*3
      CHARACTER TITLE2*50,BOOK2*3
      BOOK(M)='YES'
      NBIN(M)=NBIN2(IAX,N)
      HDEL(M)=HDEL2(IAX,N)
      HMIN(M)=HMIN2(IAX,N)
      HMAX(M)=HMAX2(IAX,N)
      NHIST=MAX(NHIST,M)
      TITLE(M)=TITLE2(N)//'(PROJ)'
      DO I=1,NBIN(M)
      IF(IAX.EQ.1)      THEN
      HIST(M,I)=XPROJ(N,I)
      XHIS(M,I)=XHIS2(N,I)
      ELSEIF(IAX.EQ.2)      THEN
      HIST(M,I)=YPROJ(N,I)
      XHIS(M,I)=YHIS2(N,I)
      ENDIF
      ENDDO
      END


      SUBROUTINE MNORM2(N,X)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3
      IF(BOOK2(N).NE.'YES')RETURN
      DO 1, I=1,NBIN2(1,N)
      DO 1, J=1,NBIN2(2,N)
    1 HIST2(N,I,J)=HIST2(N,I,J)/HINT2(N)*X
      HINT2(N)=X
      END

      SUBROUTINE MPRINT2(N)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3,CTIME*8,CDATE*9
      DATA INI/0/
      IF(INI.EQ.0) THEN
c      CALL DATE(CDATE) RADEK
c      CALL TIME(CTIME) RADEK
      INI=1
      ENDIF
      IF(BOOK2(N).NE.'YES') RETURN
      WRITE(98,7) N,CDATE,CTIME(1:5)
      WRITE(98,*) TITLE2(N)
      WRITE(98,10) ((XHIS2(N,I),YHIS2(N,J),HIST2(N,I,J),J=1,NBIN2(2,N))
     *                            ,I=1,NBIN2(1,N))
      WRITE(98,20) IENT2(N),HINT2(N),IOSCORE2(N)
    7 FORMAT(4X,'HIST = ',I3,1X,A9,1X,A5/)
   10 FORMAT(4X,3E10.3)
   20 FORMAT(' ENTRIES=',I5,4X,' INTEGRAL =',E10.3,4x,
     *              'OVERSCORE=',I5,//)
      END


      SUBROUTINE MTOP2(N,BTIT,LTIT,PLOPT)
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE2*50,BOOK2*3,CTIME*8,CDATE*9
      CHARACTER*(*) LTIT,BTIT,PLOPT
      DOUBLE PRECISION RANXX
      DATA INI/0/                  
      IF(INI.EQ.0) THEN
      ISEED=2143567+2*INT(SECNDS(0.0))
c      CALL DATE(CDATE) RADEK
c      CALL TIME(CTIME) RADEK
      INI=1
      ENDIF
      IF(BOOK2(N)(:1).NE.'Y') RETURN
      WRITE(99,100) TITLE2(N),BTIT,LTIT,HMIN2(1,N),HMAX2(1,N),
     *   HMIN2(2,N),HMAX2(2,N)
  100 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 9.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET FONT DUPLEX '/1X, 
     &' TITLE TOP ','"',A,'"',/1X,
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE LEFT ','"',A,'"',/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET LIMITS Y ',F10.5,' ',F10.5,/1X,
     &' PLOT AXES',/1X,
     &' SET ORDER X Y ')

C-- squares with area proportional to the bin weight
      IF(PLOPT.EQ.'BOX') THEN
      DO 1, I=1,NBIN2(1,N)
      DO 1, J=1,NBIN2(2,N)
      SCMAX=MAX(SCMAX,HIST2(N,I,J))
  1   CONTINUE
      SCMAX=0.8*SCMAX
      DO 2, I=1,NBIN2(1,N)
      DO 2, J=1,NBIN2(2,N)
      IF(HIST2(N,I,J).EQ.0.) GO TO 2      
      WRITE(99,200) XHIS2(N,I),YHIS2(N,J),HIST2(N,I,J)/SCMAX*0.3
  2   CONTINUE
  200 FORMAT(2X,'BOX',F14.8,2X,F14.8,3X,'DATA SIZE',F14.8)
      WRITE(99,300)
  300 FORMAT('   PLOT')
      WRITE(99,400) HINT2(N),IENT2(N)
     &   ,IOSCORE2(N),CDATE,CTIME(1:5)
  400 FORMAT( /1x,                               
     &' BOX 6.25 1.0 SIZE 7.5 1.'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' SET FONT DUPLEX '/1X,
     &' TITLE 2.8 1.2 "INT =',1PE10.3,'    Entries =',I8,2x
     &                                 ,'Oversc =',I6,'"',/1X,
     &' TITLE 2.8 0.8 "',A9,2X,A5,'"'/1X,
     &' SET TITLE SIZE -2')
      WRITE(99,500)
  500 FORMAT('   NEW PLOT')

C-- dots - as many as the value of HINT2(n), distributed according to wgt
      ELSEIF(PLOPT.EQ.'DOT') THEN                                        
      XINT=0.
      DO 11, I=1,NBIN2(1,N)
      DO 11, J=1,NBIN2(2,N)
      XINT=XINT+HIST2(N,I,J)    
 11   CONTINUE
      WRITE(99,*) '  SET INTENSITY 3'
      DO 3, I=1,NBIN2(1,N)          
      DO 3, J=1,NBIN2(2,N)
      IF(HIST2(N,I,J).EQ.0..OR.XINT.EQ.0.) GO TO 3
      NDOTS=INT(500*HIST2(N,I,J)/XINT)
C      RN=SNGL(RANXX(ISEED))
C      REST=HIST2(N,I,J)-FLOAT(NDOTS)
C      IF(RN.LT.REST) NDOTS=NDOTS+1
      DO ND=1,NDOTS
      X=XHIS2(N,I)+(2.*SNGL(RANXX(ISEED))-1.)*HDEL2(1,N)
      Y=YHIS2(N,J)+(2.*SNGL(RANXX(ISEED))-1.)*HDEL2(2,N)
      WRITE(99,*) X,Y
      ENDDO
      NTOT=NTOT+NDOTS
  3   CONTINUE       
      WRITE(99,*) '  SET INTENSITY 1'
      WRITE(99,300)                 
      WRITE(99,600) HINT2(N),IENT2(N),NTOT
     &   ,CDATE,CTIME(1:5)
  600 FORMAT( /1x,                               
     &' BOX 6.25 1.0 SIZE 7.5 1.'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' SET FONT DUPLEX '/1X,
     &' TITLE 2.8 1.2 "INT =',1PE10.3,'    Entries =',I8,2x
     &                                 ,'Points =',I6,'"',/1X,
     &' TITLE 2.8 0.8 "',A9,2X,A5,'"'/1X,
     &' SET TITLE SIZE -2')
      WRITE(99,500)

      ENDIF
      END


      FUNCTION RANXX(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION MINV,RANXX
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      RANXX = SEED*MINV
      END



      SUBROUTINE MULTITOP(NH,NE,N,M,BTIT,LTIT,SCA)
      PARAMETER (NMB=200)
      COMMON/HISTO/HIST(NMB,100),XHIS(NMB,100),HDEL(NMB),HMIN(NMB)
     &,HMAX(NMB),USCORE(NMB),OSCORE(NMB)
     &,NBIN(NMB),IHIS(NMB,100),IUSCORE(NMB),IOSCORE(NMB)
     &,IENT(NMB),HAVG(NMB),HINT(NMB),HSIG(NMB),BOOK(NMB),TITLE(NMB)
     &,NHIST
c MODIFICA
      common/mttini/ini
c                   
      REAL LAB0,LABS
      CHARACTER TITLE*50,BOOK*3,CTIME*8,CDATE*9,SCALE*3
      CHARACTER*(*) LTIT,BTIT,SCA
      CHARACTER*7  PLOT(4)
      DATA PLOT/'SOLID','DASHES','DOTS','DOTDASH'/
C  PLOT SIZE, CORNERS
      DATA WIDTH,HEIGHT/11.5,8.5/,XCORN,YCORN/1.5,1./
C  PLOT VERSUS TEXT FRACTION
      DATA XPFRAC,YPFRAC/0.75,0.75/,XTFRAC,YTFRAC/0.25,0.25/
C  DEFAULT SIZES
      DATA TIT0,LAB0,TIC0/-3.,-3.,0.06/
c MODIFICA      DATA INI/0/
      IF(INI.EQ.0) THEN
C      CALL DATE(CDATE) RADEK
C      CALL TIME(CTIME) RADEK
      IFRAME=0
      WRITE(99,71) CDATE,CTIME(1:5)
   71 FORMAT(4X,' ( ',A9,1X,A5/)
      INI=1
      ENDIF
      IF(SCA.EQ.'REF') THEN
        IFRAME=0 
        RETURN
      ENDIF
      IF(BOOK(NH)(:1).NE.'Y') RETURN
      IFRMAX=N*M
      IFRAME=IFRAME+1
      IF(IFRAME.GT.IFRMAX.OR.N.NE.NOLD.OR.M.NE.MOLD) THEN
        IFRAME=1
        WRITE(99,202)
        WRITE(99,1) CDATE,CTIME(1:5)
  1     FORMAT(' ( SET FONT DUPLEX',/,'  SET TITLE SIZE 2',/,
     +      ' TITLE 12.8 9 ANGLE -90 ','" MLM   ',A9,1X,A5,'"')
      ENDIF
      IF(IFRAME.EQ.1) THEN
        I=1
        J=1
      ELSEIF(IFRAME.LE.IFRMAX) THEN
        IF(I.LE.N) I=I+1
        IF(I.GT.N) THEN
                I=1
                J=J+1
        ENDIF
      ENDIF
      IF(N.EQ.NOLD) GO TO 10
      NS=N-1
      XD=WIDTH/FLOAT(N)
      SRED=SQRT(FLOAT(N*M))
      TITS=TIT0/SRED
      LABS=LAB0/SRED
      TICS=TIC0/SRED
      XTIT0=0.55*XPFRAC*XD
      NOLD=N
10    IF(M.EQ.MOLD) GO TO 20
      YD=HEIGHT/FLOAT(M)
      YTIT0=0.06*YD
      MOLD=M
20    CONTINUE
      XL=(I-1)*XD + XCORN
      YL=(M-J)*YD + YCORN
      XU=XL+XD*XPFRAC
      YU=YL+YD*YPFRAC
      IP=0
c inizio modifiche
      XMX=0.
      DO IBIN=1,NBIN(NH)
        X=HIST(NH,IBIN)
        IF(X.NE.0.) THEN
           IF(XMX.EQ.0.) THEN
              FMX = X + HIST(NE,IBIN)
              FMN = X - HIST(NE,IBIN)
           ELSE
              FMX=MAX(FMX,X + HIST(NE,IBIN))
              FMN=MIN(FMN,X - HIST(NE,IBIN))
           ENDIF
        ENDIF
        XMX=MAX(XMX,ABS(X)+ HIST(NE,IBIN))
      ENDDO
      IF(XMX.EQ.0.) GO TO 203
      SCALE=SCA
50    IF(SCALE.EQ.'LIN') THEN
        IF(FMN.GE.0.)   FMIN=0.
        IF(FMN.LT.0.)   FMIN=FMN*1.3
        IF(FMX.GT.0.)   FMAX=FMX*1.3
        IF(FMX.LT.0.)   FMAX=0.
      ELSEIF(SCALE.EQ.'LOG') THEN
        IF(FMN.LE.0.) THEN
                SCALE='LIN'
                GO TO 50
        ENDIF
        FMAX=10.**( AINT(LOG10(FMX)+1000001) - 1000000 )
        FMIN=10.**( AINT(LOG10(FMN)+1000000) - 1000000 )
      ENDIF
C fine modifiche
      WRITE(99,100) TITS,LABS,TICS,XL,XU,YL,YU
100   FORMAT(2X,'( SET FONT DUPLEX',/,
     *       2X,'SET TITLE SIZE ',F8.4,/,
     *       2X,'SET LABEL SIZE ',F8.4,/,
     *       2X,'SET TICKS TOP OFF SIZE ',F8.4,/,
     *       2X,'SET WINDOW X ',F8.4,' TO ',F8.4,/,
     *       2X,'SET WINDOW Y ',F8.4,' TO ',F8.4)
      XTIT=XL+XTIT0
      YTIT=YU+YTIT0
      WRITE(99,101) XL,YTIT,TITLE(NH)(1:40)
101   FORMAT('  TITLE ',2(F8.4,1X),'"',A,'"')
      YTIT=YTIT-2.*YTIT0
      WRITE(99,102) XTIT,YTIT,HINT(NH)
102   FORMAT('  TITLE ',2(F8.4,1X),'" INT=',1PE10.3,'"')
      YTIT=YTIT-YTIT0
      WRITE(99,103) XTIT,YTIT,IENT(NH)
103   FORMAT('  TITLE ',2(F8.4,1X),'" ENT=',I9,'"')
      YTIT=YTIT-YTIT0
      IF(USCORE(NH).NE.0.) THEN
        WRITE(99,104) XTIT,YTIT,USCORE(NH)
104     FORMAT('  TITLE ',2(F8.4,1X),'" UFL=',1PE10.3,'"')
        YTIT=YTIT-YTIT0
      ENDIF
      IF(OSCORE(NH).NE.0.) THEN
        WRITE(99,105) XTIT,YTIT,OSCORE(NH)
105     FORMAT('  TITLE ',2(F8.4,1X),'" OFL=',1PE10.3,'"')
        YTIT=YTIT-YTIT0
      ENDIF
      WRITE(99,106) XTIT,YTIT,XU,YTIT,XTIT,YTIT,XTIT,YU
106   FORMAT(2X,'SET ORD X Y ',/,2(F8.4,1X),/,2(F8.4,1X),/,
     *       2X,'JOIN TEXT',/,
     *       2X,2(F8.4,1X),/,2(F8.4,1X),/,
     *       2X,'JOIN TEXT')
      WRITE(99,108) TITS*1.5
108   FORMAT(2X,'SET TITLE SIZE ',F8.4)
      WRITE(99,107) BTIT,XL-0.75*XD*XTFRAC,YL+(YU-YL)/3.,LTIT,SCALE,
     * HMIN(NH),HMAX(NH),FMIN,FMAX
107   FORMAT(
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE ',f10.5,f10.5,' ANGLE 90 ','"',A,'"',/1X,
     &' SET SCALE Y ',A,/1X,
     &' SET TICKS TOP OFF   '/1x,
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET LIMITS Y ',1PE10.3,' ',1PE10.3,/1X,
     &' SET ORDER X Y DY')
C
C  END HEADER , FILL TOPDRAWER WITH DATA
C
      ENTRY MTFILL(NH,NE,N,M,BTIT,LTIT,SCA)
      IP=IP+1
      IF(IP.GT.4) IP=1
      WRITE(99,110) TITLE(NH),HINT(NH),IENT(NH)
110   FORMAT(' ( ',A,/,' ( INT=',1PE10.3,'  ENTRIES=',I12)
      DO 200 IBIN=1,NBIN(NH)
      IF(HIST(NH,IBIN).EQ.0..AND.HIST(NE,IBIN).EQ.0.) GO TO 200
      WRITE(99,'(3X,F10.4,2(2X,E15.4))')
     &          XHIS(NH,IBIN),HIST(NH,IBIN),HIST(NE,IBIN)
200   CONTINUE
      WRITE(99,201)  PLOT(IP)
      IF(BOOK(NE).NE.'NO')   WRITE(99,*)  '  PLOT'
201   FORMAT(2X,'HIST ',A)
202   FORMAT('   NEW PLOT',/,/)
203   RETURN
      END

      SUBROUTINE NEWPLOT
      WRITE(99,202)
202   FORMAT('   NEW PLOT',/,/)
      CALL MULTITOP(IDUM,IDUM,IDUM,IDUM,' ',' ','REF')
      END                                    

      subroutine MBKWRN(str)
      character *(*) str
      write(*,*) '********** WARNING **********'
      write(*,*) '*********  ',str,'  *********'
      end
C*******************************************************************
C     END OF THE HISTOGRAMMING PACKAGE
C*******************************************************************
