c Interface to pdflib routines, version 4.00 and beyond (PDF set
c is identified by three parameters, NPTYPE, NGROUP, NSET).
c Includes also electron distribution functions of the jet package
c Modified on 8/11/01, since pdflib version 8.04 has in one case
c 100 pdf sets for one author group
c
      subroutine mlmpdf(mode,ih,q2,x,fx,nf)
      implicit real * 8 (a-h,o-z)
      real * 4 q2,x,fx(-nf:nf)
      common/trans/nptype,ngroup,nset
c set print flag. Set iflprt=2 in order to print the content
c of the common blocks w50511, w50512, and w50513
      common/w50510/iflprt
c x*pdf(x) in PDG format
      dimension fxp(-6:6)
c used by pdfset
      dimension val(20)
c used by pdfset
      character * 20 parm(20)
      logical ini
      data ini/.true./,imode/0/
c
      if(ih.ne.5)then
c incoming particle is not an electron: use PDFLIB
        if(ini.or.imode.ne.mode) then
          ini = .false.
          imode = mode
c pass from our conventions to PDFLIB conventions for parameter settings. 
c See subroutine newmode for comments on the conventions adopted
          call newmode(ih,mode)
          parm(1) = 'NPTYPE'
          val (1) =  nptype
          parm(2) = 'NGROUP'
          val (2) =  ngroup
          parm(3) = 'NSET'
          val (3) =  nset
          call pdfset(parm,val)
        endif
        xd = dble(x)
        qd = dble(sqrt(q2))
        call pftopdg(xd,qd,fxp)
c in our conventions 1=up, 2=down; PDFLIB 1=down, 2=up. With
c f(1)<-->f(2) we mean also f(-1)<-->f(-2)
c in the following lines, deals with particles only (no antiparticles)
c proton(ih=1) ==> f(1)<-->f(2)
c neutron(ih=2) ==> no action (f(1)<-->f(2) for PDFLIB convention and
c    f(1)<-->f(2) for isospin symmetry (u_proton=d_neutron....)
c pion+(ih=3) ==> f(2)<-->f(-2), since PDFLIB has d=u=q_v+q_sea, 
c    ubar=dbar=q_sea
c photon(ih=4) ==> f(-1)<-->f(-2) and f(i)=f(-i)/2, i=1,2 since PDFLIB 
c    has f(i)=2*f(-i), and f(1)<-->f(2)
c Notice that in the jet package pions and neutrons are not used. If
c selected, they are rejected by the routine pdfpar. This routine
c is therefore a completely general interface with PDFLIB
        if(abs(ih).eq.1) then
          tmp     = fxp(1)
          fxp(1)  = fxp(2)
          fxp(2)  = tmp
          tmp     = fxp(-1)
          fxp(-1) = fxp(-2)
          fxp(-2) = tmp
        elseif(abs(ih).eq.2) then
          continue
        elseif(abs(ih).eq.3) then
          tmp = fxp(-2)
          fxp(-2) = fxp(2)
          fxp(2)  = tmp
        elseif(abs(ih).eq.4) then
          tmp     = fxp(-1)
          fxp(-1) = fxp(-2)
          fxp(-2) = tmp
          fxp(1)=fxp(-1)
          fxp(2)=fxp(-2)
c this is (p+n)/2
        elseif(ih.eq.0) then
          va  = (fxp(1)+fxp(2))/2
          sea = (fxp(-1)+fxp(-2))/2
          fxp(1)  = va
          fxp(2)  = va
          fxp(-1) = sea
          fxp(-2) = sea
        else
          write(*,*) ' ih was', ih,' instead of 0, +-1, +-2, +-3, or 4'
          stop
        endif
c for particles, ich=1, for antiparticles, ich=-1
        if(ih.lt.0) then
          ich = -1
        else
          ich = 1
        endif
c divide by x and exchange q with qbar in the case of antiparticles
        do j=-nf,nf
          fx(j) = sngl(fxp(j*ich)/xd)
        enddo
c correct for a bug in PDFLIB specific to AFG densities in the photon
        if(ih.eq.4.and.mode.eq.603)then
          do j=-nf,nf
            if(j.ne.0)fx(j)=fx(j)/2.d0
          enddo
        endif
      else
c incoming "hadron" is an electron
        if(ini.or.imode.ne.mode) then
          ini = .false.
          imode = mode
          iset=mode-50
        endif
        if(iset.eq.1) then
          call elpdf_lac1(q2,x,fx,nf)
        elseif(iset.eq.2) then
          call elpdf_grv(q2,x,fx,nf)
        elseif(iset.eq.3) then
          call elpdf_user(q2,x,fx,nf)
        else
          write(*,*)'Electron set non implemented'
          stop
        endif
      endif
      return
      end


      subroutine pdfpar(mode,ih,xlam,scheme,iret)
      implicit real * 8 (a-h,o-z)
      common/trans/nptype,ngroup,nset
      common/w50512/qcdl4,qcdl5
      character*10 pdfver(3)
      common/w505190/pdfver
      dimension val(20)
      character * 20 parm(20)
      character * 2 scheme
      logical ini
      data ini/.true./
      logical nopions
      data nopions/.false./
c iret#0 when problem occur
      iret = 0
      if(ini)then
        write(*,*)'This is an interface to PDFLIB'
        write(*,*)'PDFLIB version number:',pdfver(1)
        write(*,*)'Released:',pdfver(2)
        write(*,*)'On:',pdfver(3)
        ini = .false.
      endif
      if(nopions)then
        if(abs(ih).ne.1.and.ih.ne.4.and.ih.ne.5)then
          write(*,*) ' hadron tpye ',ih,' not implemented'
          iret=1
          return
        endif
      endif
c fake values. If kept, the main program crashes
      scheme='XX'
      xlam=0.0
c incoming particle is not an electron: use PDFLIB
      if(ih.ne.5)then
c the scheme of the PDFLIB set is not given in any common block.
c Set it by hand in the main program
        scheme = '**'
        call newmode(ih,mode)
        parm(1) = 'NPTYPE'
        val (1) =  nptype
        parm(2) = 'NGROUP'
        val (2) =  ngroup
        parm(3) = 'NSET'
        val (3) =  nset
c print the relevant parameters for the PDF set chosen
c        iflprt = 2
c set the parameters
        call pdfset(parm,val)
c Lambda_QCD_5, as given by PDFLIB
        xlam = qcdl5
      else
c incoming particle is an electron
        if(mode.eq.51)then
          scheme='MS'
          xlam=.130
        elseif(mode.eq.52)then
          scheme='DG'
          xlam=.130
        elseif(mode.eq.53)then
          scheme='**'
          xlam=.001
        else
          write(*,*)'Electron set not implemented'
          iret=1
        endif
      endif
      return
      end


      subroutine prntsf
c     prints details of the structure function sets
c
      write(*,100)                             
     #  '  For nucleons, pions and photons, enter the set number'
     # ,'  '
     # ,'              Set # = 1000*NGROUP+NSET'
     # ,'  '
     # ,'  where NGROUP and NSET are the parameters of the PDFLIB'
     # ,'  version you are using (please read the user manual of'
     # ,'  your version of PDFLIB to check the listing)'
     # ,'  '
      write(*,100)                             
     #  '  For electrons, use'
     # ,'  '
     # ,'Set #        Set'
     # ,'  '
     # ,'  51         LAC1'
     # ,'  52         GRV-HO'
     # ,'  53         user defined'
 100  format(1x,a,100(/,1x,a))
      return
      end


      subroutine newmode(ih,mode)
c This subroutine converts our conventions for the identification
c of a set of PDFs into PDFLIB conventions. We use
c
c                     JET PACKAGE             PDFLIB
c 
c Particle type 
c
c  nucleons           -2,-1,0,1,2                1
c  pions                  -3,3                   2
c  photons                  4                    3
c 
c Given the particle type, our package uses a single number (MODE)
c to identify to PDF set, while PDFLIB uses two numbers (NGROUP and NSET).
c The value of MODE which corresponds to the values (NGROUP,NSET)
c is by definition given by
c
c                        MODE=1000*NGROUP+NSET
c
c This is working as long as every group has less the 1000 PDFs sets....
      implicit real * 8 (a-h,o-z)
      common/trans/nptype,ngroup,nset
c
      if(abs(ih).le.2)then
        nptype=1
      elseif(abs(ih).eq.3)then
        nptype=2
      elseif(ih.eq.4)then
        nptype=3
      else
        write(6,*)'hadron type not implemented in PDFLIB'
        stop
      endif
      ngroup=int(dfloat(mode)/1000.e0)
      nset=mode-1000*ngroup
      return
      end


C------- ALPHA QCD -------------------------------------
c Program to calculate alfa strong with nf flavours,
c as a function of lambda with 5 flavors.
c The value of alfa is matched at the thresholds q = mq.
c When invoked with nf < 0 it chooses nf as the number of
c flavors with mass less then q.
c
      function alfas(q2,xlam,inf)
      implicit real * 8 (a-h,o-z)
      data olam/0.d0/,pi/3.14159d0/
      data xmb/5.d0/,xmc/1.5d0/
      if(xlam.ne.olam) then
        olam = xlam
        b5  = (33-2*5)/pi/12
        bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
        b4  = (33-2*4)/pi/12
        bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
        b3  = (33-2*3)/pi/12
        bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
        xlc = 2 * log(xmc/xlam)
        xlb = 2 * log(xmb/xlam)
        xllc = log(xlc)
        xllb = log(xlb)
        c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     #        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
        c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     #        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
      endif
      q   = sqrt(q2)
      xlq = 2 * log( q/xlam )
      xllq = log( xlq )
      nf = inf
      if( nf .lt. 0) then
        if( q .gt. xmb ) then
          nf = 5
        elseif( q .gt. xmc ) then
          nf = 4
        else
          nf = 3
        endif
      endif
      if    ( nf .eq. 5 ) then
        alfas = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
        alfas = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
      elseif( nf .eq. 3 ) then
        alfas = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
      else
        print *,'error in alfa: unimplemented # of light flavours',nf
        stop
      endif
      return
      end
