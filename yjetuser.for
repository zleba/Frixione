c-----------------------------------------------------------------------
c **************************** USER ROUTINES ***************************
c-----------------------------------------------------------------------
c
c Weizsacker-Williams function for photon in electron. See also
c S. Frixione, M. Mangano, P. Nason, G. Ridolfi, Phys. Lett. B319(1993)339.
c
c THIS FUNCTION IS IRRELEVANT IN PROCESSES OTHER THAN PHOTON-HADRON COLLISIONS
c
      function phdistr(x,xmu2)
      implicit real * 8 (a-z)
      parameter (alfaem=1/137.d0)
      parameter (pi=3.14159265358979312D0)
      parameter (aemo2pi=alfaem/(2*pi))
      parameter (xme = 0.511d-3)
      parameter (xme2 = xme**2)
c cuts on E_gamma/E_e. x_inf and x_sup can be freely changed,
c provided that   0 <= x_inf < x_sup <= 1  (<= means less than or equal to)
      parameter (x_inf=0.2)
      parameter (x_sup=0.8)
      common/shadr/sh
c
c q2eff is the effective Weizsaecker-Williams scale SQUARED (in GeV^2). 
c It can be set in this function or using ZGMU2. In the latter case, 
c uncomment the following line and comment the next one
c      q2eff = xmu2
      q2eff = 4.d0
      tmp=0.d0
      if(x.gt.x_inf.and.x.lt.x_sup)
     #  tmp = (1+(1-x)**2)/x * log( q2eff*(1-x)/(xme*x)**2 )
     #        +2*xme2*x*( 1/q2eff-(1-x)/(xme*x)**2 )
      phdistr=aemo2pi*tmp
      end


c****************************************************************
c Begin of the user analysis routines
c****************************************************************


c----------------------------------------------
      subroutine inijet
c initialize histograms
c----------------------------------------------
      implicit real*4(a-z)
      real*8 sh
c mxbin is the maximum # of bins per histogram: do not change
      parameter (pi=3.14160E0, mxbin=100)
c sh is the square of the total CM energy of the colliding system
      common/shadr/sh
c
      ptmax=70.e0
      ptmin=20.e0
      ptbin=2.e0
      el=5.e0
      etabin=.2e0
      phibin=pi/40.e0
c Reset histograms
      call inihist
c Here book the needed histograms, using:
c bookup(histograms number, title, binsize, xmin, xmax) (single precision!)
c !!!!!!!!! HISTOGRAM NUMBER= N*4+1, N=0,1,...49. !!!!!!!!!!!!!
c The other histograms are used internally for statistics. For example
c bookup(1,...) uses the histogram 1,2,3,4. 
c
      call bookup(1,'pt, cone [R=1]',ptbin,ptmin,ptmax)
      call bookup(5,'pt, cone [R=0.7]',ptbin,ptmin,ptmax)
      call bookup(9,'pt, ES [D=1]',ptbin,ptmin,ptmax)
c
      call bookup(13,'eta, cone [R=1]',etabin,-el,el)
      call bookup(17,'eta, cone [R=0.7]',etabin,-el,el)
      call bookup(21,'eta, ES [D=1]',etabin,-el,el)
c
      call bookup(25,'dphi, cone [R=1]',phibin,0.e0,pi)
      call bookup(29,'dphi, cone [R=0.7]',phibin,0.e0,pi)
      call bookup(33,'dphi, ES [D=1]',phibin,0.e0,pi)
c
      return
      end


c-------------------------------------
      subroutine topout
c output histograms as topdrawer file
c-------------------------------------
c Do NOT remove the following lines
      do i=1,200
        call mfinal(i)
      enddo
c Multitop creates a topdrawer file with nx*ny small plots on the
c same page. If there are more than nx*ny plots, it goes automatically
c to the next page.
c Usage: multitop(hn, hn+3, nx, ny, x label, y label, y scale)
c where:
c hn is the histogram number (hn=4*N+1, N=0,...,49)
c nx is the number of plots per row
c ny is the number of plots per column.
c x label is the label that will appear under the x axis
c y label is the label that will appear next to the y axis
c y scale='LIN' for linear scale and 'LOG' for logarithmic scale on y axis
c mtfill has the same sintax of multitop, but it acts differently:
c it adds the desired histogram to the plot defined by the last call 
c to multitop. This allows to have two (or more if further calls to mtfill are
c made) histos in the same plot. The program will automatically change the
c pattern of the new added histogram (in the sequence solid,dashes,dots,
c dotdashes, ...) 
c
      call multitop(1,4,2,2,'pt','mub/bin','LOG')
      call mtfill(5,8,2,2,'pt','mub/bin','LOG')
      call mtfill(9,12,2,2,'pt','mub/bin','LOG')
c
      call multitop(13,16,2,2,'eta','mub/bin','LIN')
      call mtfill(17,20,2,2,'eta','mub/bin','LIN')
      call mtfill(21,24,2,2,'eta','mub/bin','LIN')
c
      call multitop(25,28,2,2,'dphi','mub/bin','LIN')
      call mtfill(29,32,2,2,'dphi','mub/bin','LIN')
      call mtfill(33,36,2,2,'dphi','mub/bin','LIN')
c
      return
      end


c--------------------------------------------------
      subroutine outfun(www,x1)
c This is the user analysis routine. It is called for each
c generated event with the parameters
c www: weight of the event (microbarn)
c x1: momentum fraction of the photon relative to the maximum energy of the
c     photon beam (therefore, x1=1 if the incoming particle is a monochromatic
c     photon, and x1=E_gamma/E_e if the incoming particle is an electron). 
c     WARNING!!!! This can only be used when running the pontlike code.
c     With the hadronic code in the e-p mode, there is no information
c     available on E_gamma/E_e, since the convolution of the 
c     Weizsaecker-Williams function with the photon parton densities
c     is carried out before running this code. In the case of the
c     hadronic component, x1 is set to -1
c The kinematical variables of the final state partons are given in common 
c blocks. The user can compute here the needed physical observables and add 
c the corresponding weight to the appropriate histogram.
c--------------------------------------------------
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c !!!! Do not alter the following variables !!!!
c xsc_min2 = minimum transverse energy squared
c xlam = lambda_QCD (5 flavours)
c xmuf2 = factorization scale squared
c xmur2 = renormalization scale squared
c xmues2 = Ellis-Sexton scale squared
c xmuww2 = Weizsaecker-Williams scale squared
c zg = strong coupling = sqrt(4*pi*alfas)
c ze2 = electron electric charge squared
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
c xkt_cm(i)=modulus of the transverse momentum of parton # i (i=1,2,3) in GeV
c xeta_cm(i)=pseudorapidity of parton # i (i=1,2,3)
c xphi_cm(i)=azimuthal angle of parton # i (i=1,2,3)
c These quantities are defined in the center-of-mass frame of the
c colliding particles (proton-antiproton at the Tevatron, proton-proton
c at the LHC, photon-proton or electron-proton system at HERA). 
c
c                      ******** WARNING ********
c
c In the case of HEAVY-ION COLLISIONS, the above quantities are given
c in the LAB FRAME.
c In all case, when needed, the user can perform a boost to an 
c arbitrary frame in this routine
      common/cmkin/xkt_cm,xeta_cm,xphi_cm
c sh is the square of the total CM energy of the colliding system (in GeV^2)
      common/shadr/sh
c imtt and imtt2 are matrices used in the jet-finding algorithm
c implemented in this version of the user file. They can be removed
c when the user implements his own version of the jet definition
      common/imtt/imtt
      common/imtt2/imtt2
      dimension xkt_cm(1:3),xeta_cm(1:3),xphi_cm(1:3)
      dimension ykt(1:3),yeta(1:3),yphi(1:3)
      dimension zkt(1:3),zeta(1:3),zphi(1:3)
      integer imtt(3:4,3:5)
      integer imtt2(3:5,3:5)
c proton and electron energies for HERA ( ==> E_cm=300.33 GeV)
      data e_p/820.d0/,e_e/27.5/
      data e_cm/300.33/
c
      e = sqrt(sh)
c Energy (in the lab frame) of the particle colliding against the proton.
c In the e-p mode, this energy is the energy of the incoming electron;
c in the gamma-p mode, it is the energy of the incoming monochromatic photon
      e_gamma = sh/(4*e_p)
c     write(*,*)'RADEK', e_gamma
c Boost to the lab frame (HERA physics)
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c We have the PHOTON in the POSITIVE direction.
      ycm_lab=-log(2*e_p/e)
c In the case is which the user is not dealing with HERA physics,
c and he is not interested in boosting the system, he should just set
c
c      ycm_lab=0.d0
c
      et_ev=0.d0
      do i=1,3
        ykt(i)=xkt_cm(i)
        yeta(i)=xeta_cm(i)+ycm_lab
        yphi(i)=xphi_cm(i)
        et_ev=et_ev+ykt(i)
      enddo
c At this stage, ykt(*), yeta(*), yphi(*) are the moduli of transverse 
c momenta, the pseudorapidities and the azimuthal angles of the outgoing 
c partons in the lab frame. If the user wants to use HERA conventions
c (that is, the photon is coming from the right) he can transform
c yeta(*) --> -yeta(*); yphi(*) --> yphi(*)+pi
c
c Jet-finding algorithm; in this case, the jets are defined by the
c cone algorithm with R=1.0. 
c The user can put here his own jet-finding routine
      call cone(ykt,yeta,yphi,zkt,zeta,zphi,1.d0)
c Now zkt(*),zeta(*),zphi(*) are the jet variables. For two-jet configurations,
c one of the entries of zkt is equal to zero. To fill the histograms, use
c mfill(hn,x,weight)
c where:
c hn = histogram number
c x = x value
c weight = weight of the event
      xmax=max(zkt(1),zkt(2),zkt(3))
      xmin=min(zkt(1),zkt(2),zkt(3))
      do i=1,3
        if(zkt(i).eq.xmax)n1=i
        if(zkt(i).eq.xmin)n3=i
      enddo
      n2=imtt2(n1+2,n3+2)-2
c now et(n3)<=et(n2)<=et(n1)
      do i=1,3
c E_t single inclusive, with eta cuts
        if(zeta(i).gt.-2.d0.and.zeta(i).lt.1.d0)
     #    call mfill(1,sngl(zkt(i)),sngl(www))
c eta single inclusive, with E_t cut
        if(zkt(i).gt.20.d0)
     #    call mfill(13,sngl(zeta(i)),sngl(www))
      enddo
c azimuthal distance between the two hardest jets
      if(zkt(n1).gt.15.d0.and.zkt(n2).gt.10.d0)then
        dphi=dacos(cos(zphi(n1)-zphi(n2)))
        call mfill(25,sngl(dphi),sngl(www))
      endif
c
c With the same partonic kinematics (ykt(*), yeta(*), yphi(*)) look for
c jets using a different jet-finding algorithm. Here the cone prescription 
c with R=0.7 is used
      call cone(ykt,yeta,yphi,zkt,zeta,zphi,0.7d0)
      xmax=max(zkt(1),zkt(2),zkt(3))
      xmin=min(zkt(1),zkt(2),zkt(3))
      do i=1,3
        if(zkt(i).eq.xmax)n1=i
        if(zkt(i).eq.xmin)n3=i
      enddo
      n2=imtt2(n1+2,n3+2)-2
c now et(n3)<=et(n2)<=et(n1)
      do i=1,3
c E_t single inclusive, with eta cuts
        if(zeta(i).gt.-2.d0.and.zeta(i).lt.1.d0)
     #    call mfill(5,sngl(zkt(i)),sngl(www))
c eta single inclusive, with E_t cut
        if(zkt(i).gt.20.d0)
     #    call mfill(17,sngl(zeta(i)),sngl(www))
      enddo
c azimuthal distance between the two hardest jets
      if(zkt(n1).gt.15.d0.and.zkt(n2).gt.10.d0)then
        dphi=dacos(cos(zphi(n1)-zphi(n2)))
        call mfill(29,sngl(dphi),sngl(www))
      endif
c
c With the same partonic kinematics (ykt(*), yeta(*), yphi(*)) look for
c jets using a different jet-finding algorithm. Here a k_T clustering
c algorithm with D=1 is used
      call esjet(ykt,yeta,yphi,zkt,zeta,zphi,1.d0)
      xmax=max(zkt(1),zkt(2),zkt(3))
      xmin=min(zkt(1),zkt(2),zkt(3))
      do i=1,3
        if(zkt(i).eq.xmax)n1=i
        if(zkt(i).eq.xmin)n3=i
      enddo
      n2=imtt2(n1+2,n3+2)-2
c now et(n3)<=et(n2)<=et(n1)
      do i=1,3
c E_t single inclusive, with eta cuts
        if(zeta(i).gt.-2.d0.and.zeta(i).lt.1.d0)
     #    call mfill(9,sngl(zkt(i)),sngl(www))
c eta single inclusive, with E_t cut
        if(zkt(i).gt.20.d0)
     #    call mfill(21,sngl(zeta(i)),sngl(www))
      enddo
c azimuthal distance between the two hardest jets
      if(zkt(n1).gt.15.d0.and.zkt(n2).gt.10.d0)then
        dphi=dacos(cos(zphi(n1)-zphi(n2)))
        call mfill(33,sngl(dphi),sngl(www))
      endif
c
      return
      end


c-------------------------------------------------------------------------
      function zgmu2()
c sets the mass scales and returns the strong coupling squared
c-------------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c xsc_min2 = minimum transverse energy squared
c xlam = lambda_QCD (5 flavours)
c xmuf2 = factorization scale squared
c xmur2 = renormalization scale squared
c xmues2 = Ellis-Sexton scale squared
c xmuww2 = Weizsaecker-Williams scale squared
c q2ww = Weizsaecker-Williams scale squared
c zg = strong coupling = sqrt(4*pi*alfas)
c ze2 = electron electric charge squared
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/q2ww/q2ww
c sclf2=xmuf2/(reference scale)**2
c sclr2=xmur2/(reference scale)**2
c scles2=xmues2/(reference scale)**2
      common/scalef/sclf2,sclr2,scles2
c xkt=modulus of the transverse momentum of parton # i (i=3,4,5)
c xeta=pseudorapidity of parton # i (i=3,4,5)
c xphi=azimuthal angle of parton # i (i=3,4,5)
c Define the reference scale in terms of these quantities
c WITH A BOOST-INVARIANT PRESCRIPTION
c NOTICE that the index i runs over the values 3,4,5 and NOT over 1,2,3
c as in the subroutine OUTFUN. xkt(*), xeta(*) and xphi(*) also have
c entry # 6. DO NOT USE IT
      common/partkin/xkt,xeta,xphi
c icls=0 ==> final state parton are well separated in the phase space
c icls=1 ==> parton # 3 and parton # 4 are close to each other.
c Two partons are close to each other in the sense of the definition
c adopted for the P-functions. Although not completely accurate, one
c can assume that, when a jet is obtained by the merging of two partons,
c these partons are labeled, in this function ZGMU2, as 3 and 4
      common/icls/icls
c number of light flavours
      common/nl/nl
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
c
c one possible definition
      xpt = (xkt(3)+xkt(4)+xkt(5))/2.d0
c another possible definition
c      if(icls.eq.0)then
c        xpt = max(xkt(3),xkt(4),xkt(5))
c      else
c        xpt = max(xkt(3)+xkt(4),xkt(5))
c      endif
c reference scale squared
      xmu2 = xpt**2
c factorization scale squared
      xmuf2 = sclf2 * xmu2
c renormalization scale squared
      xmur2 = sclr2 * xmu2
c Ellis-Sexton scale squared
      xmues2 = scles2 * xmu2
c Weiszacker-Williams scale squared
      xmuww2 = xmu2
      q2ww = xmuww2
c alpha_s at xmur2
      as = alfas(xmur2,xlam,nl)
      zgmu2 = 4.d0*pi*as
      zg = sqrt(zgmu2)
      end


c-----------------------------------------------------------------------
c ********************* JET-FINDING ROUTINES ***************************
c-----------------------------------------------------------------------
c
c Cone algorithm
c
      subroutine cone(ykt,yeta,yphi,zkt,zeta,zphi,r_cone)
c This subroutine returns jet variables (zkt, zeta, zphi) defined
c in terms of partonic variables (ykt, yeta, yphi) in the cone
c algorithm. 
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/rcone/rcone
      common/imtt2/imtt2
      dimension ykt(1:3),yeta(1:3),yphi(1:3)
      dimension zkt(1:3),zeta(1:3),zphi(1:3)
      dimension xr(1:2,2:3),xp(1:2,2:3)
      integer imtt2(3:5,3:5)
c
      rcone=r_cone
      i2bd=1
      do i=1,3
        if(ykt(i).eq.0.d0)i2bd=0
      enddo
      if(i2bd.eq.0)then
c 2-body kinematics
        do i=1,3
          zkt(i)=ykt(i)
          zeta(i)=yeta(i)
          zphi(i)=yphi(i)
        enddo
      else
c 3-body kinematics
        xr(1,2)=rdist(yeta(1),yphi(1),yeta(2),yphi(2))
        xp(1,2)=ptdist(ykt(1),ykt(2))
        xr(1,3)=rdist(yeta(1),yphi(1),yeta(3),yphi(3))
        xp(1,3)=ptdist(ykt(1),ykt(3))
        xr(2,3)=rdist(yeta(2),yphi(2),yeta(3),yphi(3))
        xp(2,3)=ptdist(ykt(2),ykt(3))
        imerge=0
        do i=1,2
          do j=i+1,3
            if(xr(i,j).lt.xp(i,j))then
              if(imerge.eq.1)then
                write(6,*)'Fatal error in subroutine cone'
                stop
              endif
              imerge=1
              n1=i
              n2=j
            endif
          enddo
        enddo
        if(imerge.eq.0)then
c no merging
          do i=1,3
            zkt(i)=ykt(i)
            zeta(i)=yeta(i)
            zphi(i)=yphi(i)
          enddo
        else
          n3=imtt2(n1+2,n2+2)-2
          zkt(n3)=ykt(n3)
          zeta(n3)=yeta(n3)
          zphi(n3)=yphi(n3)
          call xmerge(ykt(n1),ykt(n2),tmpkt,
     #                yeta(n1),yeta(n2),tmpeta,
     #                yphi(n1),yphi(n2),tmpphi)
          zkt(n1)=tmpkt
          zeta(n1)=tmpeta
          zphi(n1)=tmpphi
          zkt(n2)=0.d0
          zeta(n2)=1.d8
          zphi(n2)=0.d0
        endif
      endif
      return
      end


      function rdist(eta1,phi1,eta2,phi2)
c
c Square root of eq. (2.22)
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
      dteta=eta1-eta2
      dtphi=phi1-phi2
      dtphi=dacos(dcos(dtphi))
      rdist=dsqrt(dteta**2+dtphi**2)
      return
      end


      function ptdist(pt1,pt2)
c
c (pt1+pt2)*rcone/max(pt1,pt2), for the cone algorithm
c
      implicit real * 8 (a-h,o-z)
      common/rcone/rcone
c
      tmp=(pt1+pt2)/max(pt1,pt2)
      ptdist=tmp*rcone
      return
      end
c
c End of the cone algorithm
c
c 
c k_T clustering algorithm
c
      subroutine esjet(ykt,yeta,yphi,zkt,zeta,zphi,d_jet)
c This subroutine returns jet variables (zkt, zeta, zphi) defined
c in terms of partonic variables (ykt, yeta, yphi) in the Ellis-Soper
c algorithm. The parameter D must be set in the function ddist (in this file)
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/djet/djet
      common/imtt/imtt
      common/imtt2/imtt2
      dimension ykt(1:3),yeta(1:3),yphi(1:3)
      dimension zkt(1:3),zeta(1:3),zphi(1:3)
      dimension xd(1:3),xd2(1:2,2:3)
      integer imtt(3:4,3:5)
      integer imtt2(3:5,3:5)
c
      djet=d_jet
      i2bd=1
      do i=1,3
        if(ykt(i).eq.0.d0)i2bd=0
      enddo
      if(i2bd.eq.0)then
c 2-body kinematics
        do i=1,3
          zkt(i)=ykt(i)
          zeta(i)=yeta(i)
          zphi(i)=yphi(i)
        enddo
      else
c 3-body kinematics
        xd(1)=ykt(1)**2
        xd(2)=ykt(2)**2
        xd(3)=ykt(3)**2
        xd2(1,2)=ddist(ykt(1),yeta(1),yphi(1),ykt(2),yeta(2),yphi(2))
        xd2(1,3)=ddist(ykt(1),yeta(1),yphi(1),ykt(3),yeta(3),yphi(3))
        xd2(2,3)=ddist(ykt(2),yeta(2),yphi(2),ykt(3),yeta(3),yphi(3))
        xmin=min(xd(1),xd(2),xd(3),xd2(1,2),xd2(1,3),xd2(2,3))
        imerge=-1
        do i=1,3
          if(xd(i).eq.xmin)then
            imerge=0
            n1=i
          endif
        enddo
        if(imerge.eq.-1)then
          do i=1,2
            do j=i+1,3
              if(xd2(i,j).eq.xmin)then
                imerge=1
                n1=i
                n2=j
              endif
            enddo
          enddo
        endif
        if(imerge.eq.0)then
c no merging
          n2=imtt(3,n1+2)-2
          n3=imtt(4,n1+2)-2
          if(xd2(n2,n3).lt.min(xd(n2),xd(n3)))then
            write(6,*)'This is S2 contribution'
            stop
          endif
          do i=1,3
            zkt(i)=ykt(i)
            zeta(i)=yeta(i)
            zphi(i)=yphi(i)
          enddo
        elseif(imerge.eq.1)then
          n3=imtt2(n1+2,n2+2)-2
          zkt(n3)=ykt(n3)
          zeta(n3)=yeta(n3)
          zphi(n3)=yphi(n3)
          call xmerge(ykt(n1),ykt(n2),tmpkt,
     #                yeta(n1),yeta(n2),tmpeta,
     #                yphi(n1),yphi(n2),tmpphi)
          zkt(n1)=tmpkt
          zeta(n1)=tmpeta
          zphi(n1)=tmpphi
          zkt(n2)=0.d0
          zeta(n2)=1.d8
          zphi(n2)=0.d0
        else
          write(6,*)'Fatal error in subroutine esjet'
          write(6,*)'imerge=',imerge
          stop
        endif
      endif
      return
      end


      function rdist2(eta1,phi1,eta2,phi2)
c
c Eq. (2.22)
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
      dteta=eta1-eta2
      dtphi=phi1-phi2
      dtphi=dacos(dcos(dtphi))
      rdist2=dteta**2+dtphi**2
      return
      end


      function ddist(pt1,eta1,phi1,pt2,eta2,phi2)
c
c Eq. (2.23)
c
      implicit real * 8 (a-h,o-z)
      common/djet/djet
c
      tmp=min(pt1**2,pt2**2)*rdist2(eta1,phi1,eta2,phi2)
      ddist=tmp/djet**2
      return
      end
c
c End of the k_T clustering algorithm
c
      subroutine xmerge
     #    (xkt1,xkt2,xkt3,xeta1,xeta2,xeta3,xphi1,xphi2,xphi3)
c
c Defines the jet variables (xkt3,xeta3,xphi3) in terms of the
c partonic variables when two partons are merged
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
      xkt3=xkt1+xkt2
      xeta3=(xkt1*xeta1+xkt2*xeta2)/xkt3
      if(abs(xphi1-xphi2).lt.pi)then
        xphi3=(xkt1*xphi1+xkt2*xphi2)/xkt3
      else
        if(xphi2.lt.xphi1)then
          xphi3=(xkt1*xphi1+xkt2*(xphi2+2*pi))/xkt3
        else
          xphi3=(xkt1*xphi1+xkt2*(xphi2-2*pi))/xkt3
        endif
      endif
      xphi3=atan2(sin(xphi3),cos(xphi3))
      return
      end
