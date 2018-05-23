      implicit real*8(a-h,o-z)
      parameter (nlf=5)
      parameter (narg=2,kfun=nlf+1)
      parameter (ix=100,iq=10)
      parameter (ih=4)
      common/cwwcphyspar/alfaem,xme,xme2
      common/cwwcparton/ip,ndns
      common/cwwckine/x,qq
      common/cwwcitype/itype
      common/cwwciscale/iscale
      common/cwwcthcpar/e_e,thc
      common/cwwcxmu2/xmu2
      common/cwwcthrs/zmin,zmax
      character * 2 scheme
      real*4 q2(iq),xco(ix),xqdum2(ix,iq,kfun,2),xqdum1(ix,iq,kfun)
      external fww_conv
c
      alfaem=1/137.d0
      xme=0.511d-3
      xme2=xme**2
      write(6,*)'We use: 1/alpha_em=',1/alfaem
      write(6,*)'        m_e=',xme,' GeV'
      write(6,*)'  '
      write(6,*)'enter 0 to use WW, log term only'
      write(6,*)'      1 to use WW, log and nonlog terms'
      write(6,*)'      2 to use WW, cut on angle'
      read(5,*)itype
      if(itype.eq.0)then
        iscmax=2
      elseif(itype.eq.1)then
        iscmax=1
        write(6,*)'enter the effective WW scale in GeV (upper limit'
        write(6,*)'of the absolute value of the photon virtuality)'
        read(5,*)xmu
        xmu2=xmu**2
      elseif(itype.eq.2)then
        iscmax=1
        write(6,*)'enter electron energy, theta_cut'
        read(5,*)e_e,thc
      else
        write(6,*)'non implemented'
        stop
      endif
      write(6,*)'enter 0 to integrate WW over 0<x_gamma<1'
      write(6,*)'      1 otherwise'
      read(5,*)itmp
      if(itmp.ne.0)then
        write(6,*)'enter x_gamma(min),x_gamma(max)'
        read(5,*)xg_min,xg_max
        xww_min=log(xg_min)
        xww_max=log(xg_max)
      else
        xww_min=-1.d10
        xww_max=0.d0
        xg_max=1.d0
      endif
1     write(*,*) 'enter set # for photon, n<0 to list PDFs'
      read(*,*)ndns
      if(ndns.lt.0)then
        call prntsf
        goto 1
      endif
      call pdfpar(ndns,ih,xlam,scheme,iret)
      if(iret.ne.0)goto 1
c
      err=.1d-7
c minimum and maximum value allowed for the argument of the photon PDFs
      zmin=0.d0
      zmax=1.d0
c to set the x grid, use the following
      xmin=5.d-5
      xmax=min(1.d0,xg_max)
      xthrs=xmax
      if(xmax.eq.1.d0)xmax=0.99d0
c ymax is a parameter which is only used for sampling close to xmax
      ymax=xmax*1.01d0
      tmin=log(xmin/(ymax-xmin))
      tmax=log(xmax/(ymax-xmax))
c
      q2fix=5.d0
      do iscale=1,iscmax
c loop on the x points
        do jx=1,ix
          t = tmin + (tmax-tmin)*dfloat(jx-1)/dfloat(ix-1)
          z = exp(t)
          x = ymax*z/(1+z)
          xco(jx)=sngl(x)
          xl=max(log(x),xww_min)
          xh=min(log(1-1.d-6),xww_max)
c loop on the q points
          do jq=1,iq
            q2l10 = (q2fix*jq)/iq
            qq=10.d0**q2l10
            q2(jq)=sngl(qq)
c loop on the flavors: gluon --> bottom (nlf is set equal to 5)
            do ip=0,nlf
              jfun=ip+1
              if(xl.lt.xh)then
                val=dgauss(fww_conv,xl,xh,err)
                if(iscmax.eq.2)then
                  xqdum2(jx,jq,jfun,iscale)=sngl(log((val*x)+1))
                else
                  xqdum1(jx,jq,jfun)=sngl(log((val*x)+1))
                endif
              else
                if(iscmax.eq.2)then
                  xqdum2(jx,jq,jfun,iscale)=0.e0
                else
                  xqdum1(jx,jq,jfun)=0.e0
                endif
              endif
            enddo
          enddo
          qsqmin=q2(1)
          qsqmax=q2(iq)
        enddo
      enddo
c
      open(unit=10,file='el_pdf.for',status='new')
c
      write(10,'(a,i3)')'C PHOTON SET NUMBER = ',ndns
      if(itype.eq.0)then
        write(10,'(a)')'C WW, ONLY LOG TERM'
      elseif(itype.eq.1)then
        write(10,'(a)')'C WW, LOG AND NONLOG TERMS'
        write(10,'(a,e10.4)')'C Q^2=',sngl(xmu2)
      elseif(itype.eq.2)then
        write(10,'(a)')'C WW, CUT ON ANGLE'
        write(10,'(a,e10.4,a,e10.4)')'C E_EL=',sngl(e_e),
     #      ' THETA=',sngl(thc)
      else
        write(6,*)'error in main'
        stop
      endif
      if(itmp.eq.0)then
        write(10,'(a)')'C NO CUTS ON X_GAMMA'
      else
        write(10,'(a,e10.4,a,e10.4)')'C ',sngl(xg_min),
     #      ' < X_GAMMA < ',sngl(xg_max)
      endif
c
      write(10,'(a)')
     #'      SUBROUTINE ELPDF_USER(QSTAR2,X,FX,NF)'
C------------------------------------------------------------
C-------------   X is x and QSTAR2 is Q square  -------------
C------------------------------------------------------------
      write(10,'(a,i3,a,i3,a,i1,a)')
     #'      PARAMETER (IX=',ix,',IQ=',iq,',NARG=2,KFUN=',kfun,')'
C------------------------------------------------------------
C-----------Input tables are given as IX x points at IQ Q**2
C-----------positions (log spaced in Q**2 and x)
      write(10,'(a)')
     #'      DIMENSION FX(-NF:NF)',
     #'      DIMENSION ARG(NARG),NENT(NARG),ENT(IX+IQ)',
     #'      DIMENSION Q2(IQ),XCO(IX)'
      if(iscmax.eq.2)then
        write(10,'(a)')
     #'      DIMENSION XQDUM(IX,IQ,KFUN,2)',
     #'      REAL*8 Q2WW'
      else
        write(10,'(a)')
     #'      DIMENSION XQDUM(IX,IQ,KFUN)'
      endif
      write(10,'(a)')
     #'      REAL*8 XMIN,XMAX,XTHRS,QSQMIN,QSQMAX,QSQ',
     #'      REAL*8 IXMIN,IXMAX,IQSQMIN,IQSQMAX'
      if(iscmax.eq.2)then
        write(10,'(a)')
     #'      COMMON/Q2WW/Q2WW'
      endif
C------------------------------------------------------------
C----------- Next the IQ Q**2 of the input tables are defined
C------------------------------------------------------------
      write(10,'(a)')'      DATA Q2/'
      ir=1+(iq-1)/6
      do i=1,ir-1
         write(10,'(a,6(e10.4,1h,))')'     #',(q2(j),j=(i-1)*6+1,i*6)
      enddo
      write(10,'(a,6(e10.4,a))')
     #'     #',(q2(j),',',j=(ir-1)*6+1,iq-1),q2(iq),'/'
C----------------------------------------------------------------
C----------- And now the IX x positions  ------------------------
C----------------------------------------------------------------
      write(10,'(a)')'      DATA XCO/'
      ir=1+(ix-1)/6
      do i=1,ir-1
         write(10,'(a,6(e10.4,1h,))')'     #',(xco(j),j=(i-1)*6+1,i*6)
      enddo
      write(10,'(a,6(e10.4,a))')
     #'     #',(xco(j),',',j=(ir-1)*6+1,ix-1),xco(ix),'/'
C----------------------------------------------------------------
C----------- And now the array to interpolate -------------------
C----------------------------------------------------------------
      do iscale=1,iscmax
      do jfun=1,kfun
      do jq=1,iq
      if(iscmax.eq.2)then
        write(10,'(a,i3,a,i1,a,i1,a,i3,a)')
     #'      DATA (XQDUM(JX,',jq,',',jfun,',',iscale,'),JX=1,',IX,')/'
      else
        write(10,'(a,i3,a,i1,a,i3,a)')
     #'      DATA (XQDUM(JX,',jq,',',jfun,'),JX=1,',IX,')/'
      endif
c
      ir=1+(ix-1)/6
      if(iscmax.eq.2)then
        do i=1,ir-1
          write(10,'(a,6(e10.4,1h,))')
     #'     #',(xqdum2(j,jq,jfun,iscale),j=(i-1)*6+1,i*6)
        enddo
        write(10,'(a,6(e10.4,a))')
     #'     #',(xqdum2(j,jq,jfun,iscale),',',j=(ir-1)*6+1,ix-1),
     #          xqdum2(ix,jq,jfun,iscale),'/'
      else
        do i=1,ir-1
          write(10,'(a,6(e10.4,1h,))')
     #'     #',(xqdum1(j,jq,jfun),j=(i-1)*6+1,i*6)
        enddo
        write(10,'(a,6(e10.4,a))')
     #'     #',(xqdum1(j,jq,jfun),',',j=(ir-1)*6+1,ix-1),
     #          xqdum1(ix,jq,jfun),'/'
      endif
c
      enddo
      enddo
      enddo
      write(10,998)xthrs
998   format('      DATA XTHRS/',d10.5,'/')
      write(10,999)qsqmin,qsqmax
999   format('      DATA QSQMIN/',d10.5,'/,QSQMAX/',d10.5,'/')
      write(10,'(a)')
     #'      DATA XMIN/5.D-05/,XMAX/.99D0/',
     #'      DATA INI/0/',
     #'      IF(INI.EQ.0) THEN',
     #'        DO I=1,IX',
     #'          ENT(I)=LOG10(XCO(I))',
     #'        ENDDO',
     #'        NENT(1)=IX',
     #'        NENT(2)=IQ',
     #'        DO I=1,IQ',
     #'          ENT(IX+I)=LOG10(Q2(I))',
     #'        ENDDO',
     #'        ILXMIN=0',
     #'        ILXMAX=0',
     #'        ILQSQMIN=0',
     #'        ILQSQMAX=0',
     #'        INI=1',
     #'      ENDIF',
     #'      DX=DBLE(X)',
     #'      DQ=DBLE(SQRT(QSTAR2))'
      write(10,'(a)')
     #'      DO K=-NF,NF',
     #'        FX(K)=0',
     #'      ENDDO',
     #'      IF(DX.LT.XTHRS)THEN',
     #'        ARG(1)=LOG10(X)',
     #'        ARG(2)=LOG10(QSTAR2)'
      write(10,997)nlf
997   format('        JMAX=MIN(NF,',i1,')')
      if(iscmax.eq.2)then
        write(10,'(a)')
     #'        XLQ2WW = LOG(SNGL(Q2WW))',
     #'        DO J=0,JMAX',
     #'          FX(J) = FINT(NARG,ARG,NENT,ENT,XQDUM(1,1,J+1,1))',
     #'          IF(XLQ2WW.NE.0) THEN',
     #'            FX(J) = FX(J) +',
     #'     #        XLQ2WW*FINT(NARG,ARG,NENT,ENT,XQDUM(1,1,J+1,2))',
     #'          ENDIF',
     #'        ENDDO'
      else
        write(10,'(a)')
     #'        DO J=0,JMAX',
     #'          FX(J) = FINT(NARG,ARG,NENT,ENT,XQDUM(1,1,J+1))',
     #'        ENDDO'
      endif
      write(10,'(a)')
     #'        IF(DX.LT.XMIN) THEN',
     #'          IXMIN=IXMIN+1.',
     #'          IF(LOG10(IXMIN).GT.ILXMIN) THEN',
     #'            WRITE(*,*)''X < XMIN IN EL. PDFs MORE THAN 10**'',',
     #'     #    ILXMIN,'' TIMES''',
     #'            ILXMIN=ILXMIN+1',
     #'          ENDIF',
     #'        ENDIF',
     #'        IF(DX.GT.XMAX) THEN',
     #'          IXMAX=IXMAX+1.',
     #'          IF(LOG10(IXMAX).GT.ILXMAX) THEN',
     #'            WRITE(*,*)''X > XMAX IN EL. PDFs MORE THAN 10**'',',
     #'     #    ILXMAX,'' TIMES''',
     #'            ILXMAX=ILXMAX+1',
     #'          ENDIF',
     #'        ENDIF'
      write(10,'(a)')
     #'        QSQ=DQ**2',
     #'        IF(QSQ.LT.QSQMIN) THEN',
     #'          IQSQMIN=IQSQMIN+1.',
     #'          IF(LOG10(IQSQMIN).GT.ILQSQMIN) THEN',
     #'            WRITE(*,*)''Q**2 < MIN Q**2 IN EL. PDFs MORE '',',
     #'     #  ''THAN 10**'',ILQSQMIN,'' TIMES''',
     #'            ILQSQMIN=ILQSQMIN+1',
     #'          ENDIF',
     #'        ENDIF',
     #'        IF(QSQ.GT.QSQMAX) THEN',
     #'          IQSQMAX=IQSQMAX+1.',
     #'          IF(LOG10(IQSQMAX).GT.ILQSQMAX) THEN',
     #'            WRITE(*,*)''Q**2 > MAX Q**2 IN EL. PDFs MORE '',',
     #'     #  ''THAN 10**'',ILQSQMAX,'' TIMES''',
     #'            ILQSQMAX=ILQSQMAX+1',
     #'          ENDIF',
     #'        ENDIF'
      write(10,'(a)')
     #'        DO K=0,NF',
     #'          FX(K)=SNGL((EXP(DBLE(FX(K)))-1.D0)/DBLE(X))',
     #'          IF(FX(K).LT.0)FX(K)=0',
     #'        ENDDO',
     #'        DO K=1,NF',
     #'          FX(-K)=FX(K)',
     #'        ENDDO',
     #'      ENDIF',
     #'      END'
      end


      function fww_conv(eta)
      implicit real*8(a-h,o-z)
      parameter (ih=4)
      parameter(nlf=5)
      real*4 fx(-nlf:nlf),sz,sqq
      common/cwwcparton/ip,ndns
      common/cwwckine/x,qq
      common/cwwcthrs/zmin,zmax
c
c y is entering the WW function
      y=exp(eta)
c z is entering the photon PDFs
      z=x/y
      if(z.gt.zmax)z=zmax
      if(z.lt.zmin)z=zmin
      sz=sngl(z)
      sqq=sngl(qq)
      call mlmpdf(ndns,ih,sqq,sz,fx,nlf)
      fww_conv = fww_ww(y)*dble(fx(ip))
      return
      end


      function fww_ww(x)
      implicit real*8(a-h,o-z)
      parameter(pi=3.14159265358979312D0)
      common/cwwcphyspar/alfaem,xme,xme2
      common/cwwcitype/itype
      common/cwwciscale/iscale
      common/cwwcthcpar/e_e,thc
      common/cwwcxmu2/xmu2
c
      aemo2pi = alfaem/(2*pi)
      if(itype.eq.0)then
        if(iscale.eq.1) then
          fww_ww = aemo2pi*(1+(1-x)**2)/x*log( (1-x)/(x**2*xme2) )
        else
          fww_ww = aemo2pi*(1+(1-x)**2)/x
        endif
      elseif(itype.eq.1)then
        q2eff = xmu2
        tmp = (1+(1-x)**2) / x * log( q2eff*(1-x)/(xme*x)**2 )
     #        +2*xme2*x*( 1/q2eff-(1-x)/(xme*x)**2 )
        fww_ww = aemo2pi * tmp
      elseif(itype.eq.2)then
        q2eff = xme2*x**2+e_e**2*(1-x)**2*thc**2
        tmp = (1+(1-x)**2) / x * log( q2eff/(xme*x)**2 )
     #        +2*(1-x) * ( xme2*x/q2eff - 1/x )
        fww_ww = aemo2pi * tmp
      else
        write(6,*)'error in function fww_ww'
        stop
      endif
      return
      end      


      function dgauss(f,a,b,eps)
c.----------------------------------------------------------------------
c.
c.    gauss integral of the function f in interval a,b
c.    last update: 10/04/88
c.
c.----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension w(12),x(12)
      external f
      data const/1.e-12/
      data w
     &/0.101228536290376, 0.222381034453374, 0.313706645877887,
     & 0.362683783378362, 0.027152459411754, 0.062253523938648,
     & 0.095158511682493, 0.124628971255534, 0.149595988816577,
     & 0.169156519395003, 0.182603415044924, 0.189450610455069/
      data x
     &/0.960289856497536, 0.796666477413627, 0.525532409916329,
     & 0.183434642495650, 0.989400934991650, 0.944575023073233,
     & 0.865631202387832, 0.755404408355003, 0.617876244402644,
     & 0.458016777657227, 0.281603550779259, 0.095012509837637/
c--
c--   initialise
      delta=const*abs(a-b)
      dgauss=0.
      aa=a
c--
c--   iteration loop
   10 y=b-aa
c--
c--   epsilon reached ??
      if (abs(y).le.delta) return
   20 bb=aa+y
      c1=0.5*(aa+bb)
      c2=c1-aa
      s8=0.
      s16=0.
      do 30 i=1,4
        u=x(i)*c2
   30 s8=s8+w(i)*(f(c1+u)+f(c1-u))
      do 40 i=5,12
        u=x(i)*c2
   40 s16=s16+w(i)*(f(c1+u)+f(c1-u))
      s8=s8*c2
      s16=s16*c2
      if (abs(s16-s8).gt.eps*(1.0+abs(s16))) goto 50
      dgauss=dgauss+s16
      aa=bb
      goto 10
   50 y=0.5*y
      if (abs(y).le.delta) write(6,9040)
      goto 20
9040  format(1H ,'**** DGAUSS: Too high Accuracy required !!     ****')
      end
