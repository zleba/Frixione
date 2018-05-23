c**************************************************************************
c
c Classification of the partonic subprocesses:
c every process contributing to the physical cross section
c will be unambiguosly identified by the pair of integers
c
c                (jproc,jinst)
c
c Jproc identifies the unphysical matrix element 0 --> 5,
c while jinst identifies the initial state.
c We adopt the following conventions
c
c       jproc         unphysical process
c
c         1              0 --> 5g
c         2              0 --> 3g2q
c         3              0 --> 1g2q2Q
c         4              0 --> 1g4q
c
c
c       jinst         initial state
c
c         1                g g
c         2                q g
c         3                q qb
c         4                q q
c         5                q Qb
c         6                q Q
c
c
c We have therefore the following physical processes
c
c    (jproc,jinst)    physical process
c
c        (1,1)        g,g --> g,g,g
c
c        (2,1)        g,g --> q,qb,g
c        (2,2)        q,g --> q,g,g
c        (2,3)       q,qb --> g,g,g
c
c        (3,2)        q,g --> q,Q,Qb
c        (3,3)       q,qb --> Q,Qb,g
c        (3,5)       q,Qb --> q,Qb,g
c        (3,6)        q,Q --> q,Q,g
c
c        (4,2)        q,g --> q,q,qb
c        (4,3)       q,qb --> q,qb,g
c        (4,4)        q,q --> q,q,g
c
c 
c The leading order processes 0 --> 4 can be formally obtained
c from the 0 --> 5 ones by eliminating one gluon. We will mantain 
c the same classification scheme used for the 0 --> 5 processes. 
c Jproc=1 will mean the process 0 --> 4g, and so on.  The following 
c physical processes may occur in the leading order
c
c    (jproc,jinst)    physical process
c
c        (1,1)        g,g --> g,g
c
c        (2,1)        g,g --> q,qb
c        (2,2)        q,g --> q,g
c        (2,3)       q,qb --> g,g
c
c        (3,3)       q,qb --> Q,Qb
c        (3,5)       q,Qb --> q,Qb
c        (3,6)        q,Q --> q,Q
c
c        (4,3)       q,qb --> q,qb
c        (4,4)        q,q --> q,q
c
c**************************************************************************
c
c
c Begin of the 2 --> 3 routines
c
      function xmatel_five_part
     #     (s,xpp,xinv,xii,yi,xij,yj,jproc,jinst,ni,nj)
c This function, given the partonic CM energy (s) and the dot products of
c the four-momenta of the partons (xpp), returns the value of the 2 --> 3 
c matrix element (physical configuration) for the process specified 
c by (jproc,jinst)
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (tiny_s=1.d-7)
      parameter (tiny_c=1.d-8)
      common/i_type_part/i_type_part
      common/i_ap_in/i_ap_in1,i_ap_in2
      common/i_ap_out/i_ap_out
      common/jcoll_prc/j_prc_1_coll,j_prc_2_coll
      common/jcoll_ins/j_ins_1_coll,j_ins_2_coll
      common/jcoll_out_prc/j_prc_coll_out
      dimension xpp(1:5,1:5),ypp(1:5,1:5),xinv(1:5)
      integer i_type_part(1:4,1:6,3:5)
      integer i_ap_in1(1:4,1:6,3:5),i_ap_in2(1:4,1:6,3:5)
      integer i_ap_out(1:4,1:6,3:5,3:5)
      integer j_prc_1_coll(1:4,1:6,3:5),j_prc_2_coll(1:4,1:6,3:5)
      integer j_ins_1_coll(1:4,1:6,3:5),j_ins_2_coll(1:4,1:6,3:5)
      integer j_prc_coll_out(1:4,1:6,3:4,4:5)
c
      if(xii.lt.tiny_s)then
c Parton number ni is soft: use the analytical formula
c for the soft limit of the matrix element
        if(i_type_part(jproc,jinst,ni).ne.0)then
c The soft parton is not a gluon: the limit is finite, and the factor
c xii**2 in front of the matrix element gives therefore 0
          xmatel_five_part=0.d0
        else
c The soft parton is a gluon: use eq.(4.17)
          tmp=0.d0
          call x_xp_soft_yp(xpp,ypp,ni)
          neff=0
          do n=1,4
            if(n.ne.ni)then
              neff=neff+1
              meff=neff
            endif
            do m=n+1,5
              if(n.ne.ni.and.m.ne.ni)then
                meff=meff+1
                xfct1=eik_mni(xpp,xinv,yi,yj,n,m,ni,nj)
                xfct2=x_color_link(s,ypp,jproc,jinst,neff,meff)
                tmp=tmp+xfct1*xfct2
              endif
            enddo
          enddo
          xmatel_five_part=tmp/(8*pi**2)
        endif
      elseif(abs(1-yi).lt.tiny_c)then
c Parton number ni is collinear to parton 1: use the analytical formula
c for the collinear limit of the matrix element
        jproc0=j_prc_1_coll(jproc,jinst,ni)
        jinst0=j_ins_1_coll(jproc,jinst,ni)
        if(jproc0.eq.-1)then
c The splitting is not allowed: the limit is zero
          xmatel_five_part=0.d0
        else
c The splitting is allowed: use eq.(B.41)
          icode=i_ap_in1(jproc,jinst,ni)
c Momentum fraction entering the Altarelli-Parisi kernels
          x_ap=1-xii
c Reduced energy for the incoming partons
          s_red=s*x_ap
          xfct1=ap_kern(x_ap,abs(icode))
          call x_xp_collin_yp(s_red,xpp,ypp,ni,icode,1)
          xfct2=xmatel_four_part(s_red,ypp,jproc0,jinst0)
          xmatel_five_part=(4/s)*(1+yi)*xfct1*xfct2
        endif
      elseif(abs(1+yi).lt.tiny_c)then
c Parton number ni is collinear to parton 2: use the analytical formula
c for the collinear limit of the matrix element
        jproc0=j_prc_2_coll(jproc,jinst,ni)
        jinst0=j_ins_2_coll(jproc,jinst,ni)
        if(jproc0.eq.-1)then
c The splitting is not allowed: the limit is zero
          xmatel_five_part=0.d0
        else
c The splitting is allowed: use eq.(B.41)
          icode=i_ap_in2(jproc,jinst,ni)
c Momentum fraction entering the Altarelli-Parisi kernels
          x_ap=1-xii
c Reduced energy for the incoming partons
          s_red=s*x_ap
          xfct1=ap_kern(x_ap,abs(icode))
          call x_xp_collin_yp(s_red,xpp,ypp,ni,icode,-1)
          xfct2=xmatel_four_part(s_red,ypp,jproc0,jinst0)
          xmatel_five_part=(4/s)*(1-yi)*xfct1*xfct2
        endif
      elseif(abs(1-yj).lt.tiny_c)then
c Parton number ni is collinear to parton number nj: use the analytical 
c formula for the collinear limit of the matrix element
        n1=min(ni,nj)
        n2=max(ni,nj)
        jproc0=j_prc_coll_out(jproc,jinst,n1,n2)
        if(jproc0.eq.-1)then
c The splitting is not allowed: the limit is zero
          xmatel_five_part=0.d0
        else
c The splitting is allowed: use eq.(B.30)
          icode=i_ap_out(jproc,jinst,n1,n2)
          if(icode.eq.3.and.nj.lt.ni)then
            icode_ap=4
          else
            icode_ap=abs(icode)
          endif
c Momentum fraction entering the Altarelli-Parisi kernels
          x_ap=xij/(xii+xij)
          xfct1=ap_kern(x_ap,icode_ap)
          call x_xp_collout_yp(xpp,ypp,n1,n2,icode)
          xfct2=xmatel_four_part(s,ypp,jproc0,jinst)
          xmatel_five_part=(4/s)*(1/x_ap)*xfct1*xfct2
        endif
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element
        if(yi.eq.2.d0)then
c Factor multiplying the matrix element, eq.(4.65)
c The factor xij is inserted in the main program
          xfact=(1-yj)*xii**2
        elseif(yj.eq.2.d0)then
c Factor multiplying the matrix element, eq.(4.37)
          xfact=(1-yi**2)*xii**2
        else
          write(6,*)'xmatel_five_part called in the wrong way'
          stop
        endif
c Get the dot products to calculate the relevant matrix element
        call xpp_to_ypp(xpp,ypp,jproc,jinst)
        if(jproc.eq.1)then
          tmp=p5_ggggg(ypp)
        elseif(jproc.eq.2)then
          tmp=p5_qqbggg(ypp)
        elseif(jproc.eq.3)then
          tmp=p5_qpqpg(ypp)
        elseif(jproc.eq.4)then
          tmp=p5_qqqqg(ypp)
        else
          write(6,*)'Error in xmatel_five_part: unknown process number' 
          stop
        endif
        tmp=tmp/x_average(jinst)
        xmatel_five_part=xfact*tmp/(2*s)
      endif
      return
      end


      subroutine xmatel_coll(s,xpp,xi_ev,xi,y,jproc,jinst,
     #          jproc0,jinst0,ni,xmtel,xmtel_sc)
c This function, given a collinear configuration (xpp), returns
c the corresponding 2 --> 2 matrix elements times a factor depending
c upon the Altarelli-Parisi kernels and the K functions (which allow
c for a change of the factorization scheme), matching the definition
c of eq.(5.7). xmtel_sc contains the contribution K^{(reg)} (regular
c part of the K function) when xi#0, and K^{(d)} (the part of the
c K function proportional to delta(xi)) when xi=0.
      implicit real * 8 (a-h,o-z)
      character * 2 scheme
      parameter (pi=3.14159265358979312D0)
      parameter (one=1.d0)
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/cutpar/xicut,yincut,youtcut
      common/i_ap_in/i_ap_in1,i_ap_in2
      common/scheme/scheme
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      integer i_ap_in1(1:4,1:6,3:5),i_ap_in2(1:4,1:6,3:5)
c
c Momentum fraction entering the Altarelli-Parisi kernels
      x_ap=1-xi
c Reduced energy for the incoming partons
      s_red=s*x_ap
      if(y.eq.1.d0)then
c Parton number ni is collinear to parton number 1
        icode=i_ap_in1(jproc,jinst,ni)
        call x_xp_collin_yp(s_red,xpp,ypp,ni,icode,1)
      elseif(y.eq.-1.d0)then
c Parton number ni is collinear to parton number 2
        icode=i_ap_in2(jproc,jinst,ni)
        call x_xp_collin_yp(s_red,xpp,ypp,ni,icode,-1)
      else
        write(6,*)'error in xmatel_coll: wrong y value'
        stop
      endif
      xlgsdi_mu=log((s*yincut)/(2*xmuf2))
      xfct1=ap_kern(x_ap,abs(icode))
      xfct2=apprime_kern(x_ap,abs(icode))
      xfct3p=0.d0
      xfct3l=0.d0
      xfct4=xmatel_four_part(s_red,ypp,jproc0,jinst0)
      xfct5=0.d0
c
      if(scheme.eq.'DI')then
        xfct3p=xkplus(x_ap,abs(icode))
        xfct3l=xklog(x_ap,abs(icode))
        if(xi.ne.0.d0)then
          xfct5=xkreg(x_ap,abs(icode))
        else
          xfct5=xkdelta(abs(icode))
     #         +xkplus(one,abs(icode))*log(xicut)
     #         +xklog(one,abs(icode))*log(xicut)**2/2.d0
c This part contributes to sig2pr(soft), which is integrated in xi
c over the range (0,xicut). This implies the presence of a jacobian
c equal to xicut in the soft term, which has to be removed by hand
c in this case
          xfct5=xfct5/xicut
        endif
      elseif(scheme.ne.'MS')then
        write(6,*)'Error in xmatel_coll, y=',y
        write(6,*)'Factorization scheme ',scheme,' not known'
      endif
c
      xmtel=( xfct1*(xlgsdi_mu+2*log(xi_ev))-xfct2
     #       -xfct3p-xfct3l*log(xi_ev) )*xfct4
      xmtel_sc=-xfct5*xfct4
      return
      end


      function x_sgn(n,m)
c This function returns -1 if N or M is equal to 1 or 2, and 1 otherwise
      implicit real*8 (a-h,o-z)
c
      tmp=1.d0
      if(n.eq.1.or.n.eq.2)tmp=-tmp
      if(m.eq.1.or.m.eq.2)tmp=-tmp
      x_sgn=tmp
      return
      end


      function x_average(jinst)
c This function returns the quantity SIGN*OMEGA(A_1)*OMEGA(A_2), where
c the flavours A_1 and A_2 are equivalently given by JINST, using
c our classification for initial states. The quantity SIGN is -1 if JINST=2
c (<==> qg initial state), 1 otherwise, to take into account the sign coming
c from the crossing of a single fermionic line
      implicit real*8 (a-h,o-z)
      parameter (xnc=3.d0)
      parameter (vda=8.d0)
c
      xsgn=1.d0
      if(jinst.eq.1)then
        tmp=(2*vda)**2
      elseif(jinst.eq.2)then
        tmp=(2*vda)*(2*xnc)
        xsgn=-1.d0
      else
        tmp=(2*xnc)**2
      endif
      x_average=xsgn*tmp
      return
      end


      function eik_mni(xpp,xinv,yi,yj,n,m,ni,nj)
c This function evaluates the eikonal factors 
c         xp_n.xp_m/(xp_n.xp_i xp_m.xp_i)
c times the factor
c         xii**2 (1-yi**2)
c when we evaluate the limiting behaviour of eq.(4.37), or the factor
c         xii**2 (1-yj)
c when we evaluate the limiting behaviour of eq.(4.65). In the former case
c the function is called by fixing yj=2.d0; in the latter case,
c by fixing yi=2.d0
      implicit real*8 (a-h,o-z)
      dimension xpp(1:5,1:5),xinv(1:5)
c
      s=2*xpp(1,2)
      if(yj.eq.2.d0)then
        if(n.eq.1)then
          if(m.eq.2)then
            tmp=8.d0/s
          else
            tmp=8*x_sgn(n,m)*xpp(n,m)*(1+yi)/(s*sqrt(s)*xinv(m))
          endif
        elseif(n.eq.2)then
          tmp=8*x_sgn(n,m)*xpp(n,m)*(1-yi)/(s*sqrt(s)*xinv(m))
        else
          tmp=4*x_sgn(n,m)*xpp(n,m)*(1-yi**2)/(s*xinv(n)*xinv(m))
        endif
      elseif(yi.eq.2.d0)then
        xnum=x_sgn(n,m)*xpp(n,m)
        if(n.eq.nj.or.m.eq.nj)then
          tmp=8*xnum/(s*xinv(m)*xinv(n))
        else
          tmp=4*xnum*(1-yj)/(s*xinv(m)*xinv(n))
        endif
      else
        write(6,*)'Error in eik_mni: wrong yi or yj value'
        stop
      endif
      eik_mni=tmp
      return
      end


      function ap_kern(x,index)
c This function returns the quantity (1-x)*P_{ab}(x), where
c P_{ab} are the Altarelli-Parisi kernels, and the splitting partons
c {ab} are defined with the following conventions
c
c         index          ab
c
c           1            gg
c           2            qg
c           3            gq
c           4            qq
c
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
c
      if(index.eq.1)then
        ap_kern=2*vca*(x+(1-x)**2/x+x*(1-x)**2)
      elseif(index.eq.2)then
        ap_kern=vtf*(1-x)*(x**2+(1-x)**2)
      elseif(index.eq.3)then
        ap_kern=vcf*(1-x)*(1+(1-x)**2)/x
      elseif(index.eq.4)then
        ap_kern=vcf*(1+x**2)
      else
        write(6,*)'Error in ap_kern: wrong index value'
        stop
      endif
      return
      end


      function apprime_kern(x,index)
c This function returns the quantity (1-x)*P_{ab}^{prime}(x), where
c P_{ab}^{prime} is the ep-dependent part of the Altarelli-Parisi kernels, 
c and the codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
c
      if(index.eq.1)then
        apprime_kern=0.d0
      elseif(index.eq.2)then
        apprime_kern=-2*vtf*x*(1-x)**2
      elseif(index.eq.3)then
        apprime_kern=-vcf*(1-x)*x
      elseif(index.eq.4)then
        apprime_kern=-vcf*(1-x)**2
      else
        write(6,*)'Error in apprime_kern: wrong index value'
        stop
      endif
      return
      end


      function xkdelta(index)
c This function returns the quantity K^{(d)}_{ab}, relevant for
c the MS --> DIS change in the factorization scheme. 
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xkdelta=0.d0
      elseif(index.eq.2)then
        xkdelta=0.d0
      elseif(index.eq.3)then
        xkdelta=vcf*(9.d0/2.d0+pi**2/3.d0)
      elseif(index.eq.4)then
        xkdelta=-vcf*(9.d0/2.d0+pi**2/3.d0)
      else
        write(6,*)'Error in xkdelta: wrong index value'
        stop
      endif
      return
      end


      function xkplus(x,index)
c This function returns the quantity K^{(+)}_{ab}(x), relevant for
c the MS --> DIS change in the factorization scheme. Notice that
c there's NO multiplicative (1-x) factor like in the previous functions.
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xkplus=0.d0
      elseif(index.eq.2)then
        xkplus=0.d0
      elseif(index.eq.3)then
        xkplus=-vcf*(-3.d0/2.d0-(1+x**2)*log(x)+(1-x)*(3+2*x))
      elseif(index.eq.4)then
        xkplus=vcf*(-3.d0/2.d0-(1+x**2)*log(x)+(1-x)*(3+2*x))
      else
        write(6,*)'Error in xkplus: wrong index value'
        stop
      endif
      return
      end


      function xklog(x,index)
c This function returns the quantity K^{(l)}_{ab}(x), relevant for
c the MS --> DIS change in the factorization scheme. Notice that
c there's NO multiplicative (1-x) factor like in the previous functions.
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xklog=0.d0
      elseif(index.eq.2)then
        xklog=0.d0
      elseif(index.eq.3)then
        xklog=-vcf*(1+x**2)
      elseif(index.eq.4)then
        xklog=vcf*(1+x**2)
      else
        write(6,*)'Error in xklog: wrong index value'
        stop
      endif
      return
      end


      function xkreg(x,index)
c This function returns the quantity K^{(reg)}_{ab}(x), relevant for
c the MS --> DIS change in the factorization scheme. Notice that
c there's NO multiplicative (1-x) factor like in the previous functions.
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xkreg=-2*nl*vtf*( (x**2+(1-x)**2)*log((1-x)/x)+8*x*(1-x)-1 )
      elseif(index.eq.2)then
        xkreg=vtf*( (x**2+(1-x)**2)*log((1-x)/x)+8*x*(1-x)-1 )
      elseif(index.eq.3)then
        xkreg=0.d0
      elseif(index.eq.4)then
        xkreg=0.d0
      else
        write(6,*)'Error in xkreg: wrong index value'
        stop
      endif
      return
      end


      subroutine xpp_to_ypp(xpp,ypp,jproc,jinst)
c This subroutine fills the array ypp, in such a way that a 2 --> 3
c crossing invariant matrix element function called with entries ypp
c is eventually returning the value of the matrix element for the
c physical process defined by JPROC, JINST
      implicit real * 8 (a-h,o-z)
      common/i_phys/i_phys
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      integer i_phys(1:5,1:4,1:6)
c
      do i=1,5
        do j=1,5
          ypp(i,j)=xpp(i_phys(i,jproc,jinst),i_phys(j,jproc,jinst))
        enddo
      enddo
      return
      end


      subroutine xpp_to_ypp_0(xpp,ypp,jproc,jinst)
c This subroutine has the same meaning of XPP_TO_YPP, but is relevant
c for the 2 --> 2 processes
      implicit real * 8 (a-h,o-z)
      common/i0_phys/i0_phys
      dimension xpp(1:5,1:5),ypp(1:4,1:4)
      integer i0_phys(1:4,1:4,1:6)
c
      do i=1,4
        do j=1,4
          ypp(i,j)=xpp(i0_phys(i,jproc,jinst),i0_phys(j,jproc,jinst))
        enddo
      enddo
      return
      end


      subroutine x_xp_soft_yp(xpp,ypp,ni)
c This subroutine builds a 2 --> 2 kinematical configuration (ypp)
c starting from a 2 --> 3 kinematical configuration (xpp) in which the
c parton number ni is nearly soft. The relations
c
c      (12)+(13)+(14)=0,   (12)=(34),   (13)=(24),   (14)=(23)
c
c among the 2 --> 2 invariants are exploited
      implicit real * 8 (a-h,o-z)
      common/imtt/imtt
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      integer imtt(3:4,3:5)
c
      n3=imtt(3,ni)
      ypp(1,2)=xpp(1,2)
      ypp(1,3)=xpp(1,n3)
      ypp(1,4)=-ypp(1,2)-ypp(1,3)
      ypp(2,3)=ypp(1,4)
      ypp(2,4)=ypp(1,3)
      ypp(3,4)=ypp(1,2)
      do i=1,3
        do j=i+1,4
          ypp(j,i)=ypp(i,j)
        enddo
      enddo
      return
      end


      subroutine x_xp_collin_yp(s,xpp,ypp,ni,icode,iflag)
c This subroutine builds a 2 --> 2 kinematical configuration (ypp)
c starting from a 2 --> 3 kinematical configuration (xpp) which is
c nearly collinear (initial state). We use the following conventions
c
c  iflag                     icode
c
c    1   ==>   ni||1           >0   ==>   no action   
c   -1   ==>   ni||2           <0   ==>   ypp(1,*) <--> ypp(2,*)
c
      implicit real * 8 (a-h,o-z)
      common/imtt/imtt
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      integer imtt(3:4,3:5)
c
      n3=imtt(3,ni)
      ypp(1,2)=s/2
      if(iflag.eq.1)then
        if(icode.gt.0)then
          ypp(2,3)=xpp(2,n3)
          ypp(2,4)=-ypp(1,2)-ypp(2,3)
          ypp(1,3)=ypp(2,4)
          ypp(1,4)=ypp(2,3)
        else
          ypp(1,3)=xpp(2,n3)
          ypp(1,4)=-ypp(1,2)-ypp(1,3)
          ypp(2,3)=ypp(1,4)
          ypp(2,4)=ypp(1,3)
        endif
      elseif(iflag.eq.-1)then
        if(icode.gt.0)then
          ypp(1,3)=xpp(1,n3)
          ypp(1,4)=-ypp(1,2)-ypp(1,3)
          ypp(2,3)=ypp(1,4)
          ypp(2,4)=ypp(1,3)
        else
          ypp(2,3)=xpp(1,n3)
          ypp(2,4)=-ypp(1,2)-ypp(2,3)
          ypp(1,3)=ypp(2,4)
          ypp(1,4)=ypp(2,3)
        endif
      endif
      ypp(3,4)=ypp(1,2)
      do i=1,3
        do j=i+1,4
          ypp(j,i)=ypp(i,j)
        enddo
      enddo
      return
      end


      subroutine x_xp_collout_yp(xpp,ypp,ni,nj,icode)
c This subroutine builds a 2 --> 2 kinematical configuration (ypp)
c starting from a 2 --> 3 kinematical configuration (xpp) which is
c nearly collinear (final state). We use the following conventions 
c (remember that here ni<nj)
c
c   icode    splitting        action on momenta
c
c     1      g --> gg       yk(4,*)=xk(ni,*)+xk(nj,*)
c     2      g --> qqbar    yk(4,*)=xk(ni,*)+xk(nj,*)
c     3      q --> gq       yk(ni,*)=xk(ni,*)+xk(nj,*)
c
      implicit real * 8 (a-h,o-z)
      common/imtt2/imtt2
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      integer imtt2(3:5,3:5)
c
      n3=imtt2(ni,nj)
      ypp(1,2)=xpp(1,2)
      if(icode.eq.1.or.icode.eq.2)then
        ypp(1,3)=xpp(1,n3)
        ypp(1,4)=-ypp(1,2)-ypp(1,3)
      elseif(icode.eq.3)then
        nn=4*(4-ni)+3*(ni-3)
        ypp(1,nn)=xpp(1,n3)
        ypp(1,ni)=-ypp(1,2)-ypp(1,nn)
      else
        write(6,*)'Error in x_xp_collout_yp: wrong icode'
        stop
      endif
      ypp(2,3)=ypp(1,4)
      ypp(2,4)=ypp(1,3)
      ypp(3,4)=ypp(1,2)
      do i=1,3
        do j=i+1,4
          ypp(j,i)=ypp(i,j)
        enddo
      enddo
      return
      end
c
c Formulae for the cross sections
c
      function p5_ggggg(xpp)
c This function returns the 0 --> 5g matrix element. The spin and color
c average factors are NOT included
      implicit real * 8 (a-z)
      dimension xpp(1:5,1:5)
c
      x12 = xpp(1,2)
      x13 = xpp(1,3)
      x14 = xpp(1,4)
      x15 = xpp(1,5)
      x23 = xpp(2,3)
      x24 = xpp(2,4)
      x25 = xpp(2,5)
      x34 = xpp(3,4)
      x35 = xpp(3,5)
      x45 = xpp(4,5)
      xf1_3 = 1/(x12*x13*x24*x35*x45)+1/(x12*x14*x23*x35*x45)+1/(x12*x13
     1   *x25*x34*x45)+1/(x12*x15*x23*x34*x45)+1/(x13*x14*x23*x25*x45)+1
     2   /(x13*x15*x23*x24*x45)+1/(x12*x14*x25*x34*x35)+1/(x12*x15*x24*x
     3   34*x35)+1/(x13*x14*x24*x25*x35)+1/(x14*x15*x23*x24*x35)+1/(x13*
     4   x15*x24*x25*x34)+1/(x14*x15*x23*x25*x34)
      xf2 = x45**4+x35**4+x34**4+x25**4+x24**4+x23**4+x15**4+x14**4+x13*
     1   *4+x12**4
      p5_ggggg = 432*xf1_3*xf2
      return 
      end 


      function p5_qqbggg(xpp)
c This function returns the 0 --> 3g2q matrix element. The spin and color
c average factors are NOT included. The particle labeling is as follows
c
c         0 --> g(1) g(2) g(3) q(4) qb(5)
c
      implicit real*8 (a-h,o-z)
      dimension xpp(1:5,1:5)
c
      a1=xpp(4,1)                                            
      a2=xpp(4,2)                                            
      a3=xpp(4,3)                                            
      b1=xpp(5,1)                                            
      b2=xpp(5,2)                                            
      b3=xpp(5,3)                                            
      x12=xpp(1,2)                                           
      x23=xpp(2,3)                                           
      x31=xpp(3,1)                                           
      s=2.0*xpp(4,5)                                         
      f1=8.0/(4.0*81.0)                                           
      f2_3=(b3**2+a3**2)/(a1*a2*b1*b2)
     #    +(b2**2+a2**2)/(a1*a3*b1*b3)
     #    +(b1**2+a1**2)/(a2*a3*b2*b3)
      f4=(a1*b2+a2*b1)/x12+(a2*b3+a3*b2)/x23+(a3*b1+a1*b3)/x31    
      f4=9.0*(0.5*s-f4)                                           
      f5=a3*b3*(a1*b2+a2*b1)/(x23*x31)+                           
     $   a1*b1*(a2*b3+a3*b2)/(x31*x12)+                           
     $   a2*b2*(a3*b1+a1*b3)/(x12*x23)                            
      f5=f5*2.0*81.0/s                                            
      sum=f1*f2_3*(0.5*s+f4+f5)                                  
      p5_qqbggg=sum*4.0*9.0
      return                                                      
      end                                                         


      function p5_qpqpg(xpp)
c This function returns the 0 --> 1g2q2Q matrix element. The spin and color
c average factors are NOT included. The particle labeling is as follows
c
c         0 --> Q(1) Qb(2) g(3) q(4) qb(5)
c
      implicit real * 8 (a-z)
      parameter (xnc=3.d0)
      dimension xpp(1:5,1:5)
c
      s = 2*xpp(1,4)
      sp = 2*xpp(2,5)
      t = 2*xpp(1,2)
      tp = 2*xpp(4,5)
      u = 2*xpp(1,5)
      up = 2*xpp(2,4)
      x13 = xpp(1,3)
      x34 = xpp(3,4)
      x23 = xpp(2,3)
      x35 = xpp(3,5)
      xf1 = (up**2+u**2+sp**2+s**2)/(t*tp)
      xf2_3 = 2*(9*u+t)*x35/(x13*x23*x34)-2*t/(x34*x35)+2*(sp+s)/(x23*x3
     1   5)+14*u/(x13*x35)+2*(8*up+8*u+t)/(x23*x34)+2*(7*u-t)/(x13*x34)+
     2   18*u/(x13*x23)
      p5_qpqpg = (2*xnc)**2*xf1*xf2_3/54.d0
      return 
      end 


      function p5_qqqqg(xpp)
c This function returns the 0 --> 1g4q matrix element. The spin and color
c average factors are NOT included. The particle labeling is as follows
c
c         0 --> q(1) qb(2) g(3) q(4) qb(5)
c
      implicit real * 8 (a-z)
      parameter (xnc=3.d0)
      parameter (c1=64.d0/108.d0)
      parameter (c2=c1/8.d0)
      parameter (c3=10.d0/81.d0)
      parameter (c4=8.d0/81.d0)
      dimension xpp(1:5,1:5)
c
      s = 2*xpp(1,4)
      sp = 2*xpp(2,5)
      t = 2*xpp(1,2)
      tp = 2*xpp(4,5)
      u = 2*xpp(1,5)
      up = 2*xpp(2,4)
      x13 = xpp(1,3)
      x34 = xpp(3,4)
      x23 = xpp(2,3)
      x35 = xpp(3,5)
      xg12 = (up**2+u**2+sp**2+s**2)/(t*tp)
      xg34 = (tp**2+t**2+sp**2+s**2)/(u*up)
      xg56 = (sp**2+s**2)*(-u*up-t*tp+s*sp)/(t*tp*u*up)
      xf_3 = (xg56*(-c4*(2*up+tp+t+2*sp)/2.d0-c3*(2*sp+2*s)/2.d0)+
     1   c2*(sp+s)*xg34+c2*(sp+s)*xg12)/(x23*x35)+(xg56*(-c4*(up+u+tp+t)
     2   /2.d0-c3*(-2*u-2*t)/2.d0)+(c2*(-u-t)+c1*t)*xg34+(c1*u+c2*
     3   (-u-t))*xg12)/(x13*x34)+(xg56*(c3*u-c4*(up+u)/2.d0)+(c1*t-c2
     4   *u)*xg34+(c1*u-c2*u)*xg12)/(x13*x35)
      xf_3 = (xg56*(-c4*(up-u+tp+t-2*s)/2.d0-c3*u)+(c2*u+c1*t)*xg34+(
     1   c2*u+c1*u)*xg12)/(x13*x23)+(xg56*(c3*t-c4*(tp+t)/2.d0)+(c1*(
     2   tp+t)-c2*t)*xg34-c2*t*xg12)/(x34*x35)-c4*(up-u)*x34*xg56/(x13*x
     3   23*x35*2.d0)+((c4*(s-up)-c3*t)*xg56+c2*t*xg34+(c1*(up+u)+c2*
     4   t)*xg12)/(x23*x34)+x35*((c3*(-u-t)+c4*s)*xg56+c2*(u+t)*xg34+(c2
     5   *(u+t)+c1*u)*xg12)/(x13*x23*x34)+c1*t*x23*xg34/(x13*x34*x35)+xf
     6   _3
      p5_qqqqg = (2*xnc)**2*xf_3/2.d0
      return 
      end 
c
c End of the 2 --> 3 routines
c
c
c Begin of the 2 --> 2 routines
c
      function xmatel_2pv_contr(s,xpp,jproc,jinst)
c This function is relevant for the evaluation of eq.(5.5).
c It is defined such that
c
c  XMATEL_2PV_CONTR = Q({a_l}) M^{(3,0)}({a_l})
c                   + sum_{m<n} I_{mn}^{reg} M_{mn}^{(3,0)}({a_l})
c                   + M_{NS}^{(3,1)}({a_l})
c
c The sum over final state flavours is performed in the main program
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/cutpar/xicut,yincut,youtcut
      common/color_factors/xca,xgmm,xgmmprime
      common/i_type_part/i_type_part
      common/i_type_in/i_type_in
      dimension xpp(1:5,1:5),xca(0:1),xgmm(0:1),xgmmprime(0:1)
      dimension xe(1:4),xpnpm(1:3,2:4)
      integer i_type_part(1:4,1:6,3:5),i_type_in(1:6,1:2)
c
      xlgsd0_qes=log((s*youtcut)/(2*xmues2))
      xlgxi=log(xicut)
      xlgmu_qes=log(xmuf2/xmues2)
c Energy of the partons
      xe(1)=sqrt(s)/2.d0
      xe(2)=sqrt(s)/2.d0
      xe(3)=-(xpp(1,3)+xpp(2,3))/sqrt(s)
      xe(4)=-(xpp(1,4)+xpp(2,4))/sqrt(s)
      do i=1,3
        do j=i+1,4
          xpnpm(i,j)=x_sgn(i,j)*xpp(i,j)
        enddo
      enddo
      xq=0.d0
      do j=3,4
        ip=i_type_part(jproc,jinst,j)
        xq=xq+xgmmprime(ip)
     #    -xlgsd0_qes*( xgmm(ip)
     #                 -2*xca(ip)*log((2*xe(j))/(xicut*sqrt(s))))
     #    +2*xca(ip)*(log((2*xe(j))/sqrt(s))**2-xlgxi**2)
     #    -2*xgmm(ip)*log((2*xe(j))/sqrt(s))
      enddo
      i1=i_type_in(jinst,1)
      i2=i_type_in(jinst,2)
      xq=xq-xlgmu_qes*( xgmm(i1)+2*xca(i1)*xlgxi
     #                + xgmm(i2)+2*xca(i2)*xlgxi )
      xtmp1=xq*xmatel_four_part(s,xpp,jproc,jinst)
c
      xtmp2=0.d0
      do m=1,3
        do n=m+1,4
c In the case when partons n and m are back to back, xlg_4_mn is
c equal to zero, being xlg_4_mn=lim(ep-->0) [ log(2*ep)*log(1-ep/2) ],
c with cos(theta_{mn})=-1+ep ==> treat this case separately
          if((m.eq.1.and.n.eq.2).or.(m.eq.3.and.n.eq.4))then
            xlg_4_mn=0.d0
          else
            xlg_4_mn=log(4.d0-(2*xpnpm(m,n))/(xe(n)*xe(m)))
     #               *log(xpnpm(m,n)/(2*xe(n)*xe(m)))
          endif
          xinm=log(xicut**2*s/xmues2)**2/2.d0
     #        +log(xicut**2*s/xmues2)*log(xpnpm(m,n)/(2*xe(n)*xe(m)))
     #        -ddilog(xpnpm(m,n)/(2*xe(n)*xe(m)))
     #        +log((2*xpnpm(m,n))/(xe(n)*xe(m)))**2/2.d0
     #        -xlg_4_mn-2*log(2.d0)**2
          xtmp2=xtmp2+xinm*x_color_link(s,xpp,jproc,jinst,m,n)/(8*pi**2)
        enddo
      enddo
c
      xmatel_2pv_contr=xtmp1+xtmp2+x_loop_four(s,xpp,jproc,jinst)
      return
      end


      function xmatel_four_part(s,xpp,jproc,jinst)
c This function has the same meaning of xmatel_five_part, but
c it is relevant for the 2 --> 2 processes
      implicit real * 8 (a-h,o-z)
      dimension xpp(1:5,1:5),ypp(1:4,1:4)
      integer imap(1:4)
c Map present classification scheme for unphysical processes to
c the Kunszt-Soper scheme (1=2q2Q, 2=4q, 3=2g2q, 4=4g)
      data imap/4,3,1,2/
c
c Get the dot products to calculate the relevant matrix element
      call xpp_to_ypp_0(xpp,ypp,jproc,jinst)
      tmp=psitilde4(imap(jproc),ypp)/x_average(jinst)
      xmatel_four_part=tmp/(2*s)
      return
      end


      function x_color_link(s,xpp,jproc,jinst,n,m)
c This function gives the color-linked Born squared amplitudes;
c it is relevant for the soft limit of the 2 --> 3 processes. The
c factor 8 pi**2 inserted is consistent with the normalization
c of eq.(3.3)
      implicit real * 8 (a-h,o-z)
      real*8 lamtilde
      parameter (pi=3.14159265358979312D0)
      common/i_xcol/i_xcol
      dimension xpp(1:5,1:5),ypp(1:4,1:4)
      integer imap(1:4),i_xcol(1:4,1:4,1:6)
c Map present classification scheme for unphysical processes to
c the Kunszt-Soper scheme (1=2q2Q, 2=4q, 3=2g2q, 4=4g)
      data imap/4,3,1,2/
c
c Get the dot products to calculate the relevant matrix element
      call xpp_to_ypp_0(xpp,ypp,jproc,jinst)
      tmp=lamtilde(imap(jproc),i_xcol(n,jproc,jinst),
     #             i_xcol(m,jproc,jinst),ypp)/x_average(jinst)
      x_color_link=8*pi**2*tmp/(2*s)
      return
      end


      function x_loop_four(s,xpp,jproc,jinst)
c This function returns the M_{NS}^{(3,1)} part of the virtual
c contribution, as defined in eq. (3.2)
      implicit real * 8 (a-h,o-z)
      dimension xpp(1:5,1:5),ypp(1:4,1:4)
      integer imap(1:4)
c Map present classification scheme for unphysical processes to
c the Kunszt-Soper scheme (1=2q2Q, 2=4q, 3=2g2q, 4=4g)
      data imap/4,3,1,2/
c
c Get the dot products to calculate the relevant matrix element
      call xpp_to_ypp_0(xpp,ypp,jproc,jinst)
      tmp=psitilde6ns(imap(jproc),ypp)/x_average(jinst)
      x_loop_four=tmp/(2*s)
      return
      end


      function theta(flag)
c This function is the Heaviside function
      implicit real*8(a-z)
      logical flag
c
      if(flag)then
        tmp=1.d0
      else
        tmp=0.d0
      endif
      theta=tmp
      return
      end
c
c Formulae for the cross sections
c
      DOUBLE PRECISION FUNCTION PSITILDE4(PROCESS,K)
c This function is taken from the Kunszt-Soper package. The process
c number 1,2,3,4 correspond to process type A,B,C,D respectively.
c The particle labels are given in eq.(B.13) in Kunszt-Soper paper
      IMPLICIT NONE
      INTEGER PROCESS
      DOUBLE PRECISION K(4,4)
C
C      This function gives the psitilde^(4) functions from E.S.
C
      DOUBLE PRECISION S,T,U
      DOUBLE PRECISION N,V
      PARAMETER ( V = 8.0 D0)
      PARAMETER ( N = 3.0 D0)
C
      S = 2.0D0 * K(1,2)
      T = 2.0D0 * K(1,3)
      U = 2.0D0 * K(1,4)
C
      IF (PROCESS.EQ.1) THEN
       PSITILDE4 = 2.0D0 * V * (S**2 + U**2) / T**2
      ELSEIF (PROCESS.EQ.2) THEN
       PSITILDE4 = 2.0D0 * V * (S**2 + U**2) / T**2
     >           + 2.0D0 * V * (S**2 + T**2) / U**2
     >           - 4.0D0 * V / N * S**2 / U / T
      ELSEIF (PROCESS.EQ.3) THEN
       PSITILDE4 = 2.0D0 * V / N
     >             *( V / U / T - 2.0D0 * N**2 / S**2 )
     >             *( T**2 + U**2 )
      ELSEIF (PROCESS.EQ.4) THEN
       PSITILDE4 = 4.0D0 * V * N**2
     >          *( S**2 + T**2 + U**2 )
     >          /S**2 /T**2 /U**2
     >          *( S**4 + T**4 + U**4 )
      ELSE
       PSITILDE4 = 0.0D0
      ENDIF
      RETURN
      END


      DOUBLE PRECISION FUNCTION LAMTILDE(PROCESS,I1,I2,K)
c This function is taken from the Kunszt-Soper package. The process
c number 1,2,3,4 correspond to process type A,B,C,D respectively.
c The particle labels are given in eq.(B.13) in Kunszt-Soper paper
      IMPLICIT NONE
      INTEGER PROCESS
      INTEGER I1,I2
      DOUBLE PRECISION K(4,4)
C
C      This function gives the lambdatilde(a1,a2,s,t,u) functions
C
      INTEGER II1,II2
      DOUBLE PRECISION S,T,U
      DOUBLE PRECISION HAT,HAU,HB,HCT,HCU,HCS,HDS,HDT,HDU
      DOUBLE PRECISION N,V,CF
      PARAMETER( N = 3.0D0 )
      PARAMETER( V = 8.0D0 )
      PARAMETER( CF = 4.0D0 / 3.0D0 )
C
      S = 2.0D0 * K(1,2)
      T = 2.0D0 * K(1,3)
      U = 2.0D0 * K(1,4)
C
C  Make it a symmetric function of I1,I2
C
      IF (I1.LT.I2) THEN
        II1 = I1
        II2 = I2
      ELSE
        II1 = I2
        II2 = I1
      ENDIF
C
      IF (PROCESS.EQ.1) THEN
        HAT = 2.0D0 * (S**2 + U**2) /T**2
        IF     ( (II1.EQ.1).AND.(II2.EQ.2) ) THEN
          LAMTILDE = 4.0D0 * CF * HAT
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.3) ) THEN
          LAMTILDE = - 2.0D0 * CF * HAT
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 2.0D0 * CF * (N**2 - 2.0D0) * HAT
        ELSEIF ( (II1.EQ.3).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 4.0D0 * CF * HAT
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.4) ) THEN
          LAMTILDE = - 2.0D0 * CF* HAT
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.3) ) THEN
          LAMTILDE = 2.0D0 * CF * (N**2 - 2.0D0) * HAT
        ENDIF
C
      ELSEIF (PROCESS.EQ.2) THEN
        HAT = 2.0D0 * (S**2 + U**2) /T**2
        HAU = 2.0D0 * (S**2 + T**2) /U**2
        HB  = 2.0D0 * S**2 /T /U
        IF     ( (II1.EQ.1).AND.(II2.EQ.2) ) THEN
          LAMTILDE = 4.0D0 * CF /N
     >       *( - HB + N*(HAT + HAU) - N**2 * HB )
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.3) ) THEN
          LAMTILDE = 4.0D0 * CF /N
     >       *( HB - N*(HAU + 0.5D0*HAT) + 0.5D0 * N**3 * HAU )
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 4.0D0 * CF /N
     >       *( HB - N*(HAT + 0.5D0*HAU) + 0.5D0 * N**3 * HAT )
        ELSEIF ( (II1.EQ.3).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 4.0D0 * CF /N
     >       *( - HB + N*(HAT + HAU) - N**2 * HB )
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 4.0D0 * CF /N
     >       *( HB - N*(HAU + 0.5D0*HAT) + 0.5D0 * N**3 * HAU )
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.3) ) THEN
          LAMTILDE = 4.0D0 * CF /N
     >       *( HB - N*(HAT + 0.5D0*HAU) + 0.5D0 * N**3 * HAT )
        ENDIF
C
      ELSEIF (PROCESS.EQ.3) THEN
        HCT = 2.0D0 * (T**2 + U**2) /S**2 * T/U
        HCU = 2.0D0 * (T**2 + U**2) /S**2 * U/T
        HCS = 2.0D0 * (T**2 + U**2) /T /U
        IF     ( (II1.EQ.1).AND.(II2.EQ.2) ) THEN
          LAMTILDE = V * ( - HCT - HCU + HCS + HCS /N**2 )
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.3) ) THEN
          LAMTILDE = V * ( N**2 * HCU - HCS )
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.4) ) THEN
          LAMTILDE = V * ( N**2 * HCT - HCS )
        ELSEIF ( (II1.EQ.3).AND.(II2.EQ.4) ) THEN
          LAMTILDE =  V * N**2 * ( HCT + HCU )
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.4) ) THEN
          LAMTILDE = V * ( N**2 * HCU - HCS )
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.3) ) THEN
          LAMTILDE = V * ( N**2 * HCT - HCS )
        ENDIF
C
      ELSEIF (PROCESS.EQ.4) THEN
       HDS = 2.0D0 * (T**2 + U**2) /S**2 /T**2 /U**2
     >          * (S**4 + T**4 + U**4)
       HDT = 2.0D0 * (U**2 + S**2) /S**2 /T**2 /U**2
     >          * (S**4 + T**4 + U**4)
       HDU = 2.0D0 * (S**2 + T**2) /S**2 /T**2 /U**2
     >          * (S**4 + T**4 + U**4)
        IF     ( (II1.EQ.1).AND.(II2.EQ.2) ) THEN
          LAMTILDE = 2.0D0 * V * N**3 * HDS
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.3) ) THEN
          LAMTILDE = 2.0D0 * V * N**3 * HDT
        ELSEIF ( (II1.EQ.1).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 2.0D0 * V * N**3 * HDU
        ELSEIF ( (II1.EQ.3).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 2.0D0 * V * N**3 * HDS
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.4) ) THEN
          LAMTILDE = 2.0D0 * V * N**3 * HDT
        ELSEIF ( (II1.EQ.2).AND.(II2.EQ.3) ) THEN
          LAMTILDE = 2.0D0 * V * N**3 * HDU
        ENDIF
C
      ELSE
       LAMTILDE = 0.0D0
      ENDIF
      RETURN
      END


      DOUBLE PRECISION FUNCTION PSITILDE6NS(PROCESS,K)
c This function is taken from the Kunszt-Soper package. Modified on
c June, 17, 1996 by S. Frixione to take into account the new 
c structure of the common blocks
      IMPLICIT NONE
      INTEGER PROCESS
      DOUBLE PRECISION K(4,4)
C
C      This function gives the psitilde^(6)_NS functions calculated
C      by subtracting the singular pieces from from E.S.
C
      DOUBLE PRECISION S,T,U
      DOUBLE PRECISION THETA
      DOUBLE PRECISION LS,LT,LU,LMU,L2S,L2T,L2U
      DOUBLE PRECISION XSC_MIN2,XLAM,XMUF2,XMUR2,XMUES2,XMUWW2,ZG,ZE2
      DOUBLE PRECISION PI,N,V,NFLAVOR,MUUV,QES
      INTEGER NL
      COMMON/FIXVAR/XSC_MIN2,XLAM,XMUF2,XMUR2,XMUES2,XMUWW2,ZG,ZE2
      COMMON/NL/NL
      PARAMETER ( PI = 3.141592654 D0)
      PARAMETER ( V = 8.0 D0)
      PARAMETER ( N = 3.0 D0)
C
      NFLAVOR=DFLOAT(NL)
      MUUV=DSQRT(XMUF2)
      QES=DSQRT(XMUES2)
C
      S = 2.0D0 * K(1,2)
      T = 2.0D0 * K(1,3)
      U = 2.0D0 * K(1,4)
C
C  Here are the Ellis and Sexton log functions
C
      LS = DLOG(DABS(S/QES**2))
      LT = DLOG(DABS(T/QES**2))
      LU = DLOG(DABS(U/QES**2))
      LMU = DLOG(DABS(MUUV**2/QES**2))
      L2S = LS**2 - PI**2 * THETA(S.GT.0.0D0)
      L2T = LT**2 - PI**2 * THETA(T.GT.0.0D0)
      L2U = LU**2 - PI**2 * THETA(U.GT.0.0D0)
C
C  Find what type, then calculate
C
      IF (PROCESS.EQ.1) THEN
       PSITILDE6NS =
     >    V/(9.0D0*N*T**2)
     >    *(-  9*PI**2*(N**2-4 )*(S**2 - U**2)
     >      +  2*(72 + 13*N**2 - 10*N*NFLAVOR + 9*N**2*PI**2)
     >          *(S**2 + U**2)
     >      +  6*N*(11*N - 2*NFLAVOR)*(S**2 + U**2)*LMU
     >      - 36*T*U*LS
     >      +  3*(- 6*S**2 - 30*U**2 + 4*N*NFLAVOR*(S**2 + U**2)
     >                     - N**2*(7*S**2 + 3*T**2 + U**2))*LT
     >      - 18*S*T*(N**2 -2)*LU
     >      - 36*(3*S**2 + U**2)*LS*LT
     >      - 18*(N**2 -2)*(S**2 + 3*U**2)*LT*LU
     >      + 18*(S**2 - U**2)*L2S
     >      +  9*(N**2*(S**2 + 3*U**2) + 2*(3*S**2 - U**2))*L2T
     >      -  9*(N**2 - 2)*(S**2 - U**2)*L2U)
      ELSEIF (PROCESS.EQ.2) THEN
C
       PSITILDE6NS =
     >     40*N*NFLAVOR*S**2*T*U - 4*N**2*S**2*T*U*(13 + 9*PI**2)
     >   - 9*T*U*(32*S**2 + PI**2*(S**2 + 2*T*U))
     >   + 36*N*(S**2*(4 + PI**2) * (T**2 + U**2)
     >   + (4 - PI**2)  *  (T**4 + U**4))
     >   + N**3*(S**2*(26 + 9*PI**2)*(T**2 + U**2)
     >   + (26 + 27*PI**2)*(T**4 + U**4))
     >   - 20*N**2*NFLAVOR*(T**4 + U**4 + S**2*(T**2 + U**2))
C
       PSITILDE6NS =  PSITILDE6NS
     >  + 6*N*(11*N - 2*NFLAVOR)
     >     *(-2*S**2*T*U + N*(T**4 + U**4 + S**2*(T**2 + U**2))) *LMU
     >  - 36*N*T*U*(T**2 + U**2)*LS
     >  - 6*U*(- 2*N**2*S**2*T + 2*N*NFLAVOR*S**2*T
     >         - 2*N**2*NFLAVOR*U*(S**2 + U**2)
     >         + 6*S*T*(T + 2*U) + N**3*(-3*T**3 + 2*T**2*U
     >         + 7*T*U**2 + 4*U**3)
     >         + 3*N*(2*T**3 + 3*T**2*U + 2*T*U**2 + 6*U**3)) *LT
     >  - 6*T*(- 2*N**2*S**2*U + 2*N*NFLAVOR*S**2*U
     >         - 2*N**2*NFLAVOR*T*(S**2 + T**2)
     >         + 6*S*U*(2*T + U)
     >         + N**3*(4*T**3 + 7*T**2*U + 2*T*U**2 - 3*U**3)
     >         + 3*N*(6*T**3 + 2*T**2*U + 3*T*U**2 + 2*U**3)) * LU
C
       PSITILDE6NS =  PSITILDE6NS
     >  - 36*U*(- S**2*T - N**2*S**2*T + N*U*(3*S**2 + U**2)) *LS*LT
     >  - 36*T*(- S**2*U - N**2*S**2*U + N*T*(3*S**2 + T**2)) *LS*LU
     >  - 18*(  2*N*(-2 + N**2)*(T**2 - T*U + U**2)
     >          *(2*T**2 + 3*T*U + 2*U**2)
     >        + T*U*(3*T**2 + 4*T*U + 3*U**2)) * LT*LU
     >  + 18*N*T*U*(S**2 + T**2 + U**2) * L2S
     >  + 9*U*(- 2*N**2*S**2*T - 4*N*(S**3 - S*T**2 - T**3)
     >         - T*(3*T**2 + 8*T*U + 3*U**2)
     >         - 2*N**3*(T**3 - T*U**2 - 2*U**3)) * L2T
     >  + 9*T*(- 2*N**2*S**2*U - 4*N*(S**3 - S*U**2 - U**3)
     >         - U*(3*T**2 + 8*T*U + 3*U**2)
     >         + 2*N**3*(2*T**3 + T**2*U - U**3)) * L2U
C
       PSITILDE6NS = V/(9*N**2*T**2*U**2) * PSITILDE6NS
C
      ELSEIF (PROCESS.EQ.3) THEN
C
       PSITILDE6NS =
     >   3*S**2*( -7*(T**2 + U**2) + PI**2*(3*S**2 -2*T*U) )
     >   + 3*N**2*( 14*(T**4 + U**4) + T*U*(13*S**2 + 4*T*U )
     >            + PI**2*(T**3*U + 2*T**2*U**2 + T*U**3)    )
     >   - 3*N**4*( 7*(T**4 + U**4)+ T*U*(S**2 + 10*T*U)
     >            - PI**2*(T - U)**2*(T**2 + U**2)  )
     >   + 2*N*(11*N - 2*NFLAVOR)*(-S**2 + N**2*(T**2 + U**2))
     >         *(T**2 + U**2) * LMU
     >   - 3*T*U*(  4*T**2 + 8*T*U + 4*U**2
     >         + N**2*(-1 + N**2)*(T**2 - 10*T*U + U**2) ) * LS
     >   + 3*S*U*(S + N**2*U)*(2*T + 3*U + N**2*(2*T - 3*U)) * LT
     >   + 3*S*T*(S + N**2*T)*(3*T + 2*U + N**2*(-3*T + 2*U)) * LU
C
       PSITILDE6NS = PSITILDE6NS +
     >   - 6*(  N**4*U*(U - T)*(2*T**2 + T*U + U**2)
     >     + S*(S - N**2*T)*(2*T**2 + 2*T*U + U**2)  ) * LS*LT
     >   - 6*(  N**4*T*(T - U)*(T**2 + T*U + 2*U**2)
     >   + S*(S - N**2*U)*(T**2 + 2*T*U + 2*U**2)  ) * LS*LU
     >   + 12*N**2*S**2*(T**2 + U**2) * LT*LU
     >   + 3*( 2*S**4 - 2*N**4*T*U*(T**2 + U**2)
     >     + N**2*(T**2 + T*U + 2*U**2)*(2*T**2 + T*U + U**2) ) *L2S
     >   - 3*S*(- N**4*U*(2*T**2 - T*U + U**2)
     >       + (N**2*T - S)*(2*T**2 + 2*T*U + U**2) ) *L2T
     >   - 3*S*(- N**4*T*(2*U**2 - T*U + T**2 )
     >       + (N**2*U - S)*( 2*U**2 + 2*T*U + T**2)) *L2U
C
       PSITILDE6NS = V/(3.0D0*N**2*S**2*T*U) * PSITILDE6NS
C
      ELSEIF (PROCESS.EQ.4) THEN
C
       PSITILDE6NS =
     >  2*N**2*NFLAVOR*( - (66 + 27*PI**2)*S**2*T**2*U**2
     >           + 40*(S**6 + T**6 + U**6))
     >   + 2*N**3*(  6*(125 - 27*PI**2)*S**2*T**2*U**2
     >         - 4*(67 - 9*PI**2)*(S**6 + T**6 + U**6))
     >   + 6*N**2*(11*N - 2*NFLAVOR)*(S**2 + T**2 + U**2)**3 *LMU
     >   + 6*N**2*T*U*(S**2 + T**2 + U**2)
     >       *(  NFLAVOR *(5*T**2 + 2*T*U + 5*U**2)
     >         - 2*N*(7*T**2 - 8*T*U + 7*U**2) ) *LS
C
       PSITILDE6NS = PSITILDE6NS +
     >   + 6*N**2*S*U*(S**2 + T**2 + U**2)
     >       *(  NFLAVOR *(5*S**2 + 2*S*U + 5*U**2)
     >         - 2*N*(7*S**2 - 8*S*U + 7*U**2) ) *LT
     >   + 6*N**2*S*T*(S**2 + T**2 + U**2)
     >       *(  NFLAVOR* (5*S**2 + 2*S*T + 5*T**2)
     >         - 2*N*(7*S**2 - 8*S*T + 7*T**2) ) *LU
     >   - 36*N**2*U**2*( NFLAVOR*S*T*(S**2 + T**2)
     >             + 2*N*(2*S**4 + 2*S**3*T + 3*S**2*T**2 + 2*S*T**3
     >             + 2*T**4))*LS*LT
C
       PSITILDE6NS = PSITILDE6NS +
     >   - 36*N**2*S**2*(  NFLAVOR*T*U*(T**2 + U**2)
     >             + 2*N*(2*T**4 + 2*T**3*U + 3*T**2*U**2 + 2*T*U**3
     >             + 2*U**4))*LT*LU
     >   - 36*N**2*T**2*(  NFLAVOR*S*U*(S**2 + U**2)
     >             + 2*N*(2*S**4 + 2*S**3*U + 3*S**2*U**2 + 2*S*U**3
     >             + 2*U**4))*LS*LU
     >   + 18*N**2*S**2*T*U*(4*N*(T**2 + U**2)
     >             - NFLAVOR*(T**2 + 3*T*U + U**2))*L2S
     >   + 18*N**2*S*T**2*U*(4*N*(S**2 + U**2)
     >             - NFLAVOR*(S**2 + 3*S*U + U**2))*L2T
     >   + 18*N**2*S*T*U**2*(4*N*(S**2 + T**2)
     >             - NFLAVOR*(S**2 + 3*S*T + T**2))*L2U
C
       PSITILDE6NS = V/(9.0D0*S**2*T**2*U**2) * PSITILDE6NS
C
      ELSE
       PSITILDE6NS = 0.0D0
      ENDIF
      RETURN
      END
c
c End of the 2 --> 2 routines
c
