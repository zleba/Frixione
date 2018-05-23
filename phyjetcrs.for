c**************************************************************************
c
c Classification of the partonic subprocesses for photon-hadron 
c collisions: every process contributing to the physical cross section
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
c         1              0 --> 1p2g2q
c         2              0 --> 1p2q2Q
c         3              0 --> 1p4q
c
c
c       jinst         initial state
c
c         1                p g
c         2                p q
c
c
c We have therefore the following physical processes
c
c    (jproc,jinst)    physical process
c
c        (1,1)        p,g --> q,qb,g
c        (1,2)        p,q --> q,g,g
c
c        (2,2)        p,q --> q,Q,Qb
c
c        (3,2)        p,q --> q,q,qb
c
c 
c The leading order processes 0 --> 4 can be formally obtained
c from the 0 --> 5 ones by eliminating one gluon. We will mantain 
c the same classification scheme used for the 0 --> 5 processes. 
c Jproc=1 will mean the process 0 --> 1p1g2q, and so on.  Furthermore,
c we need also 0 --> 4 processes which do not involve a photon, since
c they factorize in the collinear limit of the 0 --> 5 processes.
c We adopt the following conventions
c
c       jproc         unphysical process
c
c         1              0 --> 1p1g2q
c         4              0 --> 2g2q
c         5              0 --> 2q2Q
c         6              0 --> 4q
c
c
c       jinst         initial state
c
c         1                p g
c         2                p q
c         3                q g
c         4                q qb
c         5                q q
c         6                q Qb
c         7                q Q
c
c The following physical processes will therefore eventually occur
c
c    (jproc,jinst)    physical process
c
c        (1,1)        p,g --> q,qb
c        (1,2)        p,q --> q,g
c
c        (4,3)        q,g --> q,g
c        (4,4)        q,qb --> gg
c
c        (5,4)        q,qb --> Q,Qb
c        (5,6)        q,Qb --> q,Qb
c        (5,7)        q,Q --> q,Q
c
c        (6,4)        q,qb --> q,qb
c        (6,5)        q,q --> q,q
c
c**************************************************************************
c
c
c Begin of the 2 --> 3 routines
c
      subroutine xmatel_five_part
     # (s,xpp,xinv,xii,yi,xij,yj,jproc,jinst,ni,nj,p5e1q,p5e2q,p5e1e2)
c This subroutine, given the partonic CM energy (s) and the dot products of
c the four-momenta of the partons (xpp), returns the value of the 2 --> 3 
c matrix element (physical configuration) for the process specified 
c by (jproc,jinst)
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (tiny_s=1.d-8)
      parameter (tiny_c=1.d-8)
      common/i_type_part/i_type_part
      common/i_ap_in/i_ap_in1,i_ap_in2
      common/i_ap_out/i_ap_out
      common/jcoll_prc/j_prc_1_coll,j_prc_2_coll
      common/jcoll_ins/j_ins_1_coll,j_ins_2_coll
      common/jcoll_out_prc/j_prc_coll_out
      dimension xpp(1:5,1:5),ypp(1:5,1:5),xinv(1:5)
      integer i_type_part(1:3,1:2,3:5)
      integer i_ap_in1(1:3,1:2,3:5),i_ap_in2(1:3,1:2,3:5)
      integer i_ap_out(1:3,1:2,3:5,3:5)
      integer j_prc_1_coll(1:3,1:2,3:5),j_prc_2_coll(1:3,1:2,3:5)
      integer j_ins_1_coll(1:3,1:2,3:5),j_ins_2_coll(1:3,1:2,3:5)
      integer j_prc_coll_out(1:3,1:2,3:4,4:5)
c
      if(xii.lt.tiny_s)then
c Parton number ni is soft: use the analytical formula
c for the soft limit of the matrix element
        if(i_type_part(jproc,jinst,ni).ne.0)then
c The soft parton is not a gluon: the limit is finite, and the factor
c xii**2 in front of the matrix element gives therefore 0
          p5e1q=0.d0
          p5e2q=0.d0
          p5e1e2=0.d0
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
          p5e1q=tmp/(8*pi**2)
          p5e2q=0.d0
          p5e1e2=0.d0
        endif
      elseif(abs(1-yi).lt.tiny_c)then
c Parton number ni is collinear to parton 1: use the analytical formula
c for the collinear limit of the matrix element
        jproc0=j_prc_1_coll(jproc,jinst,ni)
        jinst0=j_ins_1_coll(jproc,jinst,ni)
        if(jproc0.eq.-1)then
c The splitting is not allowed: the limit is zero
          p5e1q=0.d0
          p5e2q=0.d0
          p5e1e2=0.d0
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
          p5e1q=(4/s)*(1+yi)*xfct1*xfct2
          p5e2q=0.d0
          p5e1e2=0.d0
          if( jproc.eq.2.and.jinst.eq.2 .and.
     #        (ni.eq.4.or.ni.eq.5) )then
            aaa=p5e1q
            p5e1q=0.d0
            p5e2q=aaa
          endif
        endif
      elseif(abs(1+yi).lt.tiny_c)then
c Parton number ni is collinear to parton 2: use the analytical formula
c for the collinear limit of the matrix element
        jproc0=j_prc_2_coll(jproc,jinst,ni)
        jinst0=j_ins_2_coll(jproc,jinst,ni)
        if(jproc0.eq.-1)then
c The splitting is not allowed: the limit is zero
          p5e1q=0.d0
          p5e2q=0.d0
          p5e1e2=0.d0
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
          p5e1q=(4/s)*(1-yi)*xfct1*xfct2
          p5e2q=0.d0
          p5e1e2=0.d0
          if( jproc.eq.2.and.jinst.eq.2 )then
            aaa=p5e1q
            p5e1q=0.d0
            p5e2q=aaa
          endif
        endif
      elseif(abs(1-yj).lt.tiny_c)then
c Parton number ni is collinear to parton number nj: use the analytical 
c formula for the collinear limit of the matrix element
        n1=min(ni,nj)
        n2=max(ni,nj)
        jproc0=j_prc_coll_out(jproc,jinst,n1,n2)
        if(jproc0.eq.-1)then
c The splitting is not allowed: the limit is zero
          p5e1q=0.d0
          p5e2q=0.d0
          p5e1e2=0.d0
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
          p5e1q=(4/s)*(1/x_ap)*xfct1*xfct2
          p5e2q=0.d0
          p5e1e2=0.d0
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
          call p5_pqqbgg(ypp,tmpe1q,tmpe2q,tmpe1e2)
        elseif(jproc.eq.2)then
          call p5_pqqbrrb(ypp,tmpe1q,tmpe2q,tmpe1e2)
        elseif(jproc.eq.3)then
          call p5_pqqbqqb(ypp,tmpe1q,tmpe2q,tmpe1e2)
        else
          write(6,*)'Error in xmatel_five_part: unknown process number' 
          stop
        endif
        xfact=xfact/x_average(jinst)/(2*s)
        p5e1q=xfact*tmpe1q
        p5e2q=xfact*tmpe2q
        p5e1e2=xfact*tmpe1e2
      endif
      return
      end


      subroutine xmatel_coll(s,xpp,xi_ev,xi,y,jproc,jinst,
     #          jproc0,jinst0,ni,p5e1q,p5e2q,p5e1q_sc,p5e2q_sc)
c This function, given a collinear configuration (xpp), returns
c the corresponding 2 --> 2 matrix elements times a factor depending
c upon the Altarelli-Parisi kernels and the K functions (which allow
c for a change of the factorization scheme), matching the definition
c of eq.(5.7). Notice that there is no e1e2 (interference) term,
c since such contribution is always zero in the collinear limit.
c p5e1q_sc and p5e2q_sc contain the contribution K^{(reg)} (regular
c part of the K function) when xi#0, and K^{(d)} (the part of the
c K function proportional to delta(xi)) when xi=0.
      implicit real * 8 (a-h,o-z)
      logical lperm
      character * 2 schhad,schpho,scheme
      parameter (pi=3.14159265358979312D0)
      parameter (one=1.d0)
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/cutpar/xicut,yincut,youtcut
      common/i_ap_in/i_ap_in1,i_ap_in2
      common/schhad/schhad
      common/schpho/schpho
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      integer i_ap_in1(1:3,1:2,3:5),i_ap_in2(1:3,1:2,3:5)
c
      lperm=.false.
c Momentum fraction entering the Altarelli-Parisi kernels
      x_ap=1-xi
c Reduced energy for the incoming partons
      s_red=s*x_ap
      if(y.eq.1.d0)then
c Parton number ni is collinear to parton number 1
        icode=i_ap_in1(jproc,jinst,ni)
        call x_xp_collin_yp(s_red,xpp,ypp,ni,icode,1)
        if( jproc.eq.2.and.jinst.eq.2 .and.
     #      (ni.eq.4.or.ni.eq.5) )lperm=.true.
        scheme=schpho
      elseif(y.eq.-1.d0)then
c Parton number ni is collinear to parton number 2
        icode=i_ap_in2(jproc,jinst,ni)
        call x_xp_collin_yp(s_red,xpp,ypp,ni,icode,-1)
        if( jproc.eq.2.and.jinst.eq.2 )lperm=.true.
        scheme=schhad
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
      p5e1q=( xfct1*(xlgsdi_mu+2*log(xi_ev))-xfct2
     #       -xfct3p-xfct3l*log(xi_ev) )*xfct4
      p5e2q=0.d0
      p5e1q_sc=-xfct5*xfct4
      p5e2q_sc=0.d0
      if(lperm)then
        aaa=p5e1q
        p5e1q=0.d0
        p5e2q=aaa
        aaa=p5e1q_sc
        p5e1q_sc=0.d0
        p5e2q_sc=aaa
      endif
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
c This function returns the quantity SIGN*OMEGA(GAMMA)*OMEGA(A), where
c the GAMMA and A are equivalently given by JINST, using our classification 
c for initial states. The quantity SIGN is -1 if JINST=2 (<==> pq initial 
c state), 1 otherwise, to take into account the sign coming from the 
c crossing of a single fermionic line
      implicit real*8 (a-h,o-z)
      parameter (xnc=3.d0)
      parameter (vda=8.d0)
c
      if(jinst.eq.1)then
        tmp=2*(2*vda)
        xsgn=1.d0
      elseif(jinst.eq.2)then
        tmp=2*(2*xnc)
        xsgn=-1.d0
      elseif(jinst.eq.3)then
        tmp=(2*xnc)*(2*vda)  
        xsgn=-1.d0
      else
        tmp=(2*xnc)**2
        xsgn=1.d0
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
c           5            qp
c
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
c
      if(index.eq.1)then
        ap_kern=2*vca*(x+(1-x)**2/x+x*(1-x)**2)
      elseif(index.eq.2)then
        ap_kern=vtf*(1-x)*(x**2+(1-x)**2)
      elseif(index.eq.3)then
        ap_kern=vcf*(1-x)*(1+(1-x)**2)/x
      elseif(index.eq.4)then
        ap_kern=vcf*(1+x**2)
      elseif(index.eq.5)then
        ap_kern=xnc*(1-x)*(x**2+(1-x)**2)
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
      parameter (xnc=3.d0)
c
      if(index.eq.1)then
        apprime_kern=0.d0
      elseif(index.eq.2)then
        apprime_kern=-2*vtf*x*(1-x)**2
      elseif(index.eq.3)then
        apprime_kern=-vcf*(1-x)*x
      elseif(index.eq.4)then
        apprime_kern=-vcf*(1-x)**2
      elseif(index.eq.5)then
        apprime_kern=-2*xnc*x*(1-x)**2
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
      elseif(index.eq.5)then
        xkdelta=0.d0
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
      elseif(index.eq.5)then
        xkplus=0.d0
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
      elseif(index.eq.5)then
        xklog=0.d0
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
      elseif(index.eq.5)then
        xkreg=xnc*( (x**2+(1-x)**2)*log((1-x)/x)+8*x*(1-x)-1 )
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
      integer i_phys(1:5,1:3,1:2)
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
      integer i0_phys(1:4,1:6,1:7)
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
c Formulae for the cross sections: from
c
c  P. Aurenche et al., Nucl. Phys. B286(87)553
c
      subroutine p5_pqqbgg(xpp,p5e1q,p5e2q,p5e1e2)
c This function returns the 0 --> pqqbgg matrix element. The spin and color
c average factors and the electric charge of the quarks are NOT included.
c The particle labeling is as follows
c
c   0 --> p(1) + q(2) + qb(3) + g(4) + g(5)
c
      implicit real * 8 (a-z)
      parameter (vcf=4.d0/3.d0)
      parameter (xnc=3.d0)
      dimension xpp(1:5,1:5)
c
      a1=xpp(1,3)
      a2=xpp(3,4)
      a3=xpp(3,5)
      a4=xpp(2,3)
      a5=xpp(1,2)
      a6=xpp(2,4)
      a7=xpp(2,5)
      a8=xpp(1,5)
      a9=xpp(4,5)
      a10=xpp(1,4)
      tmp=( 2*(vcf-xnc/2.d0)*a4+xnc*(a2*a7+a3*a6)/a9)
     #   *( (a1**2+a5**2)/(a2*a3*a6*a7)
     #     +(a2**2+a6**2)/(a1*a3*a5*a7)
     #     +(a3**2+a7**2)/(a1*a2*a5*a6) )
      p5e1q=4*vcf*xnc*tmp
      p5e2q=0.d0
      p5e1e2=0.d0
      return 
      end 


      subroutine p5_pqqbqqb(xpp,p5e1q,p5e2q,p5e1e2)
c This function returns the 0 --> pqqbqqb matrix element. The spin and color
c average factors and the electric charge of the quarks are NOT included. 
c The particle labeling is as follows
c
c   0 --> p(1) + q(2) + qb(3) + q(4) + qb(5)
c
      implicit real * 8 (a-z)
      parameter (vcf=4.d0/3.d0)
      parameter (xnc=3.d0)
      dimension xpp(1:5,1:5)
c
      a1=xpp(3,5)
      a2=xpp(2,5)
      a3=xpp(1,5)
      a4=xpp(4,5)
      a5=xpp(3,4)
      a6=xpp(2,4)
      a7=xpp(1,4)
      a8=xpp(1,3)
      a9=xpp(1,2)
      a10=xpp(2,3)
      esse=-a1/(a3*a8)+a2/(a3*a9)+a4/(a3*a7)
     #     +a5/(a7*a8)-a6/(a7*a9)+a10/(a8*a9)
      tmp=2*esse*(a1**2+a2**2+a5**2+a6**2)/(a4*a10)
     #   +2*esse*(a1**2+a4**2+a6**2+a10**2)/(a2*a5)
     #   -2*esse/xnc*(a1**2+a6**2)*(a1*a6-a2*a5-a4*a10)
     #              /(a2*a4*a5*a10)
      p5e1q=2*vcf*xnc*tmp
      p5e2q=0.d0
      p5e1e2=0.d0
      return
      end


      subroutine p5_pqqbrrb(xpp,p5e1q,p5e2q,p5e1e2)
c This function returns the 0 --> pqqbrrb matrix element. The spin and color
c average factors and the electric charge of the quarks are NOT included. 
c The particle labeling is as follows
c
c   0 --> p(1) + q(2) + qb(3) + Q(4) + Qb(5)
c
      implicit real * 8 (a-z)
      parameter (vcf=4.d0/3.d0)
      parameter (xnc=3.d0)
      dimension xpp(1:5,1:5)
c
      a1=xpp(3,5)
      a2=xpp(2,5)
      a3=xpp(1,5)
      a4=xpp(4,5)
      a5=xpp(3,4)
      a6=xpp(2,4)
      a7=xpp(1,4)
      a8=xpp(1,3)
      a9=xpp(1,2)
      a10=xpp(2,3)
      xfact=2*(a1**2+a2**2+a5**2+a6**2)/(a4*a10)
      tmpe1q=a10/(a8*a9)
      tmpe2q=a4/(a3*a7)
      tmpe1e2=-a1/(a3*a8)+a2/(a3*a9)+a5/(a7*a8)-a6/(a7*a9)
      p5e1q=2*vcf*xnc*xfact*tmpe1q
      p5e2q=2*vcf*xnc*xfact*tmpe2q
      p5e1e2=2*vcf*xnc*xfact*tmpe1e2
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
      dimension xpp(1:5,1:5),xca(0:2),xgmm(0:2),xgmmprime(0:2)
      dimension xe(1:4),xpnpm(1:3,2:4)
      integer i_type_part(1:3,1:2,3:5),i_type_in(1:2,1:2)
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
      integer imap(1:6)
c Map present classification scheme for unphysical processes to
c the Kunszt-Soper scheme (1=2q2Q, 2=4q, 3=2g2q, 4=4g)
      data imap/0,0,0,3,1,2/
c
c Get the dot products to calculate the relevant matrix element
      call xpp_to_ypp_0(xpp,ypp,jproc,jinst)
      if(jproc.eq.1)then
        tmp=p4_pqqbg(ypp)/x_average(jinst)
      elseif(jproc.ge.4)then
        tmp=psitilde4(imap(jproc),ypp)/x_average(jinst)
      else
        write(6,*)'Error in xmatel_four_part'
        stop
      endif
      xmatel_four_part=tmp/(2*s)
      return
      end


      function x_color_link(s,xpp,jproc,jinst,n,m)
c This function gives the color-linked Born squared amplitudes;
c it is relevant for the soft limit of the 2 --> 3 processes. The
c factor 8 pi**2 inserted is consistent with the normalization
c of eq.(3.3)
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vca=3.d0)
      common/i_xcol/i_xcol
      dimension xpp(1:5,1:5),ypp(1:4,1:4)
      integer i_xcol(1:4,1:3,1:2)
c
c Get the dot products to calculate the relevant matrix element
      call xpp_to_ypp_0(xpp,ypp,jproc,jinst)
      if(jproc.ne.1)then
        tmp=0.d0
      else
        n1=min(i_xcol(n,jproc,jinst),i_xcol(m,jproc,jinst))
        m1=max(i_xcol(n,jproc,jinst),i_xcol(m,jproc,jinst))
        if(n1.eq.1)then
          tmp=0.d0
        elseif(n1.eq.2.and.m1.eq.3)then
          tmp=(2*vcf-vca)*p4_pqqbg(ypp)
        elseif(n1.eq.2.and.m1.eq.4)then
          tmp=vca*p4_pqqbg(ypp)
        elseif(n1.eq.3.and.m1.eq.4)then
          tmp=vca*p4_pqqbg(ypp)
        else
          write(*,*)'Error in x_color_link'
          stop
        endif
      endif
      x_color_link=8*pi**2*tmp/x_average(jinst)/(2*s)
      return
      end


      function x_loop_four(s,xpp,jproc,jinst)
c This function returns the M_{NS}^{(3,1)} part of the virtual
c contribution, as defined in eq. (3.2)
      implicit real * 8 (a-h,o-z)
      dimension xpp(1:5,1:5),ypp(1:4,1:4)
c
c Get the dot products to calculate the relevant matrix element
      call xpp_to_ypp_0(xpp,ypp,jproc,jinst)
      if(jproc.ne.1)then
        tmp=0.d0
      else
        tmp=p4_loop_pqqbg(ypp)/x_average(jinst)
      endif
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
c Formulae for the cross sections: from
c
c   D. Boedeker, Z. Phys. C59(93)501 
c
      function p4_loop_pqqbg(xpp)
c This function returns the 0 --> pqqbg matrix element (virtual corrections). 
c The spin and color average factors and the electric charge of the quarks 
c are NOT included. The particle labeling is as follows
c
c  0 --> p(1) + q(2) + qb(3) + g(4) 
c
      implicit real * 8 (a-z)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vca=3.d0)
      integer nl
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/nl/nl
      dimension xpp(1:4,1:4)
c
      s=2*xpp(1,4)
      t=2*xpp(1,2)
      u=2*xpp(1,3)
      ls = dlog(dabs(s/xmues2))
      lt = dlog(dabs(t/xmues2))
      lu = dlog(dabs(u/xmues2))
      lmu = dlog(dabs(xmuf2/xmues2))
      l2s = ls**2 - pi**2 * theta(s.gt.0.0d0)
      l2t = lt**2 - pi**2 * theta(t.gt.0.0d0)
      l2u = lu**2 - pi**2 * theta(u.gt.0.0d0)
c
      tmp=(2*vcf-vca)*( 2*s**2*l2s+(t**2+s**2)*l2t+(u**2+s**2)*l2u
     #                 -2*(t**2+s**2)*ls*lt-2*(u**2+s**2)*ls*lu
     #                 -4*t*u*ls+u*(3*u+2*t)*lt+t*(3*t+2*u)*lu
     #                 +(2*s**2+u**2+t**2)*pi**2-7*(t**2+u**2) )
     #   +vca*( 3*(t**2*lu+u**2*lt)+(pi**2-7-2*lt*lu)*(t**2+u**2) )
      tmp=tmp/(2*u*t)+(11*vca-2.d0*nl)/6.d0*lmu*(t/u+u/t)
      p4_loop_pqqbg=32.d0*tmp
      return
      end


      function p4_pqqbg(xpp)
c This function returns the 0 --> pqqbg matrix element (Born). The spin and 
c color average factors and the electric charge of the quarks are NOT included.
c The particle labeling is as follows
c
c  0 --> p(1) + q(2) + qb(3) + g(4) 
c
      implicit real * 8 (a-z)
      dimension xpp(1:4,1:4)
c
      s=2*xpp(1,4)
      t=2*xpp(1,2)
      u=2*xpp(1,3)
      tmp=t/u+u/t
      p4_pqqbg=32.d0*tmp
      return
      end


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
