c This file must be always used with jetpdf in the jet package.
c It is derived from hvqpdfpho_nuc, and it includes the routines relevant to
c the PDFs in nuclei. Any given parton density in a proton embedded in a
c nucleus is either returned as in the case of free protons, or it is 
c expressed as follows
c    f_A = R_A * f_p
c where f_p is the parton density in the free proton. Upon calling MLMPDF_NUC,
c which is analogous to MLMPDF, one gets f_A. PDFPAR_NUC has to be called
c before MLMPDF_NUC, to initialize the relevant parameters
c
C-------------------------------------------------------------------------
      SUBROUTINE PRNTSF_RAT
C     prints details of the nuclear structure function sets
C-------------------------------------------------------------------------
      WRITE(*,100)
     #  ' Set     Authors       f_A or R_A=f_A/f_p'
     # ,'   1     ---- Returns R_A=1, for consistency checks ----'
     # ,'   2     EKS98                R_A'
 100  FORMAT(1X,A,100(/,1X,A))
      END


      SUBROUTINE PDFPAR_RAT(J,IH,IRATF,XLAM,SCHE,IRET)
c J is the PDF identification number
c IH is the particle type
c IRATF=0 if the corresponding set returns f_A
c      =1 if the corresponding set returns R_A in terms of valence and sea
c      =2 if the corresponding set returns R_A in terms of q and qbar
c XLAM is Lambda_QCD relevant to the set. If R_A is returned, it is set to -1
c SCHE is the scheme relevant to the set. If R_A is returned, it is set to **
c IRET=0 if the call has been successful, 1 otherwise
      PARAMETER (NPDF=2)
      IMPLICIT REAL * 8 (A-H,O-Z)
      CHARACTER * 2 SCHE,SCH(NPDF)
      DIMENSION XLA(NPDF)
      INTEGER IRATF0(NPDF)
      DATA IRATF0/
c fake set
     # 2,
c EKS98
     # 1/
      DATA XLA/
c fake set
     # -1.e0,
c EKS98
     # -1.e0/
      DATA SCH/
c fake set
     # '**',
c EKS98
     # '**'/
c
      IRET=0                             
      IF(J.LT.1.OR.J.GT.NPDF) THEN
        WRITE(*,*) ' NUCLEAR SET ',J,' NOT EXISTING'
        IRET=1
        RETURN
      ENDIF
      IF(ABS(IH).GT.2)THEN
         WRITE(*,*) ' PARTICLE TYPE ',IH,' MEANINGLESS'
         IRET=1
         RETURN
      ENDIF
      IF(
     # J.NE.1
C It is not the dummy set
     # .AND. J.NE.2 )
C It is not EKS98
     # THEN
         WRITE(*,*) ' PDF SET ',J,' NOT AVAILABLE'
         IRET=1
         RETURN
      ENDIF
      IRATF=IRATF0(J)
      XLAM=XLA(J)
      SCHE=SCH(J)
      END


      SUBROUTINE PDFPAR_NUC(JA,JP,IH,IRATF,XLAM,SCHE,IRET)
c Calls PDFPAR_RAT and PDFPAR. See the comments in PDFPAR_RAT
c for an explanation of the entries
      IMPLICIT REAL * 8 (A-H,O-Z)
      CHARACTER * 2 SCHE,SCHE0
c
      IRET=0
      CALL PDFPAR_RAT(JA,IH,IRATF0,XLAM0,SCHE0,IRET0)
      IF(IRATF.EQ.0)CALL PDFPAR(JP,IH,XLAM0,SCHE0,IRET0)
      IRATF=IRATF0
      XLAM=XLAM0
      SCHE=SCHE0
      END


C--------------------------------------------------
C- STRUCTURE FUNCTION MAIN PROGRAM
C--------------------------------------------------
      SUBROUTINE MLMPDF_NUC(IRATF,ANUC,NDNSA,NDNSP,IH,Q2,X,FX,NF)
c Analogous to MLMPDF, but relevant to nuclear pdfs. If IRATF=0,
c the nuclear PDFs are given as f_A, and a call to MLMPDF is not
c necessary. If IRAT=1 or 2, the nuclear PDFs are given as R_A * f_p, 
c and a call to MLMPDF is therefore necessary. When IRATF=1 or 2, NDNSA 
c is the id number of the R_A set, and NDNSP is the id number of the f_p 
c set. When IRATF=0, NDNSA is the id number of f_A, and NDNSP is unused
c When IRATF=1(2), R_A is returned in terms of valence and sea (q and qbar)
      REAL FX(-NF:NF),FXP(-NF:NF),RA(-NF:NF)
C Fix to prevent undefined math operations for x=1.
C Assumes that all structure functions vanish for x=1.
      IF(1-X.EQ.0) THEN
         DO J=-NF,NF
            FX(J) = 0.E0
         ENDDO
         RETURN
      ENDIF
C
      IH0=IH
      IF(IH.EQ.0) IH0=1
      IF(ABS(IH).EQ.2) IH0=ISIGN(1,IH)
      IF(NDNSA.EQ.1) THEN
         IF(IRATF.NE.2)THEN
            WRITE(*,*)'INCONSISTENT PARAMETERS IN MLMPDF_NUC',
     #                NDNSA,IRATF
            STOP
         ENDIF
         DO J =-NF,NF
            RA(J) = 1.E0
         ENDDO
      ELSEIF(NDNSA.LE.2) THEN
C--EKS98
         ISET=NDNSA-1
         CALL EKS98_INT(ANUC,ISET,IH0,Q2,X,RA,NF)
      ELSE
         WRITE(*,*)'NUCLEAR SET NOT IMPLEMENTED: NDNSA=',NDNSA
         STOP
      ENDIF
      IF(IRATF.EQ.0)THEN
         DO J=-NF,NF
            FX(J) = RA(J) 
         ENDDO
      ELSE
         CALL MLMPDF(NDNSP,IH0,Q2,X,FXP,NF)
         IF(IRATF.EQ.1)THEN
            FX(0) = FXP(0) * RA(0) 
            FX(-IH0) = FXP(-IH0) * RA(-IH0) 
            FX(-2*IH0) = FXP(-2*IH0) * RA(-2*IH0) 
            FX(IH0) = (FXP(IH0)-FXP(-IH0))*(RA(IH0)-RA(-IH0))+FX(-IH0)
            FX(2*IH0) = (FXP(2*IH0)-FXP(-2*IH0))*
     #                  (RA(2*IH0)-RA(-2*IH0))+FX(-2*IH0)
            DO J=3,NF
               FX(J) = FXP(J) * RA(J) 
            ENDDO
            DO J=-NF,-3
               FX(J) = FX(-J)
            ENDDO
         ELSE
            DO J=-NF,NF
               FX(J) = FXP(J) * RA(J) 
            ENDDO
         ENDIF
      ENDIF
      IF(IH.EQ.0) THEN
        FX(1)  = 0.5 * ( FX(1)+FX(2) )
        FX(-1) = 0.5 * ( FX(-1)+FX(-2) )
        FX(2)  = FX(1)
        FX(-2) = FX(-1)
      ELSEIF(ABS(IH).EQ.2) THEN
        T=FX(1)
        FX(1)=FX(2)
        FX(2)=T
        T=FX(-1)
        FX(-1)=FX(-2)
        FX(-2)=T
      ENDIF
      END


      SUBROUTINE EKS98_INT(ANUC,ISET,IH,Q2,X,RA,NF)
c Interface to eks98 routine
      REAL RA(-NF:NF)
      REAL*8 DX,DQ,DANUC,UPV,DOV,USEA,DSEA,STR,CHR,BOT,TOP,GLU
      REAL*8 XMIN,XMAX,QSQMIN,QSQMAX,QSQ
      REAL*8 IXMIN,IXMAX,IQSQMIN,IQSQMAX
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-6,1.D0,2.25D0,1.D4/
      DATA INI/0/
      IF(INI.GT.0) GO TO 1
        ILXMIN=0                    
        ILXMAX=0
        ILQSQMIN=0
        ILQSQMAX=0
        INI=1
1     CONTINUE
      IH0=IH
      IF(ABS(IH).NE.1) THEN
         WRITE(*,*)'ERROR IN EKS98_INT: IH=',IH
         STOP
      ENDIF
      Q=SQRT(Q2)
      DQ=DBLE(Q)
      DX=DBLE(X)
      DANUC=DBLE(ANUC)
      IF(DX.LT.XMIN) THEN
        IXMIN=IXMIN+1.
        IF(LOG10(IXMIN).GT.ILXMIN) THEN
      WRITE(*,*)' X < XMIN IN NUC STR. FUNCTIONS MORE THAN 10**',
     +  ILXMIN,' TIMES'                          
          ILXMIN=ILXMIN+1
        ENDIF
      ENDIF
      IF(DX.GT.XMAX) THEN
        IXMAX=IXMAX+1.
        IF(LOG10(IXMAX).GT.ILXMAX) THEN
      WRITE(*,*)' X > XMAX IN NUC STR. FUNCTIONS MORE THAN 10**',
     +  ILXMAX,' TIMES'
          ILXMAX=ILXMAX+1
        ENDIF
      ENDIF
      QSQ=DQ**2
      IF(QSQ.LT.QSQMIN) THEN
        IQSQMIN=IQSQMIN+1.
        IF(LOG10(IQSQMIN).GT.ILQSQMIN) THEN
      WRITE(*,*)'Q**2 < MIN Q**2 IN NUC STR. FUNCTIONS MORE THAN 10**',
     +  ILQSQMIN,' TIMES'
          ILQSQMIN=ILQSQMIN+1
        ENDIF
      ENDIF
      IF(QSQ.GT.QSQMAX) THEN
        IQSQMAX=IQSQMAX+1.
        IF(LOG10(IQSQMAX).GT.ILQSQMAX) THEN
      WRITE(*,*)'Q**2 > MAX Q**2 IN NUC STR. FUNCTIONS MORE THAN 10**',
     +  ILQSQMAX,' TIMES'
          ILQSQMAX=ILQSQMAX+1
        ENDIF
      ENDIF
      IF(ISET.LE.1) THEN
         CALL EKS98(DX,DQ,DANUC,UPV,DOV,USEA,DSEA,STR,CHR,BOT,TOP,GLU)
         RA(0)=SNGL(GLU)
         RA(-IH0)=SNGL(USEA)
         RA(-2*IH0)=SNGL(DSEA)
         RA(IH0)  =SNGL(UPV+USEA)
         RA(2*IH0)=SNGL(DOV+DSEA)
         IF(NF.GE.3) RA(3)=SNGL(STR)
         IF(NF.GE.4) RA(4)=SNGL(CHR)
         IF(NF.GE.5) RA(5)=SNGL(BOT)
         IF(NF.eq.6) RA(6)=0
      ENDIF
      DO I=3,NF
        RA(-I)=RA(I)
      ENDDO
      END


C***************************************************************************
C
C		 	eks98.f
C
C An interface for calculating the SCALE DEPENDENT NUCLEAR RATIOS
C 		R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 
C where f_A is the distribution of parton flavour f in a proton of a 
C nucleus A, and f_p is the corresponding parton distribution in the 
C free proton.
C  
C When you are using this interface, please REFER TO:
C K.J. Eskola, V.J. Kolhinen and C.A. Salgado, 
C "The scale dependent nuclear effects in parton distributions for 
C practical applications", Eur. Phys. J. C9 (1999) 61,
C JYFL-8/98, US-FT/14-98, hep-ph/9807297.
C
C The detailed formulation of our approach is given in 
C K.J. Eskola, V.J. Kolhinen and P.V. Ruuskanen,
C "Scale evolution of nuclear parton distributions"
C Nucl. Phys. B535 (1998) 351, CERN-TH/97-345, JYFL-2/98, hep-ph/9802350,
C so please refer also to this paper.
C
C The ratios R_f^A are to a good approximation independent of the choice
C of the parton distribution set for the free proton, so the absolute
C distributions of parton flavour f in a proton of a nucleus A can be 
C obtained simply by:
C f_A(x,Q) = R_f^A(x,Q) * f_p(x,Q), 
C where f_p is from any modern (lowest order) set of parton distributions.
C The corresponding distributions in a neutron of the nucleus can be 
C obtained through the isospin symmetry (=an approximation for non-isoscalar 
C nuclei)
C
C Questions & comments to:
C   	salgado@fpaxp1.usc.es
C   	vesa.kolhinen@phys.jyu.fi
C   	kari.eskola@phys.jyu.fi
C 
C August 4, 1998 / April 10, 2000 (new references added)
C-------------------------------------------------------------------
C
C INSTRUCTIONS:
C
C call eks98(x,Q,A,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
C
C Returns the nuclear corrections R_f^A(x,Q) in double precision for f=
C	
C	u_valence: ruv
C	d_valence: rdv (=ruv)
C	u_sea: ru
C	d_sea: rd (=ru)
C	s: rs
C	c: rc
C	b: rb
C	t: rt (always set to 1)
C	glue: rg
C
C For x, Q (Q is in GeV) and atomic number A. 
C x, Q and A are in DOUBLE PRECISION
C 
C No initialization is needed.
C
C This program needs data files par0.all and parxQA.all.
C They must be located in current working directory.
C
C This parametrization should only be applied at
C 1e-6 < x < almost 1,     1.5 < Q< 100 GeV
C Warning: No warning is given if the above kinematic region
C           in x&Q is exceeded.
C If A<=2, the function returns 1.
C
C
C
C-------------------------------------------------------------------

      subroutine eks98(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
      implicit double precision (a-h,o-z)
      dimension pqq(3)
      dimension r(5)
      common/eks983/qq0(3),x1,x116,aa1,aa8
      common/eks984/ptm(10)
      common/eks985/p0,p1,p2
      common/eks981/pa(3,10,3,8)
      common/eks982/pk0(3,180,8)
      common/eks986/indx,indpi,ikpt
      common/eks987/nm(5), kpt(5)
      common/eks988/ v1, v2
      data readFR/0/
      if (readFR.ne.1) then
         call eksinit
         readFR=1
      endif
      rt=1.d0
      if ((a.le.2.d0).or.(1.d0.le.x)) then
         ruv=1.d0
         rdv=1.d0
         ru=1.d0
         rd=1.d0
         rs=1.d0
         rc=1.d0
         rb=1.d0
         rg=1.d0
         return
      endif
      v1=dlog(a/aa1)
      v2=dlog(x/x1)
      if (x.le.0.6d0) then
         vtem=dlog10(0.6d0)+6.d0
         indx=1+dint((dlog10(x)+6.d0)/vtem*
     #        149.d0+1.d-7)
         if (indx.lt.2) indx=2
         else 
            indx=150+dint(75.d0*(x-0.6d0)+1.d-7)
      endif
      do 33 inum=1,5
         aqq=q*q
         nmin=nm(inum)
         ikpt=kpt(inum)
         if (aqq.lt.qq0(nmin)) aqq=qq0(nmin)
         qlq=dlog(aqq/qq0(nmin))
         q1=eksar0(a)
         indx=indx+1
         q2=eksar0(a)
         indx=indx-1
         if (x.le.0.6d0) then
            xxp=10.d0**(-6.d0+(indx-1)*vtem/149.d0)
            xxu=10.d0**(-6.d0+(indx)*vtem/149.d0)
            else
               xxp=0.6d0+(indx-150)*0.4d0/30.d0
               xxu=0.6d0+(indx+1-150)*0.4d0/30.d0
         endif
         r0=q1+(q2-q1)/(xxp-xxu)*(xxp-x)
         do 32 jk=1,3
            indpi=jk
            do 31 kl=1,10
               p0=pa(1,kl,indpi,ikpt)
               p1=pa(2,kl,indpi,ikpt)
               p2=pa(3,kl,indpi,ikpt)
               ptm(kl)=eksara(a)
31          continue
         pqq(jk)=eksarp(x)
32       continue
         r(inum)=r0+pqq(1)*qlq+pqq(2)*qlq*qlq+pqq(3)*dsqrt(qlq)
         if (inum.eq.2) then
            do 34 jk=1,3
               indpi=jk
               do 35 kl=1,10
                  p0=pa(1,kl,indpi,5)
                  p1=pa(2,kl,indpi,5)
                  p2=pa(3,kl,indpi,5)
                  ptm(kl)=eksara(a)
35             continue
            pqq(jk)=eksarp(x)
34          continue
            rs=r0+pqq(1)*qlq+pqq(2)*qlq*qlq+pqq(3)*dsqrt(qlq)
         endif
33    continue
      ruv=r(1)
      rdv=ruv
      ru=r(2)
      rd=ru
      rc=r(3)
      rb=r(4)
      rg=r(5)
      return
      end

      function eksar0(aa)
      implicit double precision (a-h,o-z)
      common/eks986/indx,indpi,kpt
      common/eks982/pk0(3,180,8)
      common/eks988/ v1, v2
      z=v1
      eksar0=pk0(1,indx,kpt)+pk0(2,indx,kpt)*z
     #+pk0(3,indx,kpt)*z*z
      return
      end

      function eksara(aa)
      implicit double precision (a-h,o-z)
      common/eks985/yy1,p1,p2
      common/eks988/ v1, v2
      z=v1
      eksara=yy1+p1*z+p2*z*z
      return
      end

      function eksarp(x)
      implicit double precision (a-h,o-z)
      common/eks986/indx,indpi,kpt
      common/eks983/qq0(3),x1,x116,aa1,aa8
      common/eks984/yy1,p1,p2,p3,p4,p5,p6,p7,p8,p9
      common/eks988/ v1, v2
      z=v2
      xx=x-x1
      if (x.le.x116) then
         eksarp=yy1+p1*z+p2*z**2+p3*xx+p4*xx**2
         else
            z1=dlog(x116/x1)
            xx1=x116-x1
            ff0=yy1+p1*z1+p2*z1**2+p3*xx1+p4*xx1**2
            xx16=x-x116
            z16=dlog(x/x116)
            qexp=19.d0+(indpi-1)*(indpi-2)*8.d0/2.d0
            eksarp=ff0+p5*xx16+p6*xx16**2+p7*xx16**3
     #             +p8*xx16**qexp+p9*z16
      endif
11    return
      end

      subroutine eksinit
      implicit double precision (a-h,o-z)
      common/eks981/pa(3,10,3,8)
      common/eks982/pk0(3,180,8)
      common/eks983/qq0(3),x1,x116,aa1,aa8
      common/eks987/nm(5), kpt(5)
      data qq0 /2.25d0, 2.54958d0, 21.3474d0/, 
     #     x1 /1.d-6/, x116 /.263553d-01/,
     #     aa1 /4.d0/, aa8 /208.d0/  
      data nm/1, 1, 2, 3, 1/,
     #     kpt/1, 3, 6, 7, 8/
      data readFR2/0/
      if (readFR2.ne.1) then
         open(11,file='parxQA.all',status='UNKNOWN')
         do 30 i=1,8
            do 40 j=1,3
               do 50 k=1,10
50             read(11,*),pa(1,k,j,i),pa(2,k,j,i),pa(3,k,j,i)
40       continue
30       continue
         close(11)
         open(11,file='par0.all',status='UNKNOWN')
         do 10 i=1,8
            do 20 j=1,180
20             read(11,137),pk0(1,j,i),pk0(2,j,i),pk0(3,j,i)
10       continue
         close(11)
         readFR2=1
      endif
      return
137   format(3e15.8)
      end
      


