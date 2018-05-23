       PROGRAM PHYJETDIFF
c This program calculates fully differential distributions at 
c next-to-leading order in QCD for infrared-safe quantities in photon-hadron
c collisions, according to the description given in:
c
c S. Frixione, Z. Kunszt and A. Signer, Nucl. Phys. B467(96)399, 
c                                       hep-ph/9512328.
c S. Frixione, Nucl. Phys. B507(97)295, hep-ph/9706545.
c S. Frixione and G. Ridolfi, Nucl. Phys. B507(97)315, hep-ph/9707345.
c
      implicit real*8 (a-h,o-z)
      character * 2 schhad,schpho
      character * 25 fname
      character * 80 runstr
      character * 25 pref,prefn
      character * 7 newver
      logical isave,fl_lo,fl_nlo
      external sig0wb,sig2pvwb,sig2prwb,sig1awb,sig1bwb,inijet
      external sig0mon,sig2pvmon,sig2prmon,sig1amon,sig1bmon
      parameter (pi=3.14159265358979312D0)
c use newver='NEW' for vaxes, 'UNKNOWN' in other machines
c This variable is assigned in the subroutine sysdep.
      common/newver/newver
c----------------------------------------------------------
c Common blocks for the fixed variables
c
c     xsc_min2 = minimum transverse energy squared
c     xlam = lambda_QCD (5 flavours)
c     xmuf2 = factorization scale squared
c     xmur2 = renormalization scale squared
c     xmues2 = Ellis-Sexton scale squared
c     xmuww2 = Weizsaecker-Williams scale squared
c     zg = strong coupling = sqrt(4*pi*alfas)
c     ze2 = electron electric charge squared
c
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
c Scale factors: 
c
c     sclf2=xmuf2/(reference scale)**2
c     sclr2=xmur2/(reference scale)**2
c     scles2=xmues2/(reference scale)**2
c
      common/scalef/sclf2,sclr2,scles2
c Number of flavours
      common/nl/nl
c ndns2 = pdf type
c ih2 = hadron type 
      common/pdf/ih2,ndns2
c schhad,schpho = 'DI' for deep inelastic,  'MS' for msbar scheme
      common/schhad/schhad
      common/schpho/schpho
c sh = (p_had1 + p_had2)^2
      common/shadr/sh
c P-function algorithm
      common/i_type_sfun/i_type_sfun
c merging definition (for the P-function algorithm)
      common/i_type_merge/i_type_merge
c f and D**2 parameters (f is unused. If the cone algorithm is
c used, D**2=R is set by definition)
      common/ffact/ffact,djet2
c----------------------------------------------------------
c Integration control parameters
      common/samp/nsamp
      common/cutpar/xicut,yincut,youtcut
c----------------------------------------------------------
c Run identification string
      common/runstr/runstr
c----------------------------------------------------------
c Common blocks for sum over processes and initial states
      common/iinst/iinst
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/jphysproc/jproc,jinst
      common/irange/imax_inst,imin_inst
      common/i0range/imax0_inst,imin0_inst
c----------------------------------------------------------
c Arrays for multiple list of energy, scalefactors, minimum E_T
      dimension ecmlst(100),sclstf(100),xetminlst(100)
c----------------------------------------------------------
c Arrays for the sum over initial states
      integer imax_inst(1:3),imin_inst(1:3)
      integer imax0_inst(1:3),imin0_inst(1:3)
C------------------------------------------------------------------------
C                             START                                     -
C------------------------------------------------------------------------
c iinput=1 ==> inquire for P-functions, xi_cut,....
      iinput=0
c   i_type_merge=1 ==> weighted sum of eta and phi
      i_type_merge=1
c Set system dependent parameters
      call sysdep
c----- vegas prints nothing
c      call nopr(0)
c-----
c Open the file collecting all the input parameter. This file is meant 
c to be converted in a command file in a subsequent run
      open(unit=11,file='phjetlog',status=newver)
c Enter identification string for this run. It will be written
c in the save files and on the default output device 
      write(*,*)
     # 'Enter id. string (eg. ''test'') for this run (< 81 characters)'
      read(*,*) runstr
      write(11,*) ''''//runstr(1:istrl(runstr))//''''
c
      write(*,*)' '
      write(*,'(1x,a)') 
     #'Enter 1 if you want restart files and topdrawer file'
      read (*,*) irest
      write(11,'(1x,i1,32x,a)') irest,'! 1 for restart and .top files'
      if(irest.eq.1) then
         isave = .true.
      else
         isave = .false.
      endif
c 
      write(*,*)' '
      write(*,*)
     # 'Enter prefix for name of generated files (< 11 chars.)'
      read (*,*) pref
      write(11,*) ''''//pref(1:istrl(pref))//'''','  ! prefix for files'
c----------------------------------------------------------
c Parameters of the run
c
      alphaem=1/137.d0
      ze2=4*pi*alphaem
c
c Enter energy, scale factors and minimum transverse energy. Multiple 
c inputs (up to 100 values) are accepted by the program
      write(*,*)' '
      write(*,*)
     # 'Enter E_cm, scf, E_T(min) sequence, where:'
      write(*,*)
     # 'E_cm=CM energy of the colliding particles (in GeV),'
      write(*,*)
     # 'scf=mu_fac/mu_0=mu_ren/mu_0,'
      write(*,*)
     # 'E_T(min)=minimum transverse energy (in GeV).'
      write(*,*)
     # 'mu_0 is the reference scale, set in the user file.'
      write(*,*)
     # 'Enter all the entries < 0 to end the sequence'
      jecm = 0
      read(*,*) xecm,xscf,xetmin
      write(11,'(3(1x,d10.4),1x,a)') xecm,xscf,xetmin
     #,'! energy, scalefactor, E_T(min)'
      dowhile(jecm.lt.100.and.xecm.gt.0.and.xscf.gt.0.and.xetmin.gt.0)
         jecm=jecm+1
         ecmlst(jecm)=xecm
         sclstf(jecm)=xscf
         xetminlst(jecm)=xetmin
         read(*,*) xecm,xscf,xetmin
         write(11,'(3(1x,d10.4))') xecm,xscf,xetmin
      enddo
      if(jecm.eq.100.and.xecm.gt.0.and.xscf.gt.0.and.xetmin.gt.0)then
         write(*,*) 'no more than 100 values'
         stop
      endif
c Enter the P-function algorithm 
      if(iinput.eq.1)then
      write(*,*)' '
      write(*,*)'Enter the P-function algorithm: 1 for ES'
      write(*,*)'                                2 for cone'
      write(*,*)'                                3 for mod. cone'
      read(*,*)i_type_sfun
      write(11,'(1x,i2,31x,a)') i_type_sfun,
     # '! 1=ES, 2=cone, 3=mod. cone'
      else
      i_type_sfun=2
      endif
c Enter number of flavours, 1<=nl<=5
      write(*,*)' '
      write(*,*)'Enter number of flavours'
      read(*,*) nl
      write(11,'(1x,i2,31x,a)') nl,'! # of flavours'
c Set the parameters for the run
      call setpar
c Enter hadron colliding with the photon or the electron
      write(*,*)' '
      write(*,*)'Enter hadron type'
      write(*,*)'(0=(p+n)/2, 1=p, -1=pbar, 2=n, -2=nbar)'
      read(*,*) ih2
      write(11,'(1x,i2,31x,a)') ih2,'! Hadron type'
c
      write(*,*)' '
      call getset
     # ('Enter parton distribution set for the hadron ',ndns2,ih2)
      write(11,'(1x,i3,30x,a)') ndns2,'! PDF set for the hadron'
c
      call pdfpar(ndns2,ih2,xlam,schhad,iret)
c                              
      write(*,*)' '
      write(*,*)'Enter lambda_QCD_5, 0 for default'
      write(*,*)
     #  'Notice: this lambda MUST BE the one for FIVE FLAVOURS,'
      write(*,*)'regardless of the number of active flavours nl given'
      write(*,*)'before. Alpha_QCD will be correctly calculated with'
      write(*,*)'nl flavours'
      read(*,*) tmp
      write(11,'(1x,d10.4,23x,a)') tmp,'! Lambda_5, 0 for default'
c xlam is already fixed to default value: modify if entered value >0
      if(tmp.gt.1.d-15) xlam=tmp
      write(*,*) 'Lambda_5=',xlam,' GeV'
c assign hadron scheme, if not already given (this happens when the program
c is linked with a PDF library which is not included in this package)
      if(schhad.ne.'MS'.and.schhad.ne.'DI')then
        dowhile(schhad.ne.'MS'.and.schhad.ne.'DI')
          write(*,*)' '
          write(*,*)'Enter now scheme for hadron: ''DI'' or ''MS'''
          read(*,*) schhad
        enddo
        write(11,'(1x,a4,29x,a)')
     #   ''''//schhad//'''','! Scheme for hadron'
      endif
c assign photon scheme
      schpho='**'
      dowhile(schpho.ne.'MS'.and.schpho.ne.'DI')
         write(*,*)' '
         write(*,*)'Enter now scheme for photon: ''DI'' or ''MS'''
         write(*,*)
     # '(use DI when the PDFs in the photon in the hadronic'
          write(*,*)
     # 'component are defined in the DIS_gamma scheme, MS otherwise)'
         read(*,*) schpho
      enddo
      write(11,'(1x,a4,29x,a)')
     #   ''''//schpho//'''','! Scheme for photon'
      write(*,*) 'Scheme for hadron=',schhad
      write(*,*) 'Scheme for photon=',schpho
c---------------------------------------------------------------
c Inquire for monochromatic or broad band beam in
c photoproduction. 
      imonoch=-1
      dowhile(imonoch.ne.0.and.imonoch.ne.1)
        write(*,*)' '
        write(*,*)
     # 'Enter 0 to compute the cross sections with a monochromatic'
        write(*,*)
     # 'photon beam, enter 1 to use the photon distribution function'
        write(*,*)
     # 'in the user file'
        read(*,*) imonoch
      enddo
      write(11,'((1x,i1),32x,a)') imonoch,
     # '! 0 for monochromatic photon, 1 for broad band'
c-----------------------------------------------------------------
c Enter the parameters for the P-function algorithm (see the paper for details)
      if(iinput.eq.1)then
      if(i_type_sfun.eq.1)then
        write(*,*)' '
        write(*,*) 'Ellis-Soper algorithm'
        write(*,*) 'Enter the parameter D (D<2*pi/3)'
        read(5,*)djet
        djet2=djet**2
        ffact=0.d0
        write(11,'(1x,d10.4,23x,a)')djet,'! D'
        if(djet.ge.2*pi/3.d0)then
          write(*,*)'Please enter a smaller value for D'
          stop
        endif
      elseif(i_type_sfun.eq.2)then
        write(*,*)' '
        write(*,*) 'Cone algorithm'
        write(*,*) 'Enter the parameter R (R<pi/3)'
        read(5,*)djet
        djet2=djet
        ffact=0.d0
        write(11,'(1x,d10.4,23x,a)')djet,'! R'
        if(djet.ge.pi/3.d0)then
          write(*,*)'Please enter a smaller value for R'
          stop
        endif
      elseif(i_type_sfun.eq.3)then
        write(*,*)' '
        write(*,*) 'Modified cone algorithm'
        write(*,*) 'Enter the parameter R (R<pi/3)'
        read(5,*)djet
        djet2=djet
        ffact=0.d0
        write(11,'(1x,d10.4,23x,a)')djet,'! R'
        if(djet.ge.pi/3.d0)then
          write(*,*)'Please enter a smaller value for R'
          stop
        endif
      else
        write(*,*)' '
        write(*,*) 'Error: algorithm not implemented'
        stop 
      endif
      else
      djet2=1.d0
      ffact=0.d0
      endif
c-----------------------------------------------------------------
c tau generated according to a flat distribution in (1/tau)**nsamp
      nsamp = 1
c See the paper for details about xi_cut, y_I and y_O
      if(iinput.eq.1)then
      write(*,*)' '
      write(*,*) 'Enter integration control parameters'
      write(*,*) '  0<xi_cut<=1'
      write(*,*) '  0<y_I<=2'
      write(*,*) '  0<y_O<=2'
      read(5,*)xicut_0,yincut,youtcut
      write(11,'(3(1x,d8.2),7x,a)')xicut_0,yincut,youtcut
     #,'! xi_cut, y_I, y_O'
      else
      xicut_0=0.1d0
      yincut=0.1d0
      youtcut=0.1d0
      endif
c---------------------------------------------------------------
c Sum over the processes: user can select few of them
      write(*,*)' '
      write(*,*) 'Select the allowed partonic processes:'
      write(*,*) '              NLO             LO'
      write(*,*) 'enter 1 for 0 --> 1p2g2q    0 --> 1p1g2q'
      write(*,*) 'enter 2 for 0 --> 1p2q2Q                '
      write(*,*) 'enter 3 for 0 --> 1p4q                  '
      write(*,*) 'enter 0 to select all the processes'
      write(*,*) 'enter -1 to exclude 1p2g2q, and so on'
      read(*,*) iproc
      write(11,'(1x,i2,31x,a)') iproc
     #,'! 1=1p2g2q,2=1p2q2Q,3=1p4q,0=all'
      min_proc = 1
      max_proc = 3
      if(iproc.gt.0) then
c- do only process iproc
         min_proc=iproc
         max_proc=iproc
      endif
c---------------------------------------------------------------
c Sum over the initial states: user can select few of them
      write(*,*)' '
      write(*,*) 'Select the allowed initial states:'
      write(*,*) 'enter 1 for pg, 2 for pq'
      write(*,*) 'enter 0 to select all the initial states'
      write(*,*) 'enter -1 to exclude pg, -2 to exclude pq'
      read(*,*) iinst
      write(11,'(1x,i2,31x,a)') iinst
     #,'! 1=pg,2=pq,0=all'
      call x_inst_sum(iinst)
      call x_proc_sum(fl_lo,fl_nlo)
c---------------------------------------------------------------
c Enter Vegas parameters
      write(*,*)' '
      write(*,*) 'Enter number of iterations'
      call integrms1
      write(*,*) 'n1 and n2 for LO and NLO contribution'
      read(*,*) n0,n3
      write(11,'(2(1x,i4),24x,a)') n0,n3,
     #  '! # of iterations'
      call integrms2
      write(*,*) 'i1 and i2 for LO and NLO contribution'
      read(*,*) inew0,inew3
      write(11,'(2(1x,i2),28x,a)') inew0,inew3,
     #  '! 0 to exclude, 1 for new run, 2 to restart'
c
      write(*,*)' '
      write(*,*)
     # 'Enter number of calls for vegas ncl2,ncl3,<0 for defaults'
      read(*,*)ncl2,ncl3
      if(ncl2.lt.0)ncl2=80000
      if(ncl3.lt.0)ncl3=400000
      write(11,'(2(1x,i9),14x,a)')ncl2,ncl3,
     # '! # of calls for vegas'
c---- close logfile
      close(11)
c----------------------------------------------------------------
c  *********************  START INTEGRATION *********************
c----------------------------------------------------------------
      do jloop=1,jecm
c Main loop (over energies and scale factors)
        avtot = 0.d0
        dtot = 0.d0
        avtott0 = 0.d0
        dtott0 = 0.d0
        prefn = pref
        if(jecm.gt.1) call strnum(prefn,jloop)
        ecm = ecmlst(jloop)
        sh = ecm**2
        xsc_min=xetminlst(jloop)
        xsc_min2=xsc_min**2
        xicut=xicut_0
        if(xicut.gt.1.d0-xsc_min2/sh)then
          xicut=1.d0-xsc_min2/sh
        endif
c
        sclff = sclstf(jloop)
        sclfes = sclff
c Common block values for scale factor
        sclf2 = sclff**2
        scles2 = sclfes**2
c In a more refined version, the possibility will be given to choose
c different values for factorization and renormalization scale.
c In this version, the renormalization scale is set equal to the
c factorization one
        sclr2 = sclf2
c
        av0=0.d0
        av2pv=0.d0
        av2pr=0.d0
        av1a=0.d0
        av1b=0.d0
        d0=0.d0
        d2pv=0.d0
        d2pr=0.d0
        d1a=0.d0
        d1b=0.d0
c
c LO contribution
c
        if(fl_lo)then
          call strcat(prefn,'_0',fname)
          if(imonoch.eq.1)then
            call integrate
     #  (inijet,sig0wb,fname,n0,inew0,3,ncl2,av0,d0,chi2a,isave)
          else
            call integrate
     #  (inijet,sig0mon,fname,n0,inew0,2,ncl2,av0,d0,chi2a,isave)
          endif
          write(*,*) 'LO    :  ',av0,' +- ',d0
        endif
        avtot0 = av0
        dtot0 = d0
c
c NLO contribution
c
        if(fl_lo)then
c 2PV contribution
          call strcat(prefn,'_2pv',fname)
          if(imonoch.eq.1)then
            call integrate
     #  (inijet,sig2pvwb,fname,n3,inew3,3,ncl2,av2pv,d2pv,chi2a,isave)
          else
            call integrate
     #  (inijet,sig2pvmon,fname,n3,inew3,2,ncl2,av2pv,d2pv,chi2a,isave)
          endif
          write(*,*) '2PV   :  ',av2pv,' +- ',d2pv
        endif
        if(fl_nlo)then
c 2PR contribution
          call strcat(prefn,'_2pr',fname)
          if(imonoch.eq.1)then
            call integrate
     #  (inijet,sig2prwb,fname,n3,inew3,4,ncl2,av2pr,d2pr,chi2a,isave)
          else
            call integrate
     #  (inijet,sig2prmon,fname,n3,inew3,3,ncl2,av2pr,d2pr,chi2a,isave)
          endif
          write(*,*) '2PR   :  ',av2pr,' +- ',d2pr
c S(0) contribution
          call strcat(prefn,'_1a',fname)
          if(imonoch.eq.1)then
            call integrate
     #  (inijet,sig1awb,fname,n3,inew3,5,ncl3,av1a,d1a,chi2a,isave)
          else
            call integrate
     #  (inijet,sig1amon,fname,n3,inew3,4,ncl3,av1a,d1a,chi2a,isave)
          endif
          write(*,*) 'NLO(A):  ',av1a,' +- ',d1a
c S(1) contribution
          call strcat(prefn,'_1b',fname)
          if(imonoch.eq.1)then
            call integrate
     #  (inijet,sig1bwb,fname,n3,inew3,5,ncl3,av1b,d1b,chi2a,isave)
          else
            call integrate
     #  (inijet,sig1bmon,fname,n3,inew3,4,ncl3,av1b,d1b,chi2a,isave)
          endif
          write(*,*) 'NLO(B):  ',av1b,' +- ',d1b
        endif
c
        avtot = av0+av2pv+av2pr+av1a+av1b
        dtot = sqrt(d0**2+d2pv**2+d2pr**2+d1a**2+d1b**2)
c Write the .top file if the option was selected
        if(isave) then
          call inijet
          call mclear
          call strcat(prefn,'phjet.top',fname)
          open(unit=99,file=fname,status=newver)
          call topout
          close(99)
        endif
c Write the total integral on the default output device
        write(*,*) 'total'
        write(*,200)ih2,ndns2,nl,xlam,schpho
        write(*,250) 'tt','tt'
        write(*,300)
     #    ecmlst(jloop),sclstf(jloop),
     #    avtot,dtot,avtot0,dtot0
        write(*,'(4(1x,/))')
c End of the main loop
      enddo
200   format(' hadron=',i2,' strfun=',i3,
     # ' nl=',i2,'  lambda5=',d10.4,'  scheme(ph)=',a)
250   format(' ecm       mu/mass   ',a,
     # '         err      ',a,'0        err0')
300   format(2(1x,1pd9.3),4(1x,0pd10.4,1x,1pd8.2))
      end


      subroutine getset(str,ndns,ih)
      implicit real*8 (a-h,o-z)
      character * (*) str, scheme * 2
 2    write(*,*) str
      write(*,*)
     # '   (< 0 for a display of the features of the various sets)'
      read(*,*) ndns
      if(ndns.lt.0) then
        call prntsf
        go to 2
      endif
      call pdfpar(ndns,ih,xlam,scheme,iret)
      if(iret.ne.0) goto 2
      end


      subroutine toend(iunit)
      ios = 0    
      dowhile(ios.eq.0)
         read(unit=iunit,fmt='(1x)',iostat=ios)
      enddo                        
      end


      subroutine strfun(x1,x2,sf)
c This subroutines returns the luminosity for a given initial state
c Assuming one quark flavour to simplify the writing, the values
c returned are (a_ = qb_, the sum over the flavour is not indicated)
c
c   index=1    -->   sf(1,1) = p(x1)*g_h2(x2)
c                    sf(1,2) = 0
c                    sf(1,3) = 0
c
c   index=2    -->   sf(2,1) = p(x1)*(q_h2(x2)+a_h2(x2))*e_q^2
c                    sf(2,2) = 0
c                    sf(2,3) = 0
c
c   index=3    -->   sf(3,1) = p(x1)*(q_h2(x2)+a_h2(x2))*e_q^2
c                    sf(3,2) = p(x1)*(q_h2(x2)+a_h2(x2))*e_Q^2
c                    sf(3,3) = p(x1)*(q_h2(x2)-a_h2(x2))*e_q*e_Q
c
c
c x1 and x2 are the Bjorken variables. The value of index is equal
c to the value of jinst, except for the case (jproc,jinst)=(2,2),
c when, by defintion, index=3
c
      implicit real*8 (a-h,o-z)
      real*4 fh2x2(-5:5)
      real*8 zel(1:6),sumzel(1:6),sumzel2(1:6)
      parameter(pi=3.14159265358979312D0)
      dimension sf(1:3,1:3)
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/pdf/ih2,ndns2
      common/jphysproc/jproc,jinst
      common/nl/nl
      common/charges/zel,sumzel,sumzel2
c
      call mlmpdf(ndns2,ih2,sngl(xmuf2),sngl(x2),fh2x2,5)
c
      do i=1,3
        do j=1,3
          sf(i,j) = 0
        enddo
      enddo
c
      if(x1.gt.1.5) then
        px1 = 1.d0
      else
        px1 = phdistr(x1,xmuww2)
      endif
c
c index=1 ==> jinst=1
      sf(1,1) = px1*dble(fh2x2(0))*sumzel2(nl)
c index=2 ==> jinst=2,jproc#2
      do i=1,nl
        sf(2,1) = sf(2,1) + px1*dble(fh2x2(i)+fh2x2(-i))*zel(i)**2
      enddo
c index=3 ==> jinst=2,jproc=2
      do i=1,nl
        sf(3,1) = sf(3,1) + 
     #         px1*dble(fh2x2(i)+fh2x2(-i))*zel(i)**2*(nl-1)
        sf(3,2) = sf(3,2) + 
     #         px1*dble(fh2x2(i)+fh2x2(-i))*(sumzel2(nl)-zel(i)**2)
        sf(3,3) = sf(3,3) + 
     #         px1*dble(fh2x2(i)-fh2x2(-i))*zel(i)*(sumzel(nl)-zel(i))
      enddo
      return
      end
c
c Cross section routines
c
c
c Begin of leading-order contribution
c
      function sig0wb(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(3)
      sig0wb=sig0(xx,peso,0)
      return
      end


      function sig0mon(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(2)
      sig0mon=sig0(xx,peso,1)
      return
      end


      function sig0(xx,wght,iphot)
c
      implicit real*8 (a-h,o-z)
      logical xpass,s_2_fun
      parameter (pi=3.14159265358979312D0)
      parameter (phij=0.d0)
      parameter (n1=5)
      parameter (n2=3)
      parameter (zero=0.d0)
      parameter (one=1.d0)
      dimension xx(*)
      common/shadr/sh
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/x1x2/ycm,tau
      common/samp/nsamp
      common/caso/caso(5)
      common/i0range/imax0_inst,imin0_inst
      common/i0allwd/i0_allwd_inst
      common/i0stat/i0_stat_fact
      common/isfmap/i_sfmap
      common/imtt/imtt
      common/iinst/iinst
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/jphysproc/jproc,jinst
      common/icls/icls
      dimension sf(1:3,1:3),vvev(1:3)
      dimension xpp(1:5,1:5),xinv(1:5)
      integer imax0_inst(1:3),imin0_inst(1:3)
      integer i0_allwd_inst(1:2,1:3)
      integer i0_stat_fact(1:1,1:2)
      integer i_sfmap(1:3,1:2)
      integer imtt(3:4,3:5)
c
      icls=0
      xjac=1.d0
      roh = xsc_min2/sh
c
      ximax0 = roh**(-nsamp)
      ximin0 = 1
      tmp  = ximin0 + xx(1)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
c
      if(iphot.eq.0)then
        ymax= -log(tau)/2
        ymin=  log(tau)/2
        ycm = ymin + xx(2)*(ymax-ymin)
        xjac= xjac * (ymax-ymin)
        xxtmp=xx(3)
      else
        x1 = 1.d0
        x2 = tau
        ycm=log(1/x2)/2.d0
        xxtmp=xx(2)
      endif
c 
      ro=roh/tau
      yjlim=sqrt(1.d0-ro)
      call zzchvar(xxtmp,xyj,xjac,ro)
      yj=xyj*yjlim
      xjac=xjac*yjlim
c s is the partonic CM energy
      s=sh*tau
c
      xnorm=xjac/(16*pi)
c
      sig0sum=0.d0
      vvev(1)=0.d0
      vvev(2)=0.d0
      vvev(3)=0.d0
c Generate the kinematics for the event
      call invar_in(n1,n2,s,zero,one,zero,yj,phij,xpp,xinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
      xpass=s_2_fun(n1,-1)
      if(xpass)then
c Evaluate Bjorken variables
        if(iphot.eq.0) then
          x1=sqrt(tau) * exp(ycm)
          x2=tau/x1
        else
          x1=2.d0
          x2=tau
        endif
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
        zg2=zgmu2()
c Calculate parton luminosity
        call strfun(x1,x2,sf)
c Loop over selected processes
        do jproc=min_proc,max_proc
          if(jproc.ne.-iproc) then
c Loop over selected initial states
            do i=imin0_inst(jproc),imax0_inst(jproc)
              jinst=i0_allwd_inst(i,jproc)
              if(jinst.ne.-iinst)then
c Statistical factor for this process
                xstat_fact=i0_stat_fact(jproc,jinst)/2.d0
c Value of the matrix element for the physical process
                xmtproc=xmatel_four_part(s,xpp,jproc,jinst)
c Weight for the event
                www=xmtproc*xnorm*zg2*ze2*wght*xstat_fact
c Store the results
                index=i_sfmap(jproc,jinst)
                vvev(1)=vvev(1)+sf(index,1)*www
              endif
c End of the loop over initial states
            enddo
          endif
c End of the loop over processes
        enddo
      endif
c Output
      if(vvev(1).ne.0.d0)then
        if(iphot.eq.0) then
          x1p=sqrt(tau) * exp(ycm)
        else
          x1p=1.d0
        endif
        www=1.d0
        call outall(www,tot,vvev(1),vvev(2),vvev(3),x1p)
        sig0sum=tot/wght
      endif
      sig0=sig0sum
      return
      end
c
c End of leading-order contribution
c
c
c Begin of next-to-leading-order contribution
c
      function sig2pvwb(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(3)
      sig2pvwb=sig2pv(xx,peso,0)
      return
      end


      function sig2pvmon(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(2)
      sig2pvmon=sig2pv(xx,peso,1)
      return
      end


      function sig2pv(xx,wght,iphot)
c
      implicit real*8 (a-h,o-z)
      logical xpass,s_2_fun
      parameter (pi=3.14159265358979312D0)
      parameter (phii=0.d0)
      parameter (phij=0.d0)
      parameter (n1=5)
      parameter (n2=3)
      parameter (zero=0.d0)
      parameter (one=1.d0)
      dimension xx(*)
      common/shadr/sh
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/x1x2/ycm,tau
      common/samp/nsamp
      common/caso/caso(5)
      common/i0range/imax0_inst,imin0_inst
      common/i0allwd/i0_allwd_inst
      common/i0stat/i0_stat_fact
      common/isfmap/i_sfmap
      common/imtt/imtt
      common/iinst/iinst
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/jphysproc/jproc,jinst
      common/icls/icls
      dimension sf(1:3,1:3),vv2pv(1:3)
      dimension xpp(1:5,1:5),xinv(1:5)
      integer imax0_inst(1:3),imin0_inst(1:3)
      integer i0_allwd_inst(1:2,1:3)
      integer i0_stat_fact(1:1,1:2)
      integer i_sfmap(1:3,1:2)
      integer imtt(3:4,3:5)
c
      icls=0
      xjac=1.d0
      roh = xsc_min2/sh
c 
      ximax0 = roh**(-nsamp)
      ximin0 = 1
      tmp  = ximin0 + xx(1)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
c 
      if(iphot.eq.0)then
        ymax= -log(tau)/2
        ymin=  log(tau)/2
        ycm = ymin + xx(2)*(ymax-ymin)
        xjac= xjac * (ymax-ymin)
        xxtmp=xx(3)
      else
        x1 = 1.d0
        x2 = tau
        ycm=log(1/x2)/2.d0
        xxtmp=xx(2)
      endif
c 
      ro=roh/tau
      yjlim=sqrt(1.d0-ro)
      call zzchvar(xxtmp,xyj,xjac,ro)
      yj=xyj*yjlim
      xjac=xjac*yjlim
c s is the partonic CM energy
      s=sh*tau
c
      xnorm=xjac/(16*pi)
c The factor 1/(8*pi**2) is coming from alpha_s/(2*pi), eq.(5.5)
      xnorm=xnorm/(8*pi**2)
c
      sig2pvsum=0.d0
      vv2pv(1)=0.d0
      vv2pv(2)=0.d0
      vv2pv(3)=0.d0
c Generate the kinematics for the event
      call invar_in(n1,n2,s,zero,one,phii,yj,phij,xpp,xinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
      xpass=s_2_fun(n1,-1)
      if(xpass)then
c Evaluate Bjorken variables
        if(iphot.eq.0) then
          x1=sqrt(tau) * exp(ycm)
          x2=tau/x1
        else
          x1=2.d0
          x2=tau
        endif
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
        zg2=zgmu2()
c Calculate parton luminosity
        call strfun(x1,x2,sf)
c Loop over selected processes
        do jproc=min_proc,max_proc
          if(jproc.ne.-iproc) then
c Loop over selected initial states
            do i=imin0_inst(jproc),imax0_inst(jproc)
              jinst=i0_allwd_inst(i,jproc)
              if(jinst.ne.-iinst)then
c Statistical factor for this process
                xstat_fact=i0_stat_fact(jproc,jinst)/2.d0
c Value of the matrix element for the physical process
                xmtproc=xmatel_2pv_contr(s,xpp,jproc,jinst)
c Weight for the event
                www=xmtproc*xnorm*zg2**2*ze2*wght*xstat_fact
c Store the results
                index=i_sfmap(jproc,jinst)
                vv2pv(1)=vv2pv(1)+sf(index,1)*www
              endif
c End of the loop over initial states
            enddo
          endif
c End of the loop over processes
        enddo
      endif
c Output
      if(vv2pv(1).ne.0.d0)then
        if(iphot.eq.0) then
          x1p=sqrt(tau) * exp(ycm)
        else
          x1p=1.d0
        endif
        www=1.d0
        call outall(www,tot,vv2pv(1),vv2pv(2),vv2pv(3),x1p)
        sig2pvsum=tot/wght
      endif
      sig2pv=sig2pvsum
      return
      end


      function sig2prwb(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(4)
      sig2prwb=sig2pr(xx,peso,0)
      return
      end


      function sig2prmon(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(3)
      sig2prmon=sig2pr(xx,peso,1)
      return
      end


      function sig2pr(xx,wght,iphot)
c
      implicit real*8 (a-h,o-z)
      logical xpass,s_2_fun
      parameter (pi=3.14159265358979312D0)
      parameter (phii=0.d0)
      parameter (phij=0.d0)
      parameter (tiny=1.d-6)
      parameter (two=2.d0)
      parameter (n1=3)
      parameter (n2=4)
      dimension xx(*)
      common/shadr/sh
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/x1x2/ycm,tau
      common/samp/nsamp
      common/caso/caso(5)
      common/cutpar/xicut,yincut,youtcut
      common/irange/imax_inst,imin_inst
      common/iallwd/i_allwd_inst
      common/iskip/i_skip,i_skcnt
      common/istat/i_stat_fact
      common/jcoll_prc/j_prc_1_coll,j_prc_2_coll
      common/jcoll_ins/j_ins_1_coll,j_ins_2_coll
      common/isfmap/i_sfmap
      common/imtt/imtt
      common/iinst/iinst
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/jphysproc/jproc,jinst
      common/icls/icls
      dimension sf(1:3,1:3),vv2pr(1:3),ww2pr(1:3)
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      dimension xinv(1:5),yinv(1:5)
      integer imax_inst(1:3),imin_inst(1:3)
      integer i_allwd_inst(1:2,1:3)
      integer i_skip(1:3,1:2,3:5),i_skcnt(1:3,1:2,3:5)
      integer i_stat_fact(1:3,1:2)
      integer j_prc_1_coll(1:3,1:2,3:5),j_prc_2_coll(1:3,1:2,3:5)
      integer j_ins_1_coll(1:3,1:2,3:5),j_ins_2_coll(1:3,1:2,3:5)
      integer i_sfmap(1:3,1:2)
      integer imtt(3:4,3:5)
c
      icls=0
      xjac=1.d0
      roh = xsc_min2/sh
c 
      xi=(1.d0-roh)*xx(1)**2
      xjac=xjac*xx(1)*2*(1.d0-roh)
c 
      rxi=sqrt(1-xi)
      if(iphot.eq.0) then
        x1hat=1/rxi
        x2hat=1/rxi
      else
        x1hat = 1.d0
        x2hat = 1/(1-xi)
      endif
      rohxi=roh/(1-xi)
      ximax0 = rohxi**(-nsamp)
      ximin0 = (x1hat*x2hat)**(-nsamp)
      tmp  = ximin0 + xx(3)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
c 
      if(iphot.eq.0) then
        ymax= -log(tau)/2 + log(x1hat)
        ymin=  log(tau)/2 - log(x2hat)
        ycm = ymin + xx(4)*(ymax-ymin)
        xjac= xjac * (ymax-ymin)
      else
        ycm=log(1/tau)/2.d0
      endif
c 
      ro=roh/(tau*(1-xi))
      yjlim=sqrt(1.d0-ro)
      call zzchvar(xx(2),xyj,xjac,ro)
      yj=xyj*yjlim
      xjac=xjac*yjlim
c s is the partonic CM energy
      s=sh*tau
c
      xnorm=xjac/(16*pi)
c The factor 1/(8*pi**2) is coming from alpha_s/(2*pi), eq.(5.7)
      xnorm=xnorm/(8*pi**2)
c
      xsum_2prvv=0.d0
      xsum_2prww=0.d0
      vv2pr(1)=0.d0
      vv2pr(2)=0.d0
      vv2pr(3)=0.d0
      ww2pr(1)=0.d0
      ww2pr(2)=0.d0
      ww2pr(3)=0.d0
c Evaluate the normalization for non-soft configurations
      x_ns_fct=xnorm
c Evaluate the normalization for soft configurations
      x_s_fct=xnorm*(1-xi)
c
c Evaluate Bjorken variables
      if(iphot.eq.0) then
        x1=sqrt(tau) * exp(ycm)
        x2=tau/x1
      else
        x1=1.d0-1.d-8
        x2=tau
      endif
      rxi=sqrt(1-xi)
      sxi=s*(1-xi)
c
c Event (xi,y)=(xi,1)
c
      ytmp=1.d0
      if(iphot.eq.0) then
        x1t=x1/rxi
        x2t=x2*rxi
      else
        x1t=1.d0-1.d-8
        x2t=x2
      endif
      if(x1t.lt.1.and.x2t.lt.1)then
c Generate the kinematics
        call invar_in(n1,n2,s,xi,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
        xpass=s_2_fun(n1,-1)
        if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
          zg2=zgmu2()
c Calculate parton luminosity
          if(iphot.eq.0) then
            call strfun(x1t,x2t,sf)
          else
            call strfun(two,x2t,sf)
          endif
c Loop over selected processes
          do jproc=min_proc,max_proc
            if(jproc.ne.-iproc) then
c Loop over selected initial states
              do i=imin_inst(jproc),imax_inst(jproc)
                jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                do ni=3,5
                  jproc0_1=j_prc_1_coll(jproc,jinst,ni)
                  jinst0_1=j_ins_1_coll(jproc,jinst,ni)
c If jproc0_1=-1 the splitting is not allowed. Skip the contribution
                  if(jinst.ne.-iinst.and.jproc0_1.ne.-1
     #               .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                    xstat_fact=i_stat_fact(jproc,jinst)
     #                        *i_skcnt(jproc,jinst,ni)/6.d0
                    call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                    call xmatel_coll(s,xpp,xi,xi,ytmp,jproc,jinst,
     #     jproc0_1,jinst0_1,ni,xmte1q,xmte2q,xmte1q_sc,xmte2q_sc)
c Weight for the event
                    www0=x_ns_fct*zg2**2*ze2*wght*xstat_fact
                    www=www0/xi
c Store the results
                    index=i_sfmap(jproc,jinst)
                    vv2pr(1)=vv2pr(1)+sf(index,1)
     #                           *(www*xmte1q+www0*xmte1q_sc)
                    vv2pr(2)=vv2pr(2)+sf(index,2)
     #                           *(www*xmte2q+www0*xmte2q_sc)
                  endif
c End of the loop over parton number
                enddo
c End of the loop over initial states
              enddo
            endif
c End of the loop over processes
          enddo
        endif
      endif
c End of the event (xi,y)=(xi,1) contribution
c
c Event (xi,y)=(xi,-1)
c
      ytmp=-1.d0
      if(iphot.eq.0) then
        x1t=x1*rxi
        x2t=x2/rxi
      else
        x1t=1.d0-1.d-8
        x2t=x2
      endif
      if(x1t.lt.1.and.x2t.lt.1)then
c Generate the kinematics
        call invar_in(n1,n2,s,xi,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
        xpass=s_2_fun(n1,-1)
        if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
          zg2=zgmu2()
c Calculate parton luminosity
          if(iphot.eq.0) then
            call strfun(x1t,x2t,sf)
          else
            call strfun(two,x2t,sf)
          endif
c Loop over selected processes
          do jproc=min_proc,max_proc
            if(jproc.ne.-iproc) then
c Loop over selected initial states
              do i=imin_inst(jproc),imax_inst(jproc)
                jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                do ni=3,5
                  jproc0_2=j_prc_2_coll(jproc,jinst,ni)
                  jinst0_2=j_ins_2_coll(jproc,jinst,ni)
c If jproc0_2=-1 the splitting is not allowed. Skip the contribution
                  if(jinst.ne.-iinst.and.jproc0_2.ne.-1
     #               .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                    xstat_fact=i_stat_fact(jproc,jinst)
     #                        *i_skcnt(jproc,jinst,ni)/6.d0
                    call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                    call xmatel_coll(s,xpp,xi,xi,ytmp,jproc,jinst,
     #     jproc0_2,jinst0_2,ni,xmte1q,xmte2q,xmte1q_sc,xmte2q_sc)
c Weight for the event
                    www0=x_ns_fct*zg2**2*ze2*wght*xstat_fact
                    www=www0/xi
c Store the results
                    index=i_sfmap(jproc,jinst)
                    ww2pr(1)=ww2pr(1)+sf(index,1)
     #                           *(www*xmte1q+www0*xmte1q_sc)
                    ww2pr(2)=ww2pr(2)+sf(index,2)
     #                           *(www*xmte2q+www0*xmte2q_sc)
                  endif
c End of the loop over parton number
                enddo
c End of the loop over initial states
              enddo
            endif
c End of the loop over processes
          enddo
        endif
      endif
c End of the event (xi,y)=(xi,-1) contribution
c
c Soft counterevents
c
      if(xi.lt.xicut)then
        xtmp=0.d0
        if(iphot.eq.0) then
          x1t=x1*rxi
          x2t=x2*rxi
        else
          x1t=1.d0-1.d-8
          x2t=(1-xi)*x2
        endif
        if(x1t.lt.1.and.x2t.lt.1)then
c
c Counterevent (xi,y)=(0,1)
c
          ytmp=1.d0
c Generate the kinematics
          call invar_in(n1,n2,sxi,xtmp,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
          xpass=s_2_fun(n1,-1)
          if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
            zg2=zgmu2()
c Calculate parton luminosity
            if(iphot.eq.0) then
              call strfun(x1t,x2t,sf)
            else
              call strfun(two,x2t,sf)
            endif
c Loop over selected processes
            do jproc=min_proc,max_proc
              if(jproc.ne.-iproc) then
c Loop over selected initial states
                do i=imin_inst(jproc),imax_inst(jproc)
                  jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                  do ni=3,5
                    jproc0_1=j_prc_1_coll(jproc,jinst,ni)
                    jinst0_1=j_ins_1_coll(jproc,jinst,ni)
c If jproc0_1=-1 the splitting is not allowed. Skip the contribution
                    if(jinst.ne.-iinst.and.jproc0_1.ne.-1
     #                 .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                      xstat_fact=i_stat_fact(jproc,jinst)
     #                          *i_skcnt(jproc,jinst,ni)/6.d0
                      call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                      call xmatel_coll(sxi,xpp,xi,xtmp,ytmp,jproc,jinst,
     #       jproc0_1,jinst0_1,ni,xmte1q,xmte2q,xmte1q_sc,xmte2q_sc)
c Weight for the event
                      www0=x_s_fct*zg2**2*ze2*wght*xstat_fact
                      www=-www0/xi
c Store the results
                      index=i_sfmap(jproc,jinst)
                      ww2pr(1)=ww2pr(1)+sf(index,1)
     #                             *(www*xmte1q+www0*xmte1q_sc)
                      ww2pr(2)=ww2pr(2)+sf(index,2)
     #                             *(www*xmte2q+www0*xmte2q_sc)
                    endif
c End of the loop over parton number
                  enddo
c End of the loop over initial states
                enddo
              endif
c End of the loop over processes
            enddo
          endif
c End of the counterevent (xi,y)=(0,1) contribution
c
c Counterevent (xi,y)=(0,-1)
c
          ytmp=-1.d0
c Generate the kinematics
          call invar_in(n1,n2,sxi,xtmp,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
          xpass=s_2_fun(n1,-1)
          if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
            zg2=zgmu2()
c Calculate parton luminosity
            if(iphot.eq.0) then
              call strfun(x1t,x2t,sf)
            else
              call strfun(two,x2t,sf)
            endif
c Loop over selected processes
            do jproc=min_proc,max_proc
              if(jproc.ne.-iproc) then
c Loop over selected initial states
                do i=imin_inst(jproc),imax_inst(jproc)
                  jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                  do ni=3,5
                    jproc0_2=j_prc_2_coll(jproc,jinst,ni)
                    jinst0_2=j_ins_2_coll(jproc,jinst,ni)
c If jproc0_2=-1 the splitting is not allowed. Skip the contribution
                    if(jinst.ne.-iinst.and.jproc0_2.ne.-1
     #                .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                      xstat_fact=i_stat_fact(jproc,jinst)
     #                          *i_skcnt(jproc,jinst,ni)/6.d0
                      call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                      call xmatel_coll(sxi,xpp,xi,xtmp,ytmp,jproc,jinst,
     #       jproc0_2,jinst0_2,ni,xmte1q,xmte2q,xmte1q_sc,xmte2q_sc)
c Weight for the event
                      www0=x_s_fct*zg2**2*ze2*wght*xstat_fact
                      www=-www0/xi
c Store the results
                      index=i_sfmap(jproc,jinst)
                      ww2pr(1)=ww2pr(1)+sf(index,1)
     #                             *(www*xmte1q+www0*xmte1q_sc)
                      ww2pr(2)=ww2pr(2)+sf(index,2)
     #                             *(www*xmte2q+www0*xmte2q_sc)
                    endif
c End of the loop over parton number
                  enddo
c End of the loop over initial states
                enddo
              endif
c End of the loop over processes
            enddo
          endif
c End of the soft part
        endif
      endif
c 
c Output
      if(vv2pr(1).ne.0.d0.or.vv2pr(2).ne.0.d0 .or.
     #   ww2pr(1).ne.0.d0.or.ww2pr(2).ne.0.d0)then
        www=1.d0
        xtmp=0.d0
        ytmp=1.d0
        call invar_in(n1,n2,sxi,xtmp,ytmp,phii,yj,phij,ypp,yinv)
c
        if(vv2pr(1).ne.0.d0.or.vv2pr(2).ne.0.d0)then
          if(iphot.eq.0) then
            x1p=x1/rxi
            call outall(www,tot,vv2pr(1),vv2pr(2),vv2pr(3),x1p)
          else
            x1p=1.d0
            ycm0 = ycm
            ycm = ycm + log(1-xi)/2.d0 
            call outall(www,tot,vv2pr(1),vv2pr(2),vv2pr(3),x1p)
            ycm = ycm0
          endif
          xsum_2prvv=tot/wght
        endif
c
        if(ww2pr(1).ne.0.d0.or.ww2pr(2).ne.0.d0)then
          if(iphot.eq.0) then
            x1p=x1*rxi
            call outall(www,tot,ww2pr(1),ww2pr(2),ww2pr(3),x1p)
          else
            x1p=1.d0
            ycm0 = ycm
            ycm = ycm - log(1-xi)/2.d0 
            call outall(www,tot,ww2pr(1),ww2pr(2),ww2pr(3),x1p)
            ycm = ycm0
          endif
          xsum_2prww=tot/wght
        endif
      endif
c
      sig2pr=xsum_2prvv+xsum_2prww
      return
      end


      function sig1awb(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(5)
      sig1awb=sig1a(xx,peso,0)
      return
      end


      function sig1amon(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(4)
      sig1amon=sig1a(xx,peso,1)
      return
      end


      function sig1a(xx,wght,iphot)
c
      implicit real*8 (a-h,o-z)
      logical xpass,s0fun,s_2_fun,s0fun_soft
      parameter (pi=3.14159265358979312D0)
      parameter (phii=0.d0)
      parameter (tiny=1.d-6)
      parameter (two=2.d0)
      parameter (n1=3)
      parameter (n2=4)
      dimension xx(*)
      common/shadr/sh
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/x1x2/ycm,tau
      common/samp/nsamp
      common/caso/caso(5)
      common/cutpar/xicut,yincut,youtcut
      common/irange/imax_inst,imin_inst
      common/iallwd/i_allwd_inst
      common/iskip/i_skip,i_skcnt
      common/istat/i_stat_fact
      common/i0stat/i0_stat_fact
      common/jcoll_prc/j_prc_1_coll,j_prc_2_coll
      common/jcoll_ins/j_ins_1_coll,j_ins_2_coll
      common/isfmap/i_sfmap
      common/imtt/imtt
      common/iinst/iinst
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/jphysproc/jproc,jinst
      common/icls/icls
      dimension sf(1:3,1:3),vvev(1:3),vvcnt(1:3),wwcnt(1:3)
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      dimension xinv(1:5),yinv(1:5)
      integer imax_inst(1:3),imin_inst(1:3)
      integer i_allwd_inst(1:2,1:3)
      integer i_skip(1:3,1:2,3:5),i_skcnt(1:3,1:2,3:5)
      integer i_stat_fact(1:3,1:2)
      integer i0_stat_fact(1:1,1:2)
      integer j_prc_1_coll(1:3,1:2,3:5),j_prc_2_coll(1:3,1:2,3:5)
      integer j_ins_1_coll(1:3,1:2,3:5),j_ins_2_coll(1:3,1:2,3:5)
      integer i_sfmap(1:3,1:2)
      integer imtt(3:4,3:5)
c
      icls=0
      xjac=1.d0
      roh = xsc_min2/sh
      s2lim = xsc_min2/(3*sh)
c 
      xii=(1.d0-s2lim)*xx(1)**2
      xjac=xjac*xx(1)*2*(1.d0-s2lim)
c 
      zzz=1-2*xx(2)
      xjac=xjac*2
      ttt=tiny+(1-tiny)*zzz**2
      xjac=xjac*2*abs(zzz)
      if(zzz.gt.0) then
         th=ttt*pi/2
      else
         th=pi-ttt*pi/2
      endif
      xjac=xjac*pi/2
      yi=cos(th)
      xjac=xjac*sin(th)
c 
      omega=sqrt( (2-xii*(1+yi))/(2-xii*(1-yi)) )
      rxi=sqrt(1-xii)
      if(iphot.eq.0) then
        x1hat = max(1.d0,rxi/omega,1/(omega*rxi))
        x2hat = max(1.d0,omega/rxi,omega*rxi)
      else
        x1hat = 1.d0
        x2hat = 1/(1-xii)
      endif
      rohxi=max(roh,s2lim/(1-xii))
      ximax0 = rohxi**(-nsamp)
      ximin0 = (x1hat*x2hat)**(-nsamp)
      tmp  = ximin0 + xx(4)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
c 
      if(iphot.eq.0) then
        ymax= -log(tau)/2 + log(x1hat)
        ymin=  log(tau)/2 - log(x2hat)
        ycm = ymin + xx(5)*(ymax-ymin)
        xjac= xjac * (ymax-ymin)
      else
        ycm=log(1/tau)/2.d0
      endif
c 
      ro=min(1.d0-tiny,roh/(tau*(1-xii)))
      yjlim=sqrt(1.d0-ro)
      call jetchvar(xx(3),yj,yjlim,xjac,ro)
c caso(5) --> phij
      phij=2*pi*caso(5)
      xjac=xjac*2*pi
c s is the partonic CM energy
      s=sh*tau
c
      xnorm=xjac/(1024*pi**4)
c
      xsum_ev=0.d0
      xsum_cntvv=0.d0
      xsum_cntww=0.d0
      vvev(1)=0.d0
      vvev(2)=0.d0
      vvev(3)=0.d0
      vvcnt(1)=0.d0
      vvcnt(2)=0.d0
      vvcnt(3)=0.d0
      wwcnt(1)=0.d0
      wwcnt(2)=0.d0
      wwcnt(3)=0.d0
c Evaluate the normalization for non-soft configurations
      x_ns_fct=xnorm
c Evaluate the normalization for soft configurations
      x_s_fct=xnorm*(1-xii)
c
c Evaluate Bjorken variables
      if(iphot.eq.0) then
        x1=sqrt(tau) * exp(ycm)
        x2=tau/x1
      else
        x1=1.d0-1.d-8
        x2=tau
      endif
      rxi=sqrt(1-xii)
      sxi=s*(1-xii)
      omega=sqrt( (2-xii*(1+yi))/(2-xii*(1-yi)) )
c
c Event
c
      if(x1.lt.1.and.x2.lt.1)then
c Generate the kinematics
        call invar_in(n1,n2,s,xii,yi,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
        xpass=s0fun(n1)
        if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
          zg2=zgmu2()
c Calculate parton luminosity
          if(iphot.eq.0) then
            call strfun(x1,x2,sf)
          else
            call strfun(two,x2,sf)
          endif
c Factor appearing in eq.(4.37) and phase-space s
          xfact_ev=s/(2*xii)*(1.d0/(1-yi)+1.d0/(1+yi))
c Loop over selected processes
          do jproc=min_proc,max_proc
            if(jproc.ne.-iproc) then
c Loop over selected initial states
              do i=imin_inst(jproc),imax_inst(jproc)
                jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                do ni=3,5
                  if(jinst.ne.-iinst
     #              .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                    xstat_fact=i_stat_fact(jproc,jinst)
     #                        *i_skcnt(jproc,jinst,ni)/6.d0
                    call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                    call xmatel_five_part(s,xpp,xinv,xii,yi,
     #    2.d0,2.d0,jproc,jinst,ni,-1,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                    www=xfact_ev*x_ns_fct
     #                 *zg2**2*ze2*wght*xstat_fact
c Store the results
                    index=i_sfmap(jproc,jinst)
                    vvev(1)=vvev(1)+sf(index,1)*www*xmte1q
                    vvev(2)=vvev(2)+sf(index,2)*www*xmte2q
                    vvev(3)=vvev(3)+sf(index,3)*www*xmte1e2
                  endif
c End of the loop over parton number
                enddo
c End of the loop over initial states
              enddo
            endif
c End of the loop over processes
          enddo
        endif
c Output for the event
        if(vvev(1).ne.0.d0.or.vvev(2).ne.0.d0.or.vvev(3).ne.0.d0)then
          if(iphot.eq.0) then
            x1p=sqrt(tau) * exp(ycm)
          else
            x1p=1.d0
          endif
          www=1.d0
          call outall(www,tot,vvev(1),vvev(2),vvev(3),x1p)
          xsum_ev=tot/wght
        endif
      endif
c End of the event contribution
c
c Counterevent S(0) (xii,1)
c
      if(yi.gt.(1-yincut))then
        ytmp=1.d0
        if(iphot.eq.0) then
          x1t=x1*omega/rxi
          x2t=x2*rxi/omega
        else
          x1t=1.d0-1.d-8
          x2t=x2
        endif
        if(x1t.lt.1.and.x2t.lt.1)then
c Generate the kinematics
          call invar_in(n1,n2,s,xii,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
          xpass=s_2_fun(n1,-1)
          if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
            zg2=zgmu2()
c Calculate parton luminosity
            if(iphot.eq.0) then
              call strfun(x1t,x2t,sf)
            else
              call strfun(two,x2t,sf)
            endif
c Factor appearing in eq.(4.37) and phase-space s
            xfact_c1_s0=-s/(2*xii)*1.d0/(1-yi)
c Loop over selected processes
            do jproc=min_proc,max_proc
              if(jproc.ne.-iproc) then
c Loop over selected initial states
                do i=imin_inst(jproc),imax_inst(jproc)
                  jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                  do ni=3,5
                    jproc0_1=j_prc_1_coll(jproc,jinst,ni)
                    jinst0_1=j_ins_1_coll(jproc,jinst,ni)
c If jproc0_1=-1 the splitting is not allowed. Skip the contribution
                    if(jinst.ne.-iinst.and.jproc0_1.ne.-1
     #                 .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                      xstat_fact=i_stat_fact(jproc,jinst)
     #                          *i_skcnt(jproc,jinst,ni)/6.d0
                      call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                      call xmatel_five_part(s,xpp,xinv,xii,ytmp,
     #      2.d0,2.d0,jproc,jinst,ni,-1,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                      www=xfact_c1_s0*x_ns_fct
     #                   *zg2**2*ze2*wght*xstat_fact
c Store the results
                      index=i_sfmap(jproc,jinst)
                      vvcnt(1)=vvcnt(1)+sf(index,1)*www*xmte1q
                      vvcnt(2)=vvcnt(2)+sf(index,2)*www*xmte2q
                      vvcnt(3)=vvcnt(3)+sf(index,3)*www*xmte1e2
                    endif
c End of the loop over parton number
                  enddo
c End of the loop over initial states
                enddo
              endif
c End of the loop over processes
            enddo
          endif
        endif
      endif
c End of counterevent S(0) (xii,1)
c
c Counterevent S(0) (xii,-1)
c
      if(yi.lt.(-1+yincut))then
        ytmp=-1.d0
        if(iphot.eq.0) then
          x1t=x1*omega*rxi
          x2t=x2/(rxi*omega)
        else
          x1t=1.d0-1.d-8
          x2t=x2
        endif
        if(x1t.lt.1.and.x2t.lt.1)then
c Generate the kinematics
          call invar_in(n1,n2,s,xii,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
          xpass=s_2_fun(n1,-1)
          if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
            zg2=zgmu2()
c Calculate parton luminosity
            if(iphot.eq.0) then
              call strfun(x1t,x2t,sf)
            else
              call strfun(two,x2t,sf)
            endif
c Factor appearing in eq.(4.37) and phase-space s
            xfact_c2_s0=-s/(2*xii)*1.d0/(1+yi)
c Loop over selected processes
            do jproc=min_proc,max_proc
              if(jproc.ne.-iproc) then
c Loop over selected initial states
                do i=imin_inst(jproc),imax_inst(jproc)
                  jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                  do ni=3,5
                    jproc0_2=j_prc_2_coll(jproc,jinst,ni)
                    jinst0_2=j_ins_2_coll(jproc,jinst,ni)
c If jproc0_2=-1 the splitting is not allowed. Skip the contribution
                    if(jinst.ne.-iinst.and.jproc0_2.ne.-1
     #                 .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                      xstat_fact=i_stat_fact(jproc,jinst)
     #                          *i_skcnt(jproc,jinst,ni)/6.d0
                      call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                      call xmatel_five_part(s,xpp,xinv,xii,ytmp,
     #      2.d0,2.d0,jproc,jinst,ni,-1,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                      www=xfact_c2_s0*x_ns_fct
     #                   *zg2**2*ze2*wght*xstat_fact
c Store the results
                      index=i_sfmap(jproc,jinst)
                      wwcnt(1)=wwcnt(1)+sf(index,1)*www*xmte1q
                      wwcnt(2)=wwcnt(2)+sf(index,2)*www*xmte2q
                      wwcnt(3)=wwcnt(3)+sf(index,3)*www*xmte1e2
                    endif
c End of the loop over parton number
                  enddo
c End of the loop over initial states
                enddo
              endif
c End of the loop over processes
            enddo
          endif
        endif
      endif
c End of counterevent S(0) (xii,-1)
c
c Soft counterevents
c
      if(xii.lt.xicut)then
        xtmp=0.d0
        if(iphot.eq.0) then
          x1t=x1*omega*rxi
          x2t=x2*rxi/omega
        else
          x1t=1.d0-1.d-8
          x2t=(1-xii)*x2
        endif
        if(x1t.lt.1.and.x2t.lt.1)then
c
c S(0) (0,yi) counterevent
c
c Generate the kinematics
          call invar_in(n1,n2,sxi,xtmp,yi,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
          xpass=s_2_fun(n1,-1) .and. s0fun_soft(n1)
          if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
            zg2=zgmu2()
c Calculate parton luminosity
            if(iphot.eq.0) then
              call strfun(x1t,x2t,sf)
            else
              call strfun(two,x2t,sf)
            endif
c Factor appearing in eq.(4.37) and phase-space s
            xfact_s_s0=-sxi/(2*xii)*(1.d0/(1-yi)+1.d0/(1+yi))
c Loop over selected processes
            do jproc=min_proc,max_proc
              if(jproc.ne.-iproc) then
c Loop over selected initial states
                do i=imin_inst(jproc),imax_inst(jproc)
                  jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                  do ni=3,5
                    if(jinst.ne.-iinst
     #                .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                      xstat_fact=i_stat_fact(jproc,jinst)
     #                          *i_skcnt(jproc,jinst,ni)/6.d0
                      call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                      call xmatel_five_part(sxi,xpp,xinv,xtmp,yi,
     #      2.d0,2.d0,jproc,jinst,ni,-1,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                      www=xfact_s_s0*x_s_fct
     #                   *zg2**2*ze2*wght*xstat_fact
c Store the results
                      index=i_sfmap(jproc,jinst)
                      wwcnt(1)=wwcnt(1)+sf(index,1)*www*xmte1q
                      wwcnt(2)=wwcnt(2)+sf(index,2)*www*xmte2q
                      wwcnt(3)=wwcnt(3)+sf(index,3)*www*xmte1e2
                    endif
c End of the loop over parton number
                  enddo
c End of the loop over initial states
                enddo
              endif
c End of the loop over processes
            enddo
          endif
c End of counterevent S(0) (0,yi) 
c
          if(yi.gt.(1-yincut))then
c
c Counterevent S(0) (0,1)
c
            ytmp=1.d0
c Generate the kinematics
            call invar_in(n1,n2,sxi,xtmp,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
            xpass=s_2_fun(n1,-1)
            if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
              zg2=zgmu2()
c Calculate parton luminosity
              if(iphot.eq.0) then
                call strfun(x1t,x2t,sf)
              else
                call strfun(two,x2t,sf)
              endif
c Factor appearing in eq.(4.37) and phase-space s
              xfact_sc1_s0=sxi/(2*xii)*1.d0/(1-yi)
c Loop over selected processes
              do jproc=min_proc,max_proc
                if(jproc.ne.-iproc) then
c Loop over selected initial states
                  do i=imin_inst(jproc),imax_inst(jproc)
                    jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                    do ni=3,5
                      jproc0_1=j_prc_1_coll(jproc,jinst,ni)
                      jinst0_1=j_ins_1_coll(jproc,jinst,ni)
c If jproc0_1=-1 the splitting is not allowed. Skip the contribution
                      if(jinst.ne.-iinst.and.jproc0_1.ne.-1
     #                   .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                        xstat_fact=i_stat_fact(jproc,jinst)
     #                            *i_skcnt(jproc,jinst,ni)/6.d0
                        call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                        call xmatel_five_part(sxi,xpp,xinv,xtmp,ytmp,
     #        2.d0,2.d0,jproc,jinst,ni,-1,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                        www=xfact_sc1_s0*x_s_fct
     #                     *zg2**2*ze2*wght*xstat_fact
c Store the results
                        index=i_sfmap(jproc,jinst)
                        wwcnt(1)=wwcnt(1)+sf(index,1)*www*xmte1q
                        wwcnt(2)=wwcnt(2)+sf(index,2)*www*xmte2q
                        wwcnt(3)=wwcnt(3)+sf(index,3)*www*xmte1e2
                      endif
c End of the loop over parton number
                    enddo
c End of the loop over initial states
                  enddo
                endif
c End of the loop over processes
              enddo
            endif
          endif
c End of counterevent S(0) (0,1)
c
          if(yi.lt.(-1+yincut))then
c
c Counterevent S(0) (0,-1)
c
            ytmp=-1.d0
c Generate the kinematics
            call invar_in(n1,n2,sxi,xtmp,ytmp,phii,yj,phij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
            xpass=s_2_fun(n1,-1)
            if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
              zg2=zgmu2()
c Calculate parton luminosity
              if(iphot.eq.0) then
                call strfun(x1t,x2t,sf)
              else
                call strfun(two,x2t,sf)
              endif
c Factor appearing in eq.(4.37) and phase-space s
              xfact_sc2_s0=sxi/(2*xii)*1.d0/(1+yi)
c Loop over selected processes
              do jproc=min_proc,max_proc
                if(jproc.ne.-iproc) then
c Loop over selected initial states
                  do i=imin_inst(jproc),imax_inst(jproc)
                    jinst=i_allwd_inst(i,jproc)
c Loop over parton number
                    do ni=3,5
                      jproc0_2=j_prc_2_coll(jproc,jinst,ni)
                      jinst0_2=j_ins_2_coll(jproc,jinst,ni)
c If jproc0_2=-1 the splitting is not allowed. Skip the contribution
                      if(jinst.ne.-iinst.and.jproc0_2.ne.-1
     #                  .and.i_skip(jproc,jinst,ni).eq.0)then
c Statistical factor for this process
                        xstat_fact=i_stat_fact(jproc,jinst)
     #                            *i_skcnt(jproc,jinst,ni)/6.d0
                        call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,-1)
c Value of the matrix element for the physical process
                        call xmatel_five_part(sxi,xpp,xinv,xtmp,ytmp,
     #        2.d0,2.d0,jproc,jinst,ni,-1,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                        www=xfact_sc2_s0*x_s_fct
     #                     *zg2**2*ze2*wght*xstat_fact
c Store the results
                        index=i_sfmap(jproc,jinst)
                        wwcnt(1)=wwcnt(1)+sf(index,1)*www*xmte1q
                        wwcnt(2)=wwcnt(2)+sf(index,2)*www*xmte2q
                        wwcnt(3)=wwcnt(3)+sf(index,3)*www*xmte1e2
                      endif
c End of the loop over parton number
                    enddo
c End of the loop over initial states
                  enddo
                endif
c End of the loop over processes
              enddo
            endif
          endif
c End of counterevent S(0) (0,-1)
        endif
      endif
c End of the soft counterevents
c
c Output for the counterevents
      if(vvcnt(1).ne.0.d0.or.vvcnt(2).ne.0.d0.or.vvcnt(3).ne.0.d0 .or.
     #   wwcnt(1).ne.0.d0.or.wwcnt(2).ne.0.d0.or.wwcnt(3).ne.0.d0)then
        www=1.d0
        xtmp=0.d0
        ytmp=1.d0
        call invar_in(n1,n2,sxi,xtmp,ytmp,phii,yj,phij,ypp,yinv)
        if(iphot.eq.0)ycm=ycm+log(omega)
c
        if(vvcnt(1).ne.0.d0.or.vvcnt(2).ne.0.d0.or.vvcnt(3).ne.0.d0)then
          if(iphot.eq.0) then
            x1p=x1*omega/rxi
            call outall(www,tot,vvcnt(1),vvcnt(2),vvcnt(3),x1p)
          else
            x1p=1.d0
            ycm0 = ycm
            ycm = ycm + log(1-xii)/2.d0 
            call outall(www,tot,vvcnt(1),vvcnt(2),vvcnt(3),x1p)
            ycm = ycm0
          endif
          xsum_cntvv=tot/wght
        endif
c
        if(wwcnt(1).ne.0.d0.or.wwcnt(2).ne.0.d0.or.wwcnt(3).ne.0.d0)then
          if(iphot.eq.0) then
            x1p=x1*omega*rxi
            call outall(www,tot,wwcnt(1),wwcnt(2),wwcnt(3),x1p)
          else
            x1p=1.d0
            ycm0 = ycm
            ycm = ycm - log(1-xii)/2.d0
            call outall(www,tot,wwcnt(1),wwcnt(2),wwcnt(3),x1p)
            ycm = ycm0
          endif
          xsum_cntww=tot/wght
        endif
      endif
c
      sig1a=xsum_ev+xsum_cntvv+xsum_cntww
      return
      end


      function sig1bwb(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(5)
      sig1bwb=sig1b(xx,peso,0)
      return
      end


      function sig1bmon(xx,peso)
      implicit real * 8 (a-h,o-z)
      dimension xx(4)
      sig1bmon=sig1b(xx,peso,1)
      return
      end


      function sig1b(xx,wght,iphot)
c
      implicit real*8 (a-h,o-z)
      logical xpass,s1fun,s_2_fun,s1fun_soft
      parameter (pi=3.14159265358979312D0)
      parameter (phii=0.d0)
      parameter (tiny=1.d-6)
      parameter (xmslim=3.73205080756887729d0)
      parameter (n1=3)
      parameter (n2=4)
      dimension xx(*)
      common/shadr/sh
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/x1x2/ycm,tau
      common/samp/nsamp
      common/caso/caso(5)
      common/cutpar/xicut,yincut,youtcut
      common/irange/imax_inst,imin_inst
      common/iallwd/i_allwd_inst
      common/ijskip/ij_skip,ij_skcnt
      common/istat/i_stat_fact
      common/isfmap/i_sfmap
      common/imtt2/imtt2
      common/iinst/iinst
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/jphysproc/jproc,jinst
      common/icls/icls
      dimension sf(1:3,1:3),vvev(1:3),vvcnt(1:3)
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      dimension xinv(1:5),yinv(1:5)
      integer imax_inst(1:3),imin_inst(1:3)
      integer i_allwd_inst(1:2,1:3)
      integer ij_skip(1:3,1:2,3:5,3:5),ij_skcnt(1:3,1:2,3:5,3:5)
      integer i_stat_fact(1:3,1:2)
      integer i_sfmap(1:3,1:2)
      integer imtt2(3:5,3:5)
c
      icls=1
      xjac=1.d0
      roh = xsc_min2/sh
c 
      ximax0 = roh**(-nsamp)
      ximin0 = 1
      tmp  = ximin0 + xx(4)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
c 
      if(iphot.eq.0) then
        ymax= -log(tau)/2
        ymin=  log(tau)/2
        ycm = ymin + xx(5)*(ymax-ymin)
        xjac= xjac * (ymax-ymin)
      else
        ycm=log(1/tau)/2.d0
      endif
c 
      xii=(1.d0-tiny)*xx(1)**2
      xjac=xjac*xx(1)*2
c 
      ro=roh/tau
      yilim=sqrt(1.d0-ro)
      call jetchvar(xx(2),yi,yilim,xjac,ro)
c 
      zzz=pi*(tiny+(1.d0-tiny)*xx(3)**2)
      xjac=xjac*xx(3)*2*pi
      yj=cos(zzz)
      xjac=xjac*sin(zzz)
c caso(2) --> phij
      phij=2*pi*caso(2)
      xjac=xjac*2*pi
c s is the partonic CM energy
      s=sh*tau
c
      xnorm=xjac*s/(512*pi**4)
c
      xsum_ev=0.d0
      xsum_cnt=0.d0
      vvev(1)=0.d0
      vvev(2)=0.d0
      vvev(3)=0.d0
      vvcnt(1)=0.d0
      vvcnt(2)=0.d0
      vvcnt(3)=0.d0
c Evaluate the normalization for the non-soft configurations
      xphsp=2*(1-xii)/(2-xii*(1-yj))**2
      x_ns_fct=xnorm*xphsp
c Evaluate the normalization for the collinear configurations
      xphsp=(1-xii)/2
      x_c3_fct=xnorm*xphsp
c Evaluate the normalization for the soft configurations
      xphsp=1/2.d0
      x_s_fct=xnorm*xphsp
c
c Evaluate Bjorken variables
      if(iphot.eq.0) then
        x1=sqrt(tau) * exp(ycm)
        x2=tau/x1
      else
        x1=2.d0
        x2=tau
      endif
c
c Event
c
c Generate the kinematics
      call invar_out(n1,n2,s,xii,yi,phii,yj,phij,xij,ypp,yinv)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
      xpass=s1fun(n1,n2) .and. xij.gt.xii
      if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
        zg2=zgmu2()
c Calculate parton luminosity
        call strfun(x1,x2,sf)
c Factor appearing in eq.(4.65) and phase-space s
        xfact_ev=(1.d0/xii)*(1.d0/(1-yj))
c Loop over selected processes
        do jproc=min_proc,max_proc
          if(jproc.ne.-iproc) then
c Loop over selected initial states
            do i=imin_inst(jproc),imax_inst(jproc)
              jinst=i_allwd_inst(i,jproc)
c Loop over parton numbers
              do ni=3,5
              do nj=3,5
              if(ni.ne.nj)then
                if(jinst.ne.-iinst
     #             .and.ij_skip(jproc,jinst,ni,nj).eq.0)then
c Statistical factor for this process
                  xstat_fact=i_stat_fact(jproc,jinst)
     #                      *ij_skcnt(jproc,jinst,ni,nj)/6.d0
                  call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,nj)
c Value of the matrix element for the physical process
                  call xmatel_five_part(s,xpp,xinv,xii,2.d0,
     #         xij,yj,jproc,jinst,ni,nj,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                  www=xfact_ev*x_ns_fct
     #               *zg2**2*ze2*wght*xstat_fact
c Store the results
                  index=i_sfmap(jproc,jinst)
                  vvev(1)=vvev(1)+sf(index,1)*www*xmte1q
                  vvev(2)=vvev(2)+sf(index,2)*www*xmte2q
                  vvev(3)=vvev(3)+sf(index,3)*www*xmte1e2
                endif
              endif
c End of the loop over parton numbers
              enddo
              enddo
c End of the loop over initial states
            enddo
          endif
c End of the loop over processes
        enddo
      endif
c Output for the event
      if(vvev(1).ne.0.d0.or.vvev(2).ne.0.d0.or.vvev(3).ne.0.d0)then
        if(iphot.eq.0) then
          x1p=sqrt(tau) * exp(ycm)
        else
          x1p=1.d0
        endif
        www=1.d0
        call outall(www,tot,vvev(1),vvev(2),vvev(3),x1p)
        xsum_ev=tot/wght
      endif
c End of the event contribution
c
c Counterevent S(1) (xii,1)
c
      if(yj.gt.(1-youtcut))then
        ytmp=1.d0
c Generate the kinematics
        call invar_out(n1,n2,s,xii,yi,phii,ytmp,phij,xij,ypp,yinv)
c Define the variables of parton number 6(=n1+n2)
        call xjoin(n1,n2)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
        xpass=s_2_fun(n1,n2) .and. xij.gt.xii
        if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
          zg2=zgmu2()
c Calculate parton luminosity
          call strfun(x1,x2,sf)
c Factor appearing in eq.(4.65) and phase-space s
          xfact_c3=-(1.d0/xii)*(1.d0/(1-yj))
c Loop over selected processes
          do jproc=min_proc,max_proc
            if(jproc.ne.-iproc) then
c Loop over selected initial states
              do i=imin_inst(jproc),imax_inst(jproc)
                jinst=i_allwd_inst(i,jproc)
c Loop over parton numbers
                do ni=3,5
                do nj=3,5
                if(ni.ne.nj)then
                  if(jinst.ne.-iinst
     #               .and.ij_skip(jproc,jinst,ni,nj).eq.0)then
c Statistical factor for this process
                    xstat_fact=i_stat_fact(jproc,jinst)
     #                        *ij_skcnt(jproc,jinst,ni,nj)/6.d0
                    call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,nj)
c Value of the matrix element for the physical process
                    call xmatel_five_part(s,xpp,xinv,xii,2.d0,
     #           xij,ytmp,jproc,jinst,ni,nj,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                    www=xfact_c3*x_c3_fct
     #                 *zg2**2*ze2*wght*xstat_fact
c Store the results
                    index=i_sfmap(jproc,jinst)
                    vvcnt(1)=vvcnt(1)+sf(index,1)*www*xmte1q
                    vvcnt(2)=vvcnt(2)+sf(index,2)*www*xmte2q
                    vvcnt(3)=vvcnt(3)+sf(index,3)*www*xmte1e2
                  endif
                endif
c End of the loop over parton numbers
                enddo
                enddo
c End of the loop over initial states
              enddo
            endif
c End of the loop over processes
          enddo
        endif
      endif
c End of counterevent S(1) (xii,1)
      if(xii.lt.xicut)then
c
c S(1) (0,yj) counterevent
c
        xtmp=0.d0
c Generate the kinematics
        call invar_out(n1,n2,s,xtmp,yi,phii,yj,phij,xij,ypp,yinv)
c Define the variables of parton number 6(=n1+n2)
        call xjoin(n1,n2)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
        xpass=s_2_fun(n1,n2) .and. s1fun_soft(n1,n2)
        if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
          zg2=zgmu2()
c Calculate parton luminosity
          call strfun(x1,x2,sf)
c Factor appearing in eq.(4.65) and phase-space s
          xfact_s=-(1.d0/xii)*(1.d0/(1-yj))
c Loop over selected processes
          do jproc=min_proc,max_proc
            if(jproc.ne.-iproc) then
c Loop over selected initial states
              do i=imin_inst(jproc),imax_inst(jproc)
                jinst=i_allwd_inst(i,jproc)
c Loop over parton numbers
                do ni=3,5
                do nj=3,5
                if(ni.ne.nj)then
                  if(jinst.ne.-iinst
     #               .and.ij_skip(jproc,jinst,ni,nj).eq.0)then
c Statistical factor for this process
                    xstat_fact=i_stat_fact(jproc,jinst)
     #                        *ij_skcnt(jproc,jinst,ni,nj)/6.d0
                    call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,nj)
c Value of the matrix element for the physical process
                    call xmatel_five_part(s,xpp,xinv,xtmp,2.d0,
     #           xij,yj,jproc,jinst,ni,nj,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                    www=xfact_s*x_s_fct
     #                 *zg2**2*ze2*wght*xstat_fact
c Store the results
                    index=i_sfmap(jproc,jinst)
                    vvcnt(1)=vvcnt(1)+sf(index,1)*www*xmte1q
                    vvcnt(2)=vvcnt(2)+sf(index,2)*www*xmte2q
                    vvcnt(3)=vvcnt(3)+sf(index,3)*www*xmte1e2
                  endif
                endif
c End of the loop over parton numbers
                enddo
                enddo
c End of the loop over initial states
              enddo
            endif
c End of the loop over processes
          enddo
        endif
c End of counterevent S(1) (0,yj) 
c
c Counterevent S(1) (0,1)
c
        if(yj.gt.(1-youtcut))then
          ytmp=1.d0
c Generate the kinematics
          call invar_out(n1,n2,s,xtmp,yi,phii,ytmp,phij,xij,ypp,yinv)
c Define the variables of parton number 6(=n1+n2)
          call xjoin(n1,n2)
c Check whether the kinematics fulfills the cuts imposed by the
c measurement function. 
          xpass=s_2_fun(n1,n2)
          if(xpass)then
c Calculate the coupling costant using the partonic variables
c Set the factorization scale in the function zgmu2
            zg2=zgmu2()
c Calculate parton luminosity
            call strfun(x1,x2,sf)
c Factor appearing in eq.(4.65) and phase-space s
            xfact_sc3=(1.d0/xii)*(1.d0/(1-yj))
c Loop over selected processes
            do jproc=min_proc,max_proc
              if(jproc.ne.-iproc) then
c Loop over selected initial states
                do i=imin_inst(jproc),imax_inst(jproc)
                  jinst=i_allwd_inst(i,jproc)
c Loop over parton numbers
                  do ni=3,5
                  do nj=3,5
                  if(ni.ne.nj)then
                    if(jinst.ne.-iinst
     #                 .and.ij_skip(jproc,jinst,ni,nj).eq.0)then
c Statistical factor for this process
                      xstat_fact=i_stat_fact(jproc,jinst)
     #                          *ij_skcnt(jproc,jinst,ni,nj)/6.d0
                      call xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,nj)
c Value of the matrix element for the physical process
                      call xmatel_five_part(s,xpp,xinv,xtmp,2.d0,
     #             xij,ytmp,jproc,jinst,ni,nj,xmte1q,xmte2q,xmte1e2)
c Weight for the event
                      www=xfact_sc3*x_s_fct
     #                   *zg2**2*ze2*wght*xstat_fact
c Store the results
                      vvcnt(1)=vvcnt(1)+sf(jinst,1)*www*xmte1q
                      vvcnt(2)=vvcnt(2)+sf(jinst,2)*www*xmte2q
                      vvcnt(3)=vvcnt(3)+sf(jinst,3)*www*xmte1e2
                    endif
                  endif
c End of the loop over parton numbers
                  enddo
                  enddo
c End of the loop over initial states
                enddo
              endif
c End of the loop over processes
            enddo
          endif
c End of counterevent S(1) (0,1)
        endif
c End of soft and soft-collinear counterevents
      endif
c
c Output for the counterevent
      if(vvcnt(1).ne.0.d0.or.vvcnt(2).ne.0.d0.or.vvcnt(3).ne.0.d0)then
        if(iphot.eq.0) then
          x1p=sqrt(tau) * exp(ycm)
        else
          x1p=1.d0
        endif
        xtmp=0.d0
        ytmp=1.d0
        call invar_out(n1,n2,s,xtmp,yi,phii,ytmp,phij,xij,ypp,yinv)
        www=1.d0
        call outall(www,tot,vvcnt(1),vvcnt(2),vvcnt(3),x1p)
        xsum_cnt=tot/wght
      endif
c
      sig1b=xsum_ev+xsum_cnt
      return
      end
c
c End of next-to-leading-order contribution
c
c
c
c Begin of measurement functions
c
c   i_type_sfun=1 ==> Ellis-Soper iterative algorithm
c   i_type_sfun=2 ==> cone algorithm
c   i_type_sfun=3 ==> modified cone algorithm
c
c   i_type_merge=1 ==> weighted sum of eta and phi
c
c
      function s_2_fun(i,j)
c
c Eq. (2.27), for two partons in the final state
c
      implicit real * 8 (a-h,o-z)
      logical flag,s_2_fun
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/ffact/ffact,djet2
      common/partkin/xkt,xeta,xphi
      common/imtt/imtt
      common/imtt2/imtt2
      common/i_type_sfun/i_type_sfun
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt(3:4,3:5),imtt2(3:5,3:5)
c
      if(j.eq.-1)then
        n3=imtt(3,i)
        n4=imtt(4,i)
      else
        n3=imtt2(i,j)
        n4=6
      endif
      if(i_type_sfun.eq.1)then
        xra=rdel(xeta(n3),xphi(n3),xeta(n4),xphi(n4))
        xlwscl=(xkt(n3)+xkt(n4))**2
        flag=xlwscl.gt.xsc_min2 .and.
     #       xra.gt.djet2
      elseif(i_type_sfun.eq.2.or.i_type_sfun.eq.3)then
        xlwscl=(xkt(n3)+xkt(n4))**2
        flag=xlwscl.gt.xsc_min2
      else
        write(6,*)'Error in s_2_fun: algorithm not implemented'
        write(6,*)'Type #: ',i_type_sfun
        stop
      endif
      s_2_fun=flag
      return
      end


      function s0fun(i)
c
c Eq. (2.27), for three partons in the final state 
c
      implicit real * 8 (a-h,o-z)
      logical flag,s0fun
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/ffact/ffact,djet2
      common/partkin/xkt,xeta,xphi
      common/imtt/imtt
      common/i_type_sfun/i_type_sfun
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt(3:4,3:5)
c
      n3=imtt(3,i)
      n4=imtt(4,i)
      if(i_type_sfun.eq.1)then
        xmin=srcmin(3,5,-1,-1)
        xra=rdel(xeta(n3),xphi(n3),xeta(n4),xphi(n4))
        xlwscl=(xkt(n3)+xkt(n4)+xkt(i))**2
        flag=xmin.eq.xkt(i)**2   .and.
     #       xlwscl.gt.xsc_min2  .and.
     #       xra.gt.djet2
      elseif(i_type_sfun.eq.2.or.i_type_sfun.eq.3)then
        xra=rdel05(xeta(n3),xphi(n3),xeta(n4),xphi(n4))
        xda=conept(xkt(n3),xkt(n4))
        xrb=rdel05(xeta(n3),xphi(n3),xeta(i),xphi(i))
        xdb=conept(xkt(n3),xkt(i))
        xrc=rdel05(xeta(i),xphi(i),xeta(n4),xphi(n4))
        xdc=conept(xkt(i),xkt(n4))
        xlwscl=(xkt(n3)+xkt(n4)+xkt(i))**2
        xmin=min(xkt(n3),xkt(n4))
        flag=xra.gt.xda          .and.
     #       xrb.gt.xdb          .and.
     #       xrc.gt.xdc          .and.
     #       xmin.gt.xkt(i)      .and.
     #       xlwscl.gt.xsc_min2
      else
        write(6,*)'Error in s0fun: algorithm not implemented'
        write(6,*)'Type #: ',i_type_sfun
        stop
      endif
      s0fun=flag
      return
      end


      function s0fun_soft(i)
c This function, multiplied by S_2, gives the soft limit of 
c S^{(0)}_i 
      implicit real * 8 (a-h,o-z)
      logical flag,s0fun_soft
      common/ffact/ffact,djet2
      common/partkin/xkt,xeta,xphi
      common/imtt/imtt
      common/i_type_sfun/i_type_sfun
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt(3:4,3:5)
c
      n3=imtt(3,i)
      n4=imtt(4,i)
      if(i_type_sfun.eq.1)then
        xra=rdel(xeta(n3),xphi(n3),xeta(i),xphi(i))
        xrb=rdel(xeta(n4),xphi(n4),xeta(i),xphi(i))
        flag=xra.gt.djet2 .and.
     #       xrb.gt.djet2
      elseif(i_type_sfun.eq.2.or.i_type_sfun.eq.3)then
        xrb=rdel05(xeta(n3),xphi(n3),xeta(i),xphi(i))
        xrc=rdel05(xeta(i),xphi(i),xeta(n4),xphi(n4))
        flag=xrb.gt.djet2        .and.
     #       xrc.gt.djet2
      else
        write(6,*)'Error in s0fun_soft: algorithm not implemented'
        write(6,*)'Type #: ',i_type_sfun
        stop
      endif
      s0fun_soft=flag
      return
      end


      function s1fun(i,j)
c
c Eq. (2.28), for three partons in the final state 
c
      implicit real * 8 (a-h,o-z)
      logical flag,s1fun
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/ffact/ffact,djet2
      common/partkin/xkt,xeta,xphi
      common/imtt2/imtt2
      common/i_type_sfun/i_type_sfun
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt2(3:5,3:5)
c
      n3=imtt2(i,j)
      if(i_type_sfun.eq.1)then
        xmin=srcmin(3,5,-1,-1)
        xcmp=ddel(xkt(i),xeta(i),xphi(i),xkt(j),xeta(j),xphi(j))
        xlwscl=(xkt(n3)+xkt(i)+xkt(j))**2
        flag=xmin.eq.xcmp        .and.
     #       xlwscl.gt.xsc_min2
      elseif(i_type_sfun.eq.2.or.i_type_sfun.eq.3)then
        xra=rdel05(xeta(i),xphi(i),xeta(j),xphi(j))
        xda=conept(xkt(i),xkt(j))
        xrb=rdel05(xeta(n3),xphi(n3),xeta(i),xphi(i))
        xdb=conept(xkt(n3),xkt(i))
        xrc=rdel05(xeta(n3),xphi(n3),xeta(j),xphi(j))
        xdc=conept(xkt(n3),xkt(j))
        xlwscl=(xkt(n3)+xkt(i)+xkt(j))**2
        xmin=min(xkt(i),xkt(j))
        flag=xra.lt.xda          .and.
     #       xrb.gt.xdb          .and.
     #       xrc.gt.xdc          .and.
     #       xkt(n3).gt.xmin     .and.
     #       xlwscl.gt.xsc_min2
      else
        write(6,*)'Error in s1fun: algorithm not implemented'
        write(6,*)'Type #: ',i_type_sfun
        stop
      endif
      s1fun=flag
      return
      end


      function s1fun_soft(i,j)
c This function, multiplied by S_2, gives the soft limit of 
c S^{(1)}_{ij} 
      implicit real * 8 (a-h,o-z)
      logical flag,s1fun_soft
      common/ffact/ffact,djet2
      common/partkin/xkt,xeta,xphi
      common/imtt2/imtt2
      common/i_type_sfun/i_type_sfun
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt2(3:5,3:5)
c
      n3=imtt2(i,j)
      if(i_type_sfun.eq.1)then
        xra=rdel(xeta(i),xphi(i),xeta(j),xphi(j))
        xrb=rdel(xeta(n3),xphi(n3),xeta(i),xphi(i))
        flag=djet2.gt.xra .and.
     #       xrb.gt.xra
      elseif(i_type_sfun.eq.2.or.i_type_sfun.eq.3)then
        xra=rdel05(xeta(i),xphi(i),xeta(j),xphi(j))
        xrb=rdel05(xeta(n3),xphi(n3),xeta(i),xphi(i))
        flag=xra.lt.djet2        .and.
     #       xrb.gt.djet2
      else
        write(6,*)'Error in s1fun_soft: algorithm not implemented'
        write(6,*)'Type #: ',i_type_sfun
        stop
      endif
      s1fun_soft=flag
      return
      end


      function s2fun(i,j)
c
c Eq. (2.29), for three partons in the final state 
c
      implicit real * 8 (a-h,o-z)
      logical flag,s2fun
      common/fixvar/xsc_min2,xlam,xmuf2,xmur2,xmues2,xmuww2,zg,ze2
      common/ffact/ffact,djet2
      common/partkin/xkt,xeta,xphi
      common/imtt2/imtt2
      common/i_type_sfun/i_type_sfun
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt2(3:5,3:5)
c
      n3=imtt2(i,j)
      if(i_type_sfun.eq.1)then
        xmin=srcmin(3,5,-1,-1)
        xra=rdel(xeta(i),xphi(i),xeta(j),xphi(j))
        xlwscl=(xkt(n3)+xkt(i)+xkt(j))**2
        flag=xmin.eq.xkt(n3)**2  .and.
     #       xlwscl.gt.xsc_min2  .and.
     #       xra.lt.djet2
      elseif(i_type_sfun.eq.2.or.i_type_sfun.eq.3)then
        xra=rdel05(xeta(i),xphi(i),xeta(j),xphi(j))
        xda=conept(xkt(i),xkt(j))
        xrb=rdel05(xeta(n3),xphi(n3),xeta(i),xphi(i))
        xdb=conept(xkt(n3),xkt(i))
        xrc=rdel05(xeta(n3),xphi(n3),xeta(j),xphi(j))
        xdc=conept(xkt(n3),xkt(j))
        xlwscl=(xkt(n3)+xkt(i)+xkt(j))**2
        xmin=min(xkt(i),xkt(j))
        flag=xra.lt.xda          .and.
     #       xrb.gt.xdb          .and.
     #       xrc.gt.xdc          .and.
     #       xkt(n3).lt.xmin     .and.
     #       xlwscl.gt.xsc_min2
      else
        write(6,*)'Error in s2fun: algorithm not implemented'
        write(6,*)'Type #: ',i_type_sfun
        stop
      endif
      s2fun=flag
      return
      end


      subroutine xjoin(i,j)
c
c Combine the four-vectors of partons number i and j following the merging
c prescription and put the results into xkt(6), xeta(6), xphi(6)
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/partkin/xkt,xeta,xphi
      common/i_type_merge/i_type_merge
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
c
      if(i_type_merge.eq.1)then
        xkt(6)=xkt(i)+xkt(j)
        xeta(6)=(xkt(i)*xeta(i)+xkt(j)*xeta(j))/xkt(6)
        if(abs(xphi(i)-xphi(j)).lt.pi)then
          xphi(6)=(xkt(i)*xphi(i)+xkt(j)*xphi(j))/xkt(6)
        else
          if(xphi(j).lt.xphi(i))then
            xphi(6)=(xkt(i)*xphi(i)+xkt(j)*(xphi(j)+2*pi))/xkt(6)
          else
            xphi(6)=(xkt(i)*xphi(i)+xkt(j)*(xphi(j)-2*pi))/xkt(6)
          endif
        endif
        xphi(6)=atan2(sin(xphi(6)),cos(xphi(6)))
      else
        write(6,*)'Error in xjoin: algorithm not implemented'
        write(6,*)'Type #: ',i_type_merge
        stop
      endif
      return
      end


      function rdel(eta1,phi1,eta2,phi2)
c
c Eq. (2.22)
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
      dteta=eta1-eta2
      dtphi=phi1-phi2
      dtphi=dacos(dcos(dtphi))
      rdel=dteta**2+dtphi**2
      return
      end


      function rdel05(eta1,phi1,eta2,phi2)
c
c Square root of eq. (2.22)
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/i_type_sfun/i_type_sfun
c
      dteta=abs(eta1-eta2)
      if(i_type_sfun.eq.3)dteta=dteta*exp(dteta/10.d0)
      dtphi=phi1-phi2
      dtphi=dacos(dcos(dtphi))
      rdel05=dsqrt(dteta**2+dtphi**2)
      return
      end


      function conept(pt1,pt2)
c
c (pt1+pt2)*djet2/max(pt1,pt2), for the cone algorithm
c
      implicit real * 8 (a-h,o-z)
      common/ffact/ffact,djet2
c
      tmp=(pt1+pt2)/max(pt1,pt2)
      conept=tmp*djet2
      return
      end


      function ddel(pt1,eta1,phi1,pt2,eta2,phi2)
c
c Eq. (2.23)
c
      implicit real * 8 (a-h,o-z)
      common/ffact/ffact,djet2
c
      tmp=min(pt1**2,pt2**2)*rdel(eta1,phi1,eta2,phi2)
      ddel=tmp/djet2
      return
      end


      function srcmin(il,ih,iex,jex)
c
c This function evaluates the minimum among the d_i and d_{ij} of
c a given kinematical configuration, following Ellis-Soper prescription.
c The indices i and j run from il to ih, except for the values iex and jex
c
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/shadr/sh
      common/partkin/xkt,xeta,xphi
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
c
      xmin=sh
      xmin2=sh
      do i=il,ih
        if(i.ne.iex.and.(xkt(i)**2).lt.xmin)xmin=xkt(i)**2
      enddo
      do i=il,ih-1
        do j=i+1,ih
          ttt=ddel(xkt(i),xeta(i),xphi(i),xkt(j),xeta(j),xphi(j))
          if(i.ne.iex.and.j.ne.jex.and.ttt.lt.xmin2)xmin2=ttt
        enddo
      enddo
      ttt=min(xmin,xmin2)
      srcmin=ttt
      return
      end


      function s1_aux(i,j)
c Theta(d_j-d_i)
      implicit real * 8 (a-h,o-z)
      logical flag,s1_aux
      common/partkin/xkt,xeta,xphi
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
c
      flag=xkt(j).gt.xkt(i)
      s1_aux=flag
      return
      end
c
c End of measurement functions
c
c
c Begin of output routines
c
      subroutine outall(xw,tot,sf1,sf2,sf3,x1p)
c This subroutine boosts to the lab frame the protojets momenta
c for the D events. Notice that the conjugation procedure,
c which in principle would require the C events,
c is not needed due to the flavour blindness of the jets.
c Also, the symmetrization over the jet labels is carried out.
      implicit real * 8 (a-h,o-z)
      common/x1x2/ycm,tau
c GeV to microbarn conversion factor: sigma (mub) = hc2 * sigma (GeV^-2)
      parameter(hc2=3.8937966d2)
c
      w = xw * hc2
c Total for D events, to be returned to the calling function
      tot = 0
c Direct event
      if(sf1.ne.0.or.sf2.ne.0.or.sf3.ne.0) then
         dw = w * (sf1+sf2+sf3)
         tot  = tot + dw
c Boost to lab frame
         call labmom(ycm)
         call outfun(dw,x1p)
      endif
      end


      subroutine labmom(y)
c This subroutine boosts the momenta of the partons from the partonic
c CM frame to the laboratory frame
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      dimension xkt_cm(1:3),xeta_cm(1:3),xphi_cm(1:3)
      common/partkin/xkt,xeta,xphi
      common/cmkin/xkt_cm,xeta_cm,xphi_cm
c
      xkt_cm(1)=xkt(3)
      xkt_cm(2)=xkt(4)
      xkt_cm(3)=xkt(5)
      xeta_cm(1)=xeta(3)+y
      xeta_cm(2)=xeta(4)+y
      xeta_cm(3)=xeta(5)+y
      xphi_cm(1)=xphi(3)
      xphi_cm(2)=xphi(4)
      xphi_cm(3)=xphi(5)
      return
      end


      subroutine reflex(y)
c This subroutine has the same meaning of labmom, but it also
c performs the reflection
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      dimension xkt_cm(1:3),xeta_cm(1:3),xphi_cm(1:3)
      common/partkin/xkt,xeta,xphi
      common/cmkin/xkt_cm,xeta_cm,xphi_cm
c
      xkt_cm(1)=xkt(3)
      xkt_cm(2)=xkt(4)
      xkt_cm(3)=xkt(5)
      xeta_cm(1)=-xeta(3)-y
      xeta_cm(2)=-xeta(4)-y
      xeta_cm(3)=-xeta(5)-y
      xphi_cm(1)=dsign(pi,-xphi(3))+xphi(3)
      xphi_cm(2)=dsign(pi,-xphi(4))+xphi(4)
      xphi_cm(3)=dsign(pi,-xphi(5))+xphi(5)
      return
      end
c
c End of output routines
c
c
c Begin of phase-space routines
c
      subroutine invar_in(ni,nj,s,xii,yi,phii,yj,phij,xpp,xinv)
      implicit real * 8 (a-z)
      integer ni,nj,n3,i,j
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-10)
      common/partkin/xkt,xeta,xphi
      common/imtt2/imtt2
      dimension xpp(1:5,1:5),xinv(1:5)
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt2(3:5,3:5)
c
      do i=1,5
        do j=1,5
          xpp(i,j)=0.d0
        enddo
      enddo
      xpp(1,2)=s/2.d0
c
      n3=imtt2(ni,nj)
c
      xkt(ni)=0.d0
      xphi(ni)=atan2(sin(phii),cos(phii))
      if(yi.eq.1.d0)then
        xeta(ni)=1.d8
      elseif(yi.eq.-1.d0)then
        xeta(ni)=-1.d8
      else
        xkt(ni)=sqrt(s)*xii*sqrt(1-yi**2)/2.d0
        xeta(ni)=0.5d0*log((1.d0+yi)/(1.d0-yi))
      endif
      xpp(1,ni)=-s*xii*(1-yi)/4.d0
      xpp(2,ni)=-s*xii*(1+yi)/4.d0
      xpp(nj,n3)=s*(1-xii)/2.d0
c
      cpsi2_2 = 2*sqrt(1-xii)/sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
      spsi2_2 = xii*sqrt(1-yi**2)/sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
      cpsi1_2 = 2*yi*sqrt(1-xii)/sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
      spsi1_2 = (2-xii)*sqrt(1-yi**2)/
     #          sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
      cpsi_2  = -cpsi2_2
      spsi_2  = spsi2_2
c
      xpp(1,nj)=-s*(2-xii*(1-yi))*
     #      (1-spsi2_2*cos(phij)*sqrt(1-yj**2)-cpsi2_2*yj)/8.d0
      xpp(2,n3)=-s*(2-xii*(1+yi))*
     #      (1+spsi_2*cos(phij)*sqrt(1-yj**2)+cpsi_2*yj)/8.d0
      xpp(1,n3)=-s*(2-xii*(1-yi))*
     #      (1+spsi2_2*cos(phij)*sqrt(1-yj**2)+cpsi2_2*yj)/8.d0
      xpp(2,nj)=-s*(2-xii*(1+yi))*
     #      (1-spsi_2*cos(phij)*sqrt(1-yj**2)-cpsi_2*yj)/8.d0
      xpp(ni,nj)=s/4.d0*xii*
     #      (1-spsi1_2*cos(phij)*sqrt(1-yj**2)-cpsi1_2*yj)
      xpp(ni,n3)=s/4.d0*xii*
     #      (1+spsi1_2*cos(phij)*sqrt(1-yj**2)+cpsi1_2*yj)
c
      do i=1,4
        do j=i+1,5
          if(xpp(i,j).eq.0.d0)then
            xpp(i,j)=xpp(j,i)
          else
            xpp(j,i)=xpp(i,j)
          endif
        enddo
      enddo
c
      xkt(nj)=sqrt(s)/4.d0*sqrt( 
     #   ( sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
     #    -xii*cos(phij)*sqrt(1-yi**2)*sqrt(1-yj**2) )**2 
     #  - 4*yj**2*(1-xii)      )
      xkt(n3)=sqrt(s)/4.d0*sqrt( 
     #   ( sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
     #    +xii*cos(phij)*sqrt(1-yi**2)*sqrt(1-yj**2) )**2 
     #  - 4*yj**2*(1-xii)      )
      do i=3,5
        if(i.ne.ni)then
          if(xpp(2,i).eq.0.d0)then
            xeta(i)=-1.d-8
          elseif(xpp(1,i).eq.0.d0)then
            xeta(i)=1.d-8
          else
            xeta(i)=0.5d0*log(xpp(2,i)/xpp(1,i))
          endif
        endif
      enddo
c
      xc01_2=( -xii*sqrt(1-yi**2)+sqrt(1-yj**2)*cos(phij)*
     #          sqrt((2-xii*(1+yi))*(2-xii*(1-yi))) )/
     #        sqrt( ( sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
     #              -xii*cos(phij)*sqrt(1-yi**2)*sqrt(1-yj**2) )**2 
     #             - 4*yj**2*(1-xii) )
      xc02_2=( -xii*sqrt(1-yi**2)-sqrt(1-yj**2)*cos(phij)*
     #          sqrt((2-xii*(1+yi))*(2-xii*(1-yi))) )/
     #        sqrt( ( sqrt((2-xii*(1+yi))*(2-xii*(1-yi)))
     #              +xii*cos(phij)*sqrt(1-yi**2)*sqrt(1-yj**2) )**2 
     #             - 4*yj**2*(1-xii) )
      if(xc01_2.gt.1.d0-tiny)then
        xphi(nj)=0.d0
      elseif(xc01_2.lt.-1.d0+tiny)then
        xphi(nj)=pi
      else
        xphi(nj)=dacos(xc01_2)
      endif
      if(xc02_2.gt.1.d0-tiny)then
        xphi(n3)=2*pi
      elseif(xc02_2.lt.-1.d0+tiny)then
        xphi(n3)=pi
      else
        xphi(n3)=2*pi-dacos(xc02_2)
      endif
      if(phij.gt.pi)then
        xphi(nj)=2*pi-xphi(nj)
        xphi(n3)=2*pi-xphi(n3)
      endif
      xphi(nj)=atan2(sin(xphi(nj)+phii),cos(xphi(nj)+phii))
      xphi(n3)=atan2(sin(xphi(n3)+phii),cos(xphi(n3)+phii))
c
      xinv(1)=sqrt(s)*(1.d0-yi)/2.d0
      xinv(2)=sqrt(s)*(1.d0+yi)/2.d0
      xinv(ni)=0.d0
      xinv(nj)=sqrt(s)/2.d0*
     #        ( 1.d0-sqrt(1-yj**2)*cos(phij)*spsi1_2-yj*cpsi1_2 )
      xinv(n3)=sqrt(s)/2.d0*
     #        ( 1.d0+sqrt(1-yj**2)*cos(phij)*spsi1_2+yj*cpsi1_2 )
c
      return
      end


      subroutine invar_out(ni,nj,s,xii,yi,phii,yj,phij,xij,xpp,xinv)
      implicit real * 8 (a-z)
      integer ni,nj,n3,i,j
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-18)
      common/partkin/xkt,xeta,xphi
      common/imtt2/imtt2
      dimension xpp(1:5,1:5),xinv(1:5)
      dimension xkt(3:6),xeta(3:6),xphi(3:6)
      integer imtt2(3:5,3:5)
c
      n3=imtt2(ni,nj)
      xij=2*(1-xii)/(2-xii*(1-yj))
c
      xkt(nj)=0.d0
      xphi(nj)=atan2(sin(phii),cos(phii))
      if(yi.eq.1.d0)then
        xeta(nj)=1.d8
      elseif(yi.eq.-1.d0)then
        xeta(nj)=-1.d8
      else
        xkt(nj)=sqrt(s)*xij*sqrt(1-yi**2)/2.d0
        xeta(nj)=0.5d0*log((1.d0+yi)/(1.d0-yi))
      endif
c
      xx=cos(phij)*sqrt(1-yj**2)
      yy=sin(phij)*sqrt(1-yj**2)
      zz=yj
      si=sqrt(1-yi**2)
      cphi=cos(phii)
      sphi=sin(phii)
      aa1=xx*yi*cphi-yy*sphi+zz*si*cphi
      bb1=xx*yi*sphi+yy*cphi+zz*si*sphi
      cc1=-xx*si+zz*yi
      xkt(ni)=0.d0
      xphi(ni)=atan2(bb1,aa1)
      if(cc1.eq.1.d0)then
        xeta(ni)=1.d8
      elseif(cc1.eq.-1.d0)then
        xeta(ni)=-1.d8
      else
        xkt(ni)=sqrt(s)*xii*sqrt(aa1**2+bb1**2)/2.d0
        xeta(ni)=0.5d0*log((1.d0+cc1)/(1.d0-cc1))
      endif
c
      dd2=2.d0-xii-xij
      aa2=-(xij*cos(phii)*sqrt(1-yi**2)+xii*aa1)
      bb2=-(xij*sin(phii)*sqrt(1-yi**2)+xii*bb1)
      cc2=-(xij*yi+xii*cc1)
      xkt(n3)=0.d0
      xphi(n3)=atan2(bb2,aa2)
      if(dd2.eq.cc2)then
        xeta(n3)=1.d8
      elseif(dd2.eq.-cc2)then
        xeta(n3)=-1.d8
      else
        xeta(n3)=0.5d0*log((dd2+cc2)/(dd2-cc2))
        xkt(n3)=sqrt(s)*sqrt(aa2**2+bb2**2)/2.d0
      endif
c
      do i=1,5
        do j=1,5
          xpp(i,j)=0.d0
        enddo
      enddo
      xpp(1,2)=s/2.d0
      xpp(1,nj)=-s*xij*(1-yi)/4.d0
      xpp(1,ni)=-s*xii*(1-cc1)/4.d0
      xpp(2,nj)=-s*xij*(1+yi)/4.d0
      xpp(2,ni)=-s*xii*(1+cc1)/4.d0
c
      xpp(ni,nj)=s*xii*xij*(1-yj)/4.d0
      xpp(1,n3)=xpp(2,ni)+xpp(2,nj)+xpp(ni,nj)
      xpp(2,n3)=xpp(1,ni)+xpp(1,nj)+xpp(ni,nj)
      xpp(ni,n3)=xpp(1,2)+xpp(1,nj)+xpp(2,nj)
      xpp(nj,n3)=s*(1-xii)/2.d0
c
      do i=1,4
        do j=i+1,5
          if(xpp(i,j).eq.0.d0)then
            xpp(i,j)=xpp(j,i)
          else
            xpp(j,i)=xpp(i,j)
          endif
        enddo
      enddo
c
      xinv(1)=sqrt(s)*(1.d0-cc1)/2.d0
      xinv(2)=sqrt(s)*(1.d0+cc1)/2.d0
      xinv(ni)=0.d0
      xinv(nj)=sqrt(s)*xij
      xinv(n3)=sqrt(s)*( dd2-aa1*aa2-bb1*bb2-cc1*cc2 )/2.d0
c
      return
      end


      subroutine xpp_reorder(ypp,xpp,yinv,xinv,n1,n2,ni,nj)
c This subroutine performs the permutation over the invariants xpp
c defined by 
c                 n1 --> ni,      n2 --> nj.
c If nj=-1, nj is redefined as the smallest integer (>=3) not
c equal to ni. The inputs are ypp and yinv, the outputs xpp and xinv
      implicit real*8 (a-h,o-z)
      common/imtt/imtt
      common/imtt2/imtt2
      dimension xpp(1:5,1:5),ypp(1:5,1:5)
      dimension xinv(1:5),yinv(1:5)
      integer imtt(3:4,3:5)
      integer imtt2(3:5,3:5)
c
      ia=ni
      ib=nj
      if(ib.eq.-1)ib=imtt(3,ia)
      if(ia.eq.n1.and.ib.eq.n2)then
        do ii=1,5
          do jj=1,5
            xpp(ii,jj)=ypp(ii,jj)
          enddo
        enddo
        do ii=1,5
          xinv(ii)=yinv(ii)
        enddo
      else
        do ii=1,5
          do jj=1,5
            xpp(ii,jj)=0.d0
          enddo
        enddo
        ic=imtt2(ia,ib)
        n3=imtt2(n1,n2)
        xpp(1,2)=ypp(1,2)
        xpp(1,ia)=ypp(1,n1)
        xpp(1,ib)=ypp(1,n2)
        xpp(1,ic)=ypp(1,n3)
        xpp(2,ia)=ypp(2,n1)
        xpp(2,ib)=ypp(2,n2)
        xpp(2,ic)=ypp(2,n3)
        xpp(ia,ib)=ypp(n1,n2)
        xpp(ia,ic)=ypp(n1,n3)
        xpp(ib,ic)=ypp(n2,n3)
        do ii=1,4
          do jj=ii+1,5
            if(xpp(ii,jj).eq.0.d0)then
              xpp(ii,jj)=xpp(jj,ii)
            else
              xpp(jj,ii)=xpp(ii,jj)
            endif
          enddo
        enddo
        xinv(1)=yinv(1)
        xinv(2)=yinv(2)
        xinv(ia)=yinv(n1)
        xinv(ib)=yinv(n2)
        xinv(ic)=yinv(n3)
      endif
c
      return
      end


      subroutine jetchvar(parth1,cth1,cthlim,xjac,ro)
c
c Given 0<parth1<1 returns -1<cth1<1
c and multiplies xjac times the d cth1 / d parth1 jacobian
c
      implicit real * 8 (a-z)
      ylim=1-ro/4
      zz=2*parth1-1
      xjac=xjac*2
      if(zz.gt.ylim)then
        cth1 = (1-cthlim)/(1-ylim)*zz+(cthlim-ylim)/(1-ylim)
        xjac = xjac*(1-cthlim)/(1-ylim)
      elseif(zz.lt.-ylim)then
        cth1 = (1-cthlim)/(1-ylim)*zz-(cthlim-ylim)/(1-ylim)
        xjac = xjac*(1-cthlim)/(1-ylim)
      else
        bb = 1-ro**2/16
        xlgbb = log((1+bb)/(1-bb))
        yy = zz*4/(4-ro)
        xjac = xjac*4/(4-ro)
        expyy = exp(-yy*xlgbb)
        cth1 = cthlim*(1-expyy)/(1+expyy)/bb
        xjac = xjac*2*cthlim*xlgbb*expyy/(1+expyy)**2/bb
      endif
      return
      end


      subroutine jetchvar1(parth1,cth1,cthlim,xjac,ro)
c
c Given 0<parth1<1 returns -1<cth1<1
c and multiplies xjac times the d cth1 / d parth1 jacobian
c
      implicit real * 8 (a-z)
      bb = 1-ro**2/16
      xlgbb = log((1+bb)/(1-bb))
      yy = 2*parth1-1
      xjac = xjac*2
      ylim=1-ro/4
      if(yy.gt.ylim)then
        x0=ylim
        expyy=exp(-x0*xlgbb)
        y0=cthlim*(1-expyy)/(1+expyy)/bb
        y0p=2*cthlim*xlgbb*expyy/(1+expyy)**2/bb
        xa=(1-y0p*(1-x0)-y0)/(1-x0)**2
        xb=y0p-2*xa*x0
        xc=1-xa-xb
        cth1=xa*yy**2+xb*yy+xc
        xjac=xjac*(2*xa*yy+xb)
      elseif(yy.lt.-ylim)then
        x0=-ylim
        expyy=exp(-x0*xlgbb)
        y0=cthlim*(1-expyy)/(1+expyy)/bb
        y0p=2*cthlim*xlgbb*expyy/(1+expyy)**2/bb
        xa=(-1+y0p*(1+x0)-y0)/(1+x0)**2
        xb=y0p-2*xa*x0
        xc=-1-xa+xb
        cth1=xa*yy**2+xb*yy+xc
        xjac=xjac*(2*xa*yy+xb)
      else
        expyy = exp(-yy*xlgbb)
        cth1 = cthlim*(1-expyy)/(1+expyy)/bb
        xjac = xjac*2*cthlim*xlgbb*expyy/(1+expyy)**2/bb
      endif
      return
      end


      subroutine zzchvar(parth1,cth1,xjac,ro)
c
c Given 0<parth1<1 returns -1<cth1<1
c and multiplies xjac times the d cth1 / d parth1 jacobian
c
      implicit real * 8 (a-z)
      bb = 1-ro**2/16
      xlgbb = log((1+bb)/(1-bb))
      yy = ( parth1 * 2 - 1 ) * xlgbb
      xjac = xjac * 2 * xlgbb
      expyy = exp(-yy)
      cth1 = (1-expyy)/(1+expyy)/bb
      xjac = xjac * 2 * expyy/(1+expyy)**2 / bb
      return
      end
c
c End of phase-space routines
c
c
c Begin of initialization routines
c
      subroutine setpar
c In this subroutine we initialize all the quantities which will 
c be used in the program
      implicit real*8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (xnc=3.d0)
      parameter (zup=2.d0/3.d0,zdo=-1.d0/3.d0)
      common/nl/nl
      common/numofproc/num_of_proc,num0_of_proc
      common/iallwd/i_allwd_inst
      common/i0allwd/i0_allwd_inst
      common/i2pvskip/i2pv_skip
      common/iskip/i_skip,i_skcnt
      common/ijskip/ij_skip,ij_skcnt
      common/istat/i_stat_fact
      common/i0stat/i0_stat_fact
      common/i_type_part/i_type_part
      common/i_type_in/i_type_in
      common/i_ap_in/i_ap_in1,i_ap_in2
      common/i_ap_out/i_ap_out
      common/jcoll_prc/j_prc_1_coll,j_prc_2_coll
      common/jcoll_ins/j_ins_1_coll,j_ins_2_coll
      common/jcoll_out_prc/j_prc_coll_out
      common/i_phys/i_phys
      common/i0_phys/i0_phys
      common/i_xcol/i_xcol
      common/isfmap/i_sfmap
      common/imtt/imtt
      common/imtt2/imtt2
      common/color_factors/xca,xgmm,xgmmprime
      common/charges/zel,sumzel,sumzel2
      integer num_of_proc(1:3),i_allwd_inst(1:2,1:3)
      integer num0_of_proc(1:3),i0_allwd_inst(1:2,1:3)
      integer i2pv_skip(1:3,1:2)
      integer i_skip(1:3,1:2,3:5),i_skcnt(1:3,1:2,3:5)
      integer ij_skip(1:3,1:2,3:5,3:5),ij_skcnt(1:3,1:2,3:5,3:5)
      integer i_stat_fact(1:3,1:2),i0_stat_fact(1:1,1:2)
      integer i_type_part(1:3,1:2,3:5),i_type_in(1:2,1:2)
      integer i_ap_in1(1:3,1:2,3:5),i_ap_in2(1:3,1:2,3:5)
      integer i_ap_out(1:3,1:2,3:5,3:5)
      integer j_prc_1_coll(1:3,1:2,3:5),j_prc_2_coll(1:3,1:2,3:5)
      integer j_ins_1_coll(1:3,1:2,3:5),j_ins_2_coll(1:3,1:2,3:5)
      integer j_prc_coll_out(1:3,1:2,3:4,4:5)
      integer i_phys(1:5,1:3,1:2),i0_phys(1:4,1:6,1:7)
      integer i_xcol(1:4,1:3,1:2)
      integer i_sfmap(1:3,1:2)
      integer imtt(3:4,3:5),imtt2(3:5,3:5)
      dimension xca(0:2),xgmm(0:2),xgmmprime(0:2)
      dimension zel(1:6),sumzel(1:6),sumzel2(1:6)
c
c The matrix NUM_OF_PROC(IPROC) returns the number of 2 --> 3
c physical processes associated with the unphysical 0 --> 5
c configuration defined by IPROC
      num_of_proc(1)=2
      num_of_proc(2)=1
      num_of_proc(3)=1
c
c The matrix I_ALLWD_INST(N,IPROC), with 1<=N<=NUM_OF_PROC(IPROC),
c returns the classification number of the n-th initial state
c allowed for the process defined by IPROC
      do i=1,3
        do j=1,2
          i_allwd_inst(j,i)=0
        enddo
      enddo
      i_allwd_inst(1,1)=1
      i_allwd_inst(2,1)=2
      i_allwd_inst(1,2)=2
      i_allwd_inst(1,3)=2
c
c The matrix NUM0_OF_PROC(IPROC) has the same meaning of
c NUM_OF_PROC, but is relevant for the 2 --> 2 processes.
c Only the 2 --> 2 processes initiated by a photon are
c relevant in this case
      num0_of_proc(1)=2
      num0_of_proc(2)=-1
      num0_of_proc(3)=-1
c
c The matrix I0_ALLWD_INST(N,IPROC), with 1<=N<=NUM0_OF_PROC(IPROC),
c has the same meaning of I_ALLWD_INST, but is relevant for 
c the 2 --> 2 processes
      do i=1,3
        do j=1,2
          i0_allwd_inst(j,i)=0
        enddo
      enddo
      i0_allwd_inst(1,1)=1
      i0_allwd_inst(2,1)=2
c
c The matrix I2PV_SKIP(IPROC,IINST) returns -1 if process defined 
c by (IPROC,IINST) is not allowed at the leading order
      do i=1,3
        do j=1,2
          i2pv_skip(i,j)=-1
        enddo
      enddo
      i2pv_skip(1,1)=0
      i2pv_skip(1,2)=0
c
c The matrix I_SKIP(IPROC,IINST,N) returns -1 if the contribution
c of the parton number N to the process defined by (IPROC,IINST) is
c trivial, in the sense that it has been already taken into account. This 
c happens when the final state contains two or more identical particles
      do i=1,3
        do j=1,2
          do k=3,5
            i_skip(i,j,k)=0
          enddo
        enddo
      enddo
      i_skip(1,2,5)=-1
      i_skip(3,2,4)=-1
c
c The matrix I_SKCNT(IPROC,IINST,N) returns I+1, where I is the number of
c the final state partons identical to the parton number N in the 
c process defined by (IPROC,IINST)
      do i=1,3
        do j=1,2
          do k=3,5
            i_skcnt(i,j,k)=1
          enddo
        enddo
      enddo
      i_skcnt(1,2,4)=2
      i_skcnt(3,2,3)=2
c
c The matrix IJ_SKIP(IPROC,IINST,N1,N2) returns -1 if the contribution
c of the pair of partons number N1 and N2 to the process defined by
c (IPROC,IINST) is trivial, in the sense that it has been already taken 
c into account
      do i=1,3
        do j=1,2
          do k1=3,5
            do k2=3,5
              ij_skip(i,j,k1,k2)=0
            enddo
          enddo
        enddo
      enddo
      ij_skip(1,2,3,5)=-1
      ij_skip(1,2,5,3)=-1
      ij_skip(1,2,5,4)=-1
      ij_skip(3,2,4,5)=-1
      ij_skip(3,2,4,3)=-1
      ij_skip(3,2,5,4)=-1
c
c The matrix IJ_SKCNT(IPROC,IINST,N1,N2) returns I+1, where I is the 
c number of the final state parton pairs identical to the pair (N1,N2)
c in the process defined by (IPROC,IINST). The pairs (N1,N2) and (N2,N1)
c are considered equivalent if the flavour of the partons number N1
c and N2 is the same
      do i=1,3
        do j=1,2
          do k1=3,5
            do k2=3,5
              ij_skcnt(i,j,k1,k2)=1
            enddo
          enddo
        enddo
      enddo
      ij_skcnt(1,2,3,4)=2
      ij_skcnt(1,2,4,5)=2
      ij_skcnt(1,2,4,3)=2
      ij_skcnt(3,2,3,4)=2
      ij_skcnt(3,2,3,5)=2
      ij_skcnt(3,2,5,3)=2
c
c The matrix I_STAT_FACT(IPROC,IINST) returns the statistical factor
c relevant for the 2 --> 3 processes defined by (IPROC,IINST). This is
c the number of non-identical permutations over the final state 
c partons of the process (IPROC,IINST). A factor dependent upon
c the number of flavours, which is present in the hadronic code,
c is in this case included via the sum over electric charges, performed
c in the subroutine STRFUN
      do i=1,3
        do j=1,2
          i_stat_fact(i,j)=0
        enddo
      enddo
      i_stat_fact(1,1)=6
      i_stat_fact(1,2)=3
      i_stat_fact(2,2)=6
      i_stat_fact(3,2)=3
c
c The matrix I0_STAT_FACT(IPROC,IINST) has the same meaning of
c I_STAT_FACT, but is relevant for the 2 --> 2 processes
      do i=1,1
        do j=1,2
          i0_stat_fact(i,j)=0
        enddo
      enddo
      i0_stat_fact(1,1)=2
      i0_stat_fact(1,2)=2
c
c The matrix I_TYPE_PART(IPROC,IINST,NPART) returns 0 if the parton
c number NPART in the process defined by (IPROC,IINST) is a gluon,
c and returns 1 otherwise. This is strictly dependent upon the ordering
c convention for the final state partons. See the comments above the
c subroutine XMATEL_FIVE_PART in the file PHDIJETCRS.FOR
      do i=1,3
        do j=1,2
          do k=3,5
            i_type_part(i,j,k)=1
          enddo
        enddo
      enddo
      i_type_part(1,1,5)=0
      i_type_part(1,2,4)=0
      i_type_part(1,2,5)=0
c
c The matrix I_TYPE_IN(IINST,NPART) has the same meaning of I_TYPE_PART,
c but is relevant for initial state partons. If the parton is a photon,
c the value 2 is returned
      i_type_in(1,1)=2
      i_type_in(1,2)=0
      i_type_in(2,1)=2
      i_type_in(2,2)=1
c
c The matrix J_PRC_1_COLL(IPROC,IINST,NI) returns the 2 --> 2 process
c number which factorizes when the parton number NI is collinear
c to the parton number 1. The matrix J_INS_1_COLL returns the
c corresponding initial state number. The matrices J_PRC_2_COLL
c and J_INS_2_COLL have the same meaning, but for collinear
c emission from parton number 2
      do i=1,3
        do j=1,2
          do k=3,5
            j_prc_1_coll(i,j,k)=-1
            j_prc_2_coll(i,j,k)=-1
          enddo
        enddo
      enddo
      j_prc_1_coll(1,1,3)=4
      j_prc_1_coll(1,1,4)=4
      j_prc_1_coll(1,2,3)=4
      j_prc_1_coll(2,2,3)=5
      j_prc_1_coll(2,2,4)=5
      j_prc_1_coll(2,2,5)=5
      j_prc_1_coll(3,2,3)=6
      j_prc_1_coll(3,2,4)=6
      j_prc_1_coll(3,2,5)=6
c
      j_prc_2_coll(1,1,3)=1
      j_prc_2_coll(1,1,4)=1
      j_prc_2_coll(1,1,5)=1
      j_prc_2_coll(1,2,4)=1
      j_prc_2_coll(1,2,5)=1
      j_prc_2_coll(2,2,3)=1
      j_prc_2_coll(3,2,3)=1
      j_prc_2_coll(3,2,4)=1
c
      j_ins_1_coll(1,1,3)=3
      j_ins_1_coll(1,1,4)=3
      j_ins_1_coll(1,2,3)=4
      j_ins_1_coll(2,2,3)=4
      j_ins_1_coll(2,2,4)=6
      j_ins_1_coll(2,2,5)=7
      j_ins_1_coll(3,2,3)=4
      j_ins_1_coll(3,2,4)=4
      j_ins_1_coll(3,2,5)=5
c
      j_ins_2_coll(1,1,3)=2
      j_ins_2_coll(1,1,4)=2
      j_ins_2_coll(1,1,5)=1
      j_ins_2_coll(1,2,4)=2
      j_ins_2_coll(1,2,5)=2
      j_ins_2_coll(2,2,3)=1
      j_ins_2_coll(3,2,3)=1
      j_ins_2_coll(3,2,4)=1
c
c The matrix I_AP_IN1(JPROC,JINST,NI) returns the code for the 
c Altarelli-Parisi kernel relevant for the collinear emission
c of parton number NI from the parton number 1. Look at the 
c function AP_KERN for a list of the codes. Look at the function
c X_XP_COLLIN_YP for an explanation of the negative codes.
c The matrix I_AP_IN2 has the same meaning, but the emission is from
c parton number 2
      i_ap_in1(1,1,3)=5
      i_ap_in1(1,1,4)=5
      i_ap_in1(1,2,3)=-5
      i_ap_in1(2,2,3)=-5
      i_ap_in1(2,2,4)=-5
      i_ap_in1(2,2,5)=-5
      i_ap_in1(3,2,3)=-5
      i_ap_in1(3,2,4)=-5
      i_ap_in1(3,2,5)=5
c
      i_ap_in2(1,1,3)=2
      i_ap_in2(1,1,4)=2
      i_ap_in2(1,1,5)=1
      i_ap_in2(1,2,4)=4
      i_ap_in2(1,2,5)=4
      i_ap_in2(2,2,3)=3
      i_ap_in2(3,2,3)=3
      i_ap_in2(3,2,4)=3
c
c The matrix J_PRC_COLL_OUT(IPROC,IINST,NI,NJ) returns the 2 --> 2 process
c number which factorizes when the parton number NI is collinear
c to the parton number NJ. It is here assumed that NJ > NI. This is
c not restrictive, since S(NI,NJ)=S(NJ,NI)
      do i=1,3
        do j=1,2
          do k1=3,4
            do k2=k1+1,5
              j_prc_coll_out(i,j,k1,k2)=-1
            enddo
          enddo
        enddo
      enddo
      j_prc_coll_out(1,1,3,5)=1
      j_prc_coll_out(1,1,4,5)=1
      j_prc_coll_out(1,2,3,4)=1
      j_prc_coll_out(1,2,3,5)=1
      j_prc_coll_out(1,2,4,5)=1
      j_prc_coll_out(2,2,4,5)=1
      j_prc_coll_out(3,2,3,5)=1
      j_prc_coll_out(3,2,4,5)=1
c
c The matrix I_AP_OUT(JPROC,JINST,NI,NJ) returns the code for the 
c Altarelli-Parisi kernel relevant for the collinear splitting
c of partons number NI and number NJ from the parton number S(NI,NJ).
c With our conventions (see the beginning of PHDIJETCRS.FOR) the gluons
c have always parton numbers larger than those for the quarks. Therefore,
c with NJ > NI, the splitting q --> q(NI)g(NJ) is given by P_{gq} (code 3).
c In the main program we allow for NI to be larger than NJ. In that case,
c the splitting has to be described by P_{qq} (code 4). This has been
c properly dealt with in the subroutine XMATEL_FIVE_PART
      i_ap_out(1,1,3,5)=3
      i_ap_out(1,1,4,5)=3
      i_ap_out(1,2,3,4)=3
      i_ap_out(1,2,3,5)=3
      i_ap_out(1,2,4,5)=1
      i_ap_out(2,2,4,5)=2
      i_ap_out(3,2,3,5)=2
      i_ap_out(3,2,4,5)=2
c
c The matrix I_PHYS(N,IPROC,IINST) returns the parton number in the
c physical process defined by (IPROC,IINST) of the parton number N
c in the unphysical process defined by IPROC. If the unphysical
c and the physical processes for (IPROC0,IINST0) are
c
c       0 --> abcde,       b(bar)d(bar) --> aec
c 
c respectively, the matrix I_PHYS would read (the momenta are outgoing)
c
c  I_PHYS(1,IPROC0,IINST0)=3
c  I_PHYS(2,IPROC0,IINST0)=1
c  I_PHYS(3,IPROC0,IINST0)=5
c  I_PHYS(4,IPROC0,IINST0)=2
c  I_PHYS(5,IPROC0,IINST0)=4
c
c The definition of I_PHYS is therefore dependent upon the parton 
c labeling conventions in the functions which return the crossing
c invariant matrix element value
      i_phys(1,1,1)=1
      i_phys(2,1,1)=3
      i_phys(3,1,1)=4
      i_phys(4,1,1)=2
      i_phys(5,1,1)=5
      i_phys(1,1,2)=1
      i_phys(2,1,2)=3
      i_phys(3,1,2)=2
      i_phys(4,1,2)=4
      i_phys(5,1,2)=5
      i_phys(1,2,2)=1
      i_phys(2,2,2)=3
      i_phys(3,2,2)=2
      i_phys(4,2,2)=4
      i_phys(5,2,2)=5
      i_phys(1,3,2)=1
      i_phys(2,3,2)=3
      i_phys(3,3,2)=2
      i_phys(4,3,2)=4
      i_phys(5,3,2)=5
c
c The matrix I0_PHYS has the same meaning of I_PHYS, but is relevant
c for the 2 --> 2 processes
      i0_phys(1,1,1)=1
      i0_phys(2,1,1)=3
      i0_phys(3,1,1)=4
      i0_phys(4,1,1)=2
      i0_phys(1,1,2)=1
      i0_phys(2,1,2)=3
      i0_phys(3,1,2)=2
      i0_phys(4,1,2)=4
      i0_phys(1,4,3)=1
      i0_phys(2,4,3)=3
      i0_phys(3,4,3)=2
      i0_phys(4,4,3)=4
      i0_phys(1,4,4)=1
      i0_phys(2,4,4)=2
      i0_phys(3,4,4)=3
      i0_phys(4,4,4)=4
      i0_phys(1,5,4)=1
      i0_phys(2,5,4)=4
      i0_phys(3,5,4)=2
      i0_phys(4,5,4)=3
      i0_phys(1,5,6)=1
      i0_phys(2,5,6)=4
      i0_phys(3,5,6)=3
      i0_phys(4,5,6)=2
      i0_phys(1,5,7)=1
      i0_phys(2,5,7)=2
      i0_phys(3,5,7)=3
      i0_phys(4,5,7)=4
      i0_phys(1,6,4)=1
      i0_phys(2,6,4)=4
      i0_phys(3,6,4)=2
      i0_phys(4,6,4)=3
      i0_phys(1,6,5)=1
      i0_phys(2,6,5)=2
      i0_phys(3,6,5)=3
      i0_phys(4,6,5)=4
c
c The matrix I_XCOL(N,IPROC,IINST) returns the position of parton
c number N (in the classification scheme used here) in the classification
c scheme in which the crossing-invariant matrix elements are written,
c for the physical process defined by (IPROC,IINST). If the unphysical
c and the physical processes for (IPROC0,IINST0) are (the momenta are
c incoming, following Kunszt-Soper convention for the 2 --> 2 processes)
c
c       0 --> abcd,        da --> b(bar)c(bar)
c 
c respectively, the matrix I_XCOL would read
c
c I_XCOL(1,IPROC0,IINST0)=4
c I_XCOL(2,IPROC0,IINST0)=1
c I_XCOL(3,IPROC0,IINST0)=2
c I_XCOL(4,IPROC0,IINST0)=3
c
c Notice that the definition of I_XCOL is therefore strictly dependent
c upon the conventions adopted for the physical and unphysical
c matrix elements. In particular, it matters whether the momenta
c of the unphysical matrix elements are outgoing or incoming (due to the
c flavour conjugation)
      i_xcol(1,1,1)=1
      i_xcol(2,1,1)=4
      i_xcol(3,1,1)=2
      i_xcol(4,1,1)=3
      i_xcol(1,1,2)=1
      i_xcol(2,1,2)=3
      i_xcol(3,1,2)=2
      i_xcol(4,1,2)=4
c
c The matrix I_SFMAP(IPROC,IINST) returns the value IINST (for physically
c meaningful values of (IPROC,IINST)), except for the case 
c (IPROC,IINST)=(2,2), when is set equal to 3
      do i=1,3
        do j=1,2
          i_sfmap(i,j)=-1
        enddo
      enddo
      i_sfmap(1,1)=1
      i_sfmap(1,2)=2
      i_sfmap(2,2)=3
      i_sfmap(3,2)=2
c
c The matrix IMTT(I,J), with 3<=I<=4 and 3<=J<=5, is defined such
c that {IMTT(I,J_0)| I=3,4, J_0 fixed}={3,4,5}-{J_0}
      imtt(3,3)=4
      imtt(4,3)=5
      imtt(3,4)=3
      imtt(4,4)=5
      imtt(3,5)=3
      imtt(4,5)=4
c
c The matrix IMTT2(I,J), with 3<=I,J<=5, is defined such
c that IMTT2(I_0,J_0)={3,4,5}-{I_0,J_0}
      imtt2(3,4)=5
      imtt2(3,5)=4
      imtt2(4,3)=5
      imtt2(4,5)=3
      imtt2(5,3)=4
      imtt2(5,4)=3
c
c The matrix XCA(I) returns C_A for gluons (I=0) and C_F for quarks (I=1).
c The color factor for photon (I=2) is zero
      xca(0)=xnc
      xca(1)=(xnc**2-1.d0)/(2*xnc)
      xca(2)=0.d0
c
c The matrix XGMM(I) returns the quantities defined in eqs.(3.5) and (3.6)
c The color factor for photon (I=2) is zero
      xgmm(0)=(11.d0*xnc-2.d0*nl)/6.d0
      xgmm(1)=3.d0*xca(1)/2.d0
      xgmm(2)=0.d0
c
c The matrix XGMMPRIME(I) returns the quantities defined in 
c eqs.(A.12) and (A.13). The color factor for photon (I=2) is zero
      xgmmprime(0)=67.d0*xnc/9.d0-2*pi**2*xnc/3.d0-23.d0*nl/18.d0
      xgmmprime(1)=13.d0*xca(1)/2.d0-2*pi**2*xca(1)/3.d0
      xgmmprime(2)=0.d0
c
c The matrix ZEL(I) returns the electric charge of flavour I. The
c numbering convention for the flavours should match the one of
c the structure function package. In our case, 
c 1=UP, 2=DOWN, 3=STRANGE, 4=CHARM, 5=BOTTOM, 6=TOP
      zel(1)=zup
      zel(2)=zdo
      zel(3)=zdo
      zel(4)=zup
      zel(5)=zdo
      zel(6)=zup
c
c The matrix SUMZEL(NL) returns the sum of first NL entries of the
c matrix ZEL; the matrix SUMZEL2 has the same meaning, but returns
c the sum of the charges squared
      sumzel(1)=zup
      sumzel2(1)=zup**2
      do i=2,6
        sumzel(i)=sumzel(i-1)+zel(i)
        sumzel2(i)=sumzel2(i-1)+zel(i)**2
      enddo
c
c end of subroutine setpar
      return
      end


      subroutine x_inst_sum(iinst)
c In this subroutine we fill the matrices IMAX_INST(IPROC) and
c IMIN_INST(IPROC), which will define the sum over the initial 
c states for the 2 --> 3 processes, and the matrices IMAX0_INST(IPROC) 
c and IMIN0_INST(IPROC), which have the same meaning but are relevant 
c for the 2 --> 2 processes
      implicit real*8 (a-h,o-z)
      common/numofproc/num_of_proc,num0_of_proc
      common/irange/imax_inst,imin_inst
      common/i0range/imax0_inst,imin0_inst
      integer num_of_proc(1:3),num0_of_proc(1:3)
      integer imax_inst(1:3),imin_inst(1:3)
      integer imax0_inst(1:3),imin0_inst(1:3)
c
      do i=1,3
        imax_inst(i)=-1
        imin_inst(i)=1
      enddo
      do i=1,3
        imax0_inst(i)=-1
        imin0_inst(i)=1
      enddo
c IINST=0 or IINST<0: All the initial states or all the initial state
c but one are selected. Perform the sum over them all, and reject
c the non-selected one in the program
      if(iinst.le.0)then
        do i=1,3
          imax_inst(i)=num_of_proc(i)
          imin_inst(i)=1
        enddo
        do i=1,3
          imax0_inst(i)=num0_of_proc(i)
          imin0_inst(i)=1
        enddo
c IINST>0: A given initial state was selected.
      elseif(iinst.eq.1)then
        imax_inst(1)=1
        imin_inst(1)=1
c
        imax0_inst(1)=1
        imin0_inst(1)=1
      elseif(iinst.eq.2)then
        imax_inst(1)=2
        imin_inst(1)=2
        imax_inst(2)=1
        imin_inst(2)=1
        imax_inst(3)=1
        imin_inst(3)=1
c
        imax0_inst(1)=2
        imin0_inst(1)=2
        imax0_inst(2)=-1
        imin0_inst(2)=1
        imax0_inst(3)=-1
        imin0_inst(3)=1
      else
        write(6,*)'Error in x_inst_sum: unknown initial state number'
        stop
      endif
      return
      end


      subroutine x_proc_sum(fl_lo,fl_nlo)
c This subroutine checks whether the sum over processes and initial
c states selected by the user is trivial. When the sum is trivial
c (i.e., no contribution) fl_lo (for 2-body processes; fl_nlo is relevant
c for 3-body processes) is set to .false.
      implicit real*8 (a-h,o-z)
      logical fl_lo,fl_nlo
      common/iproc/iproc
      common/jprocess/min_proc,max_proc
      common/irange/imax_inst,imin_inst
      common/i0range/imax0_inst,imin0_inst
      integer imax_inst(1:3),imin_inst(1:3)
      integer imax0_inst(1:3),imin0_inst(1:3)
c
      fl_lo=.false.
      fl_nlo=.false.
      do jproc=min_proc,max_proc
        if(jproc.ne.-iproc)then
          if(imax0_inst(jproc).ge.imin0_inst(jproc))fl_lo=.true.
          if(imax_inst(jproc).ge.imin_inst(jproc))fl_nlo=.true.
        endif
      enddo
      return
      end
c
c End of initialization routines
c


      subroutine HVQWARN(str)
      character *(*) str
      write(*,*) '********** WARNING **********'
      write(*,*) '*********  ',str,'  *********'
      end
