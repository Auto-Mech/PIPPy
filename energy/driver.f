      subroutine driver

c Computes a full trajectory from start to finish.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'
c hack
      common/tmp/tmpprint
      double precision tmpprint(50)
      double precision taup(3),taup2(3)
      common/tauprint/taup,taup2
      double precision h12x,h12y
      common/pesh12/h12x,h12y

c     number of bins for binning the radial distribution function
      integer nbinrad,im
      parameter(nbinrad=60)

      integer i,j,k,ii,iturn,i1,i2,ithistraj,frusflag(mnsurf),
     &  newsurf,iramp,ncross
      double precision rminmax,ravg,ce,hrad,dot,peim,
     &  ti1,ti2,tu1,tu2
      double precision rr(3),rhor(mnsurf,mnsurf),eig(3),dum3(3),
     &  rhoi(mnsurf,mnsurf),step,timeprev,phop(mnsurf),tmp1,tmp2,tmp3,
     &  dmag,dmag1,dmag2,gvmag,rrr,rmax(mnoutcome),rmin(mnoutcome),
     &  dvec(3,mnat,mnsurf,mnsurf),gpem(3,mnat,mnsurf),
     &  arij(mnat,mnat),arij2(mnat,mnat),lind,
     &  raddist(0:nbinrad+1),egap,rtmp,rcom,tmp4,ecom,
     &  lastgap0,lastgap,plz,plz0,elz,plzx,elz0,elzx,lastcross,
     &  hlz,hlz0,hlzx,glz,glzx,glz0,r1,r10,r1x,r2,r20,r2x, 
     7  a123,a1230,a123x,r3,vlz,vlz0,vlzx,pairy(2),pairy0(2),pairyx(2)
      logical lhop

c used for FSTU
      double precision tu_xx(3,mnat,mnsurf),tu_pp(3,mnat,mnsurf),
     & tu_cre(mnsurf,mnsurf),tu_cim(mnsurf,mnsurf),tu_hstep(mnsurf),
     & tu_time(mnsurf),tu_phase(mnsurf,mnsurf),tu_maxt(mnsurf),
     & deltforw,deltback,tu_maxtime
      integer tu_istep(mnsurf)
      logical lfrust

c Andersen thermostat
      logical lhit

c stochastic decoherence
      integer nhop
      double precision stodecotime,stodecotau,electau

      integer wellindex1,wellindex2,n1,n2,ijk
      double precision xwell,rvdw,evdw
      logical lzflag
      lzflag=.false.
      if (tflag(4).lt.0) lzflag=.true.
c      rvdw=0.d0
      evdw=1000.d0

c initialize for this trajectory
      lastgap0 =0.d0 ! TEMP AJ
      lastgap =0.d0 ! TEMP AJ
      lastcross =0.d0 ! TEMP AJ
      ncross =0 ! TEMP AJ
      nhop = 0
      time = 0.d0
      stodecotime = -1.d0
      stodecotau = 4.d0/autofs
      iturn = 0
      dmag1 = 0.d0
      dmag2 = 0.d0
      nsurft = 1 ! TEMP DM
      nsurf0 = 1 ! TEMP DM
      nsurf = 1 ! TEMP DM
      repflag = 0 ! TEMP DM
      do i=1,nsurft
        phase(i) = 0.d0
        tu_time(i) = -1.d0
      enddo
      lfrust = .false.
      do i=1,mnat
      do j=1,mnat
        arij(i,j) = 0.d0
        arij2(i,j) = 0.d0
      enddo
      enddo
      atemp = 0.d0
      atemp2 = 0.d0
      stemp = 0.d0
      iramp = 0
      do i=0,nbinrad+1
        raddist(i) = 0.d0
      enddo
      outcome = -1   ! OUTCOME should get reassigned when the traj finishes. A final value of -1 indicates a failure of some kind.
      istep = 0
      istepw = 0
      step = 0.d0

c calculate and save initial values of things for later analysis
c     initial coordinates and momenta
      do i=1,nat
      do j=1,3
         xxi(j,i)=xx(j,i)
         ppi(j,i)=pp(j,i)
      enddo
      enddo
!      call gettemp(pp,mm,nat,tempi,kei)
!      call timing(ti1)
!      print *,"xx = ",xx
!      print *,"pp = ",pp
      print *,"nsurf = ",nsurf
      call getgrad(xx,pp,nsurf,pei,cre,cim,gv,gcre,gcim,nat,
     & phop,dmag,dvec,pem,gpem,phase)
!      call timing(ti2)
!      write(6,'(" CPU time in initial force calc is ",f20.10," s",/)')
!     &    (ti2-ti1)
!      call getplz(pp,mm,gpem,nat,plz,elz,hlz,glz,pairy,nsurf)  ! TEMP AJ

!      tei=kei+pei
!      call ange(xx,pp,mm,nat,eig,bigji,bigjtoti,eroti,erottoti)
!      write(6,106)"Initial overall ang momentum = ",bigjtoti," au"
!      write(6,106)"Initial overall rot energy   = ",erottoti*autoev,
!     & " eV"
!      write(6,106)"Initial kinetic energy       = ",kei*autoev," eV"
      write(6,106)"Initial potential energy     = ",pei*autoev," eV"
!      write(6,106)"Initial total energy         = ",tei*autoev," eV"
!      write(6,106)"Initial temp                 = ",tempi," K"
 106  format(1x,a,f15.5,1x,a)
 107  format(1x,a,3f12.5,1x,a)
      write(6,*)

c write initial coordinates and momenta to output
      if (lwrite(10)) then
      open(10)
      write(10,*)ithistraj,pei*autoev
      do i=1,nat
        write(10,1010)symbol(i),mm(i)/amutoau,(xx(j,i)*autoang,j=1,3)
      enddo
      endif
!      if (lwrite(11)) then
!        open(11)
!        write(11,*)ithistraj,kei*autoev
!        do i=1,nat
!          write(11,1010)symbol(i),mm(i)/amutoau,(pp(j,i),j=1,3)
!        enddo
!      endif

 1010 format(5x,a2,4f20.8)

c trajectory output
      if (termflag.eq.2) then
      write(6,'(20a15)')"          step",
     &                  "      time(fs)",
     &                  "        PE(eV)",
     &                  " |gradV|(eV/A)"
      else
      if (nsurft.eq.1.and.lwell.le.0) then
      write(6,'(20a15)')"           step",
     &                  "       time(fs)",
     &                  "           T(K)",
     &                  "         KE(eV)",
     &                  "         PE(eV)",
     &                  "         TE(eV)",
     &                  "      Lindemann",
     &                  "         <T(K)>",
     &                  "   st.dev(T(K))"
      else if (nsurft.eq.1.and.lwell.gt.0) then
      write(6,'(20a15)')"           step",
     &                  "          stepw",
     &                  "       time(fs)",
     &                  "           T(K)",
     &                  "         KE(eV)",
     &                  "         PE(eV)",
     &                  "         TE(eV)",
     &                  "      Lindemann",
     &                  "         <T(K)>",
     &                  "   st.dev(T(K))"
      else
      write(6,'(20a15)')"           step",
     &                  "       time(fs)",
     &                  "           T(K)",
     &                  "         KE(eV)",
     &                  "         PE(eV)",
     &                  "         TE(eV)",
     &                  "          Probs"
      endif
      endif

 999  continue

!      if (lwrite(16)) 
!     & write(16,'(2f13.5,2e18.5,3e18.5)')time*autofs,1.d4,
!     & phop(1)*step,rhor(1,1),
!     & (taup2(i)*autofs,i=1,3)

 1001 return
      end

