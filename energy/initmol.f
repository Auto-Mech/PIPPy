      subroutine initmol(im)

c Initialize AG #IM.  Use information computed and stored from PREMOL
c to set up each trajectory.

      implicit none
      include 'param.f'
      include 'c_traj.f'
      include 'c_sys.f'
 
      double precision xtot,ytot,ztot,mtot,timpnow,kinow,etmp
      double precision xxm(3,mnat),ppm(3,mnat),mmm(mnat),
     &  phop(mnsurf),dmag,ppj(3,mnat),ek2,ejrot
      integer im,i,j,k,ii
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & pema(mnsurf),gpema(3,mnat,mnsurf),dvec(3,mnat,mnsurf,mnsurf)
      character*2 symb(mnat)

c     local
      integer nqn,whwell,ntry
      double precision molpe,etarget,scalef,sampjj,sampkk,eig(3),etest,
     & dum3(3)
c      logical letot

c HARD-CODED --- changed to input flag letot(im)
c      letot = .true. ! fixed total energy
c      letot = .false. ! fixed vib+K2 energy

      write(6,*)"Generating initial conditions for AG ",im
      write(6,*)"-------------------------------------------"

c     store coords for this molecule in their own array for easier manipulation
      do i=1,natom(im)
       ii=i+iatom(im)
       mmm(i) = mm(ii)
       symb(i) = symbol(ii)
       do j=1,3
        xxm(j,i) = 0.d0
        ppm(j,i) = 0.d0
       enddo
      enddo
      molpe=0.d0
      ntry=0

c GENERATE INITIAL COORDINATES
 111  continue
      if (initx(im).eq.0) then
c       they were read in readin and stored in xx0
        write(6,*)"INITx = 0:  Initial coordinates are ",
     &      "specified by the user"
        do i=1,natom(im)
        do j=1,3
         ii=i+iatom(im)
         xxm(j,i) = xx0(j,ii)
        enddo
        enddo
      else
        write(6,*)"ERROR:  INITx = ",initx(im)," is not an option"
        stop
      endif

c put center of mass at the origin
!      write(6,100)"Placing CoM at origin"
 100  format(1x,a,3f12.5,1x,a)
!      write(6,*)
!      xtot=0.d0
!      ytot=0.d0
!      ztot=0.d0
!      mtot=0.d0
!      do i=1,natom(im)
!          xtot = xtot + mmm(i)*xxm(1,i)
!          ytot = ytot + mmm(i)*xxm(2,i)
!          ztot = ztot + mmm(i)*xxm(3,i)
!          mtot = mtot + mmm(i)
!      enddo
!      do i=1,natom(im)
!          xxm(1,i) = xxm(1,i) - xtot/mtot
!          xxm(2,i) = xxm(2,i) - ytot/mtot
!          xxm(3,i) = xxm(3,i) - ztot/mtot
!      enddo

c spin randomly if necessary
c HACK AJ
!      if (iorient.eq.1) call spin(xxm,ppm,natom(im))
!      if (iorient.eq.2) call spin(xxm,ppm,natom(im))
c HACK AJ

c print coordinates
c      if (.false.) then ! Minimal printing
      write(6,*)"Initial coordinates"
      write(6,105)"    # symb","x (A)","y (A)","z (A)"
 105  format(a10,3a12)
      do i=1,natom(im)
          ii=iatom(im)+i
          write(6,101)ii,symbol(ii),xxm(1,i)*autoang,
     &                xxm(2,i)*autoang,xxm(3,i)*autoang
 101      format(i5,a5,3f12.5)
      enddo
      write(6,*)
c      endif

c if we haven't calcualted the PE yet 
c      if (molpe.eq.0.d0) then
!      print *,"ldofrag = ",ldofrag
!      if (ldofrag) then
! CHECK HERE, FAILING FOR AG 1 H2O+H DM
      call getpem(xxm,natom(im),pema,pemd,gpema,gpemd,dvec,symb)
!        if (repflag.eq.1) then
!          molpe = pemd(nsurf,nsurf)
!        else
      print *,"PEM adiabatic = ",pema
      print *,"nsurf = ",nsurf
      molpe = pema(nsurf)
!        endif
c       correct for this atom groups zero of energy
      print *,"ezero = ",ezero
      print *,"ezeroim = ",ezeroim(im)
      print *,"im = ",im
      molpe = molpe + ezero - ezeroim(im)
!      else
!        write(6,*)"Initial potential energy evalulation skipped"
!        molpe = 0.d0
!      endif
c      endif

c GENERATE INITIAL MOMENTA
      if (initp(im).eq.0) then
c     random thermal sample
        write(6,*)"INITp = 0:  Generating initial velocities from a"
     & ," random thermal distribution"
!        call rantherm(ppm,mmm,natom(im),temp0im(im))
        write(6,106)"Desired temp    = ",temp0im(im)," K"
      elseif (initp(im).eq.1) then
       write(6,*)"INITp = 1:  Setting initial momenta to zero"
       do i=1,natom(im)
         do j=1,3
         ppm(j,i)=0.d0
         enddo
       enddo
!      elseif (initp(im).eq.2) then
!       write(6,*)"INITp = 2:  Initial momenta are specified by the user"
!       do i=1,natom(im)
!         do j=1,3
!         ppm(j,i)=pp0(j,i)
!         enddo
!       enddo
!      elseif (initp(im).eq.-1) then
!        write(6,*)"INITp = -1:  Inititial momenta selected along with",
!     & " initial coordinates"
!        write(6,*)
      else
        write(6,*)"ERROR:  INITp = ",initp(im)," is not an option"
        stop
 106  format(1x,a,f12.5,1x,a)
      write(6,*)
      endif

      call gettemp(ppm,mmm,natom(im),temp,ke)
      write(6,106)"Calculated temp = ",temp," K"
      write(6,106)"Total KE        = ",ke*autoev," eV"
      write(6,*)

c MANIPULATE COORDINATES AND MOMENTA
c     remove overall momenta
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,natom(im)
        xtot=xtot+ppm(1,i)
        ytot=ytot+ppm(2,i)
        ztot=ztot+ppm(3,i)
        mtot = mtot + mmm(i)
      enddo
      do i=1,natom(im)
        ppm(1,i)=ppm(1,i)-xtot*mmm(i)/mtot
        ppm(2,i)=ppm(2,i)-ytot*mmm(i)/mtot
        ppm(3,i)=ppm(3,i)-ztot*mmm(i)/mtot
      enddo
      call gettemp(ppm,mmm,natom(im),temp,ke)
      write(6,*)"Removed CoM motion"
      write(6,106)"Calculated temp = ",temp," K"
      write(6,106)"Total KE        = ",ke*autoev," eV"
      write(6,*)

c     calculate angular motion
      if (natom(im).gt.2.and.initx(im).ne.3) then
c     don't do for diatoms or for special atom-diatom initial conditions
      call ange(xxm,ppm,mmm,natom(im),eig,bigj,bigjtot,erot,erottot)
      write(6,106)"Total angular momentum      ",bigjtot," au"
      write(6,107)"Angular momentum components ",bigj," au"
      write(6,106)"Total rotational energy     ",erottot*autoev," eV"
      write(6,107)"Rotational energy components",erot(1)*autoev,
     &   erot(2)*autoev,erot(3)*autoev," eV"
      write(6,*)
 107  format(1x,a,3f12.5,1x,a)

      if (initj(im).ne.1) then
c Comment \/ this to skip removal of angular momentum
c      call noang(xxm,ppm,mmm,natom(im))
c      write(6,*)"Removed overall angular motion"
c      call ange(xxm,ppm,mmm,natom(im),eig,bigj,bigjtot,erot,erottot)
c      call gettemp(ppm,mmm,natom(im),temp,ke)
c      write(6,106)"Total angular momentum      ",bigjtot," au"
c      write(6,107)"Angular momentum components ",bigj," au"
c      write(6,106)"Total rotational energy     ",erottot*autoev," eV"
c      write(6,107)"Rotational energy components",erot(1)*autoev,
c     &   erot(2)*autoev,erot(3)*autoev," eV"
c Comment ^ this to skip removal of angular momentum
      call gettemp(ppm,mmm,natom(im),temp,ke)
      write(6,*)
      write(6,106)"Calculated temp = ",temp," K"
      write(6,106)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      endif

      if (initp(im).eq.0.and.escale0im(im).gt.0.d0) then
c       rescale momenta
       do i=1,natom(im)
         do j=1,3
         ppm(j,i)=ppm(j,i)*dsqrt(escale0im(im)/ke)
         enddo
       enddo
      call gettemp(ppm,mmm,natom(im),temp,ke)
      write(6,*)"Rescaling KE to ",escale0im(im)*autoev," eV"
      write(6,106)"Calculated temp = ",temp," K"
      write(6,106)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      endif

      if (.false.) then ! Minimal printing
      write(6,*)"Initial momenta for this AG"
      write(6,105)"    # symb","x (au)","y (au)","z (au)"
      do i=1,natom(im)
      ii=i+iatom(im)
          write(6,102)ii,symbol(ii),ppm(1,i),ppm(2,i),ppm(3,i)
      enddo
      write(6,*)
      endif

 102      format(i6,a5,3e12.5)

      else

      bigjtot = 0.d0
      erottot = 0.d0

      endif

c print coordinates
      if (.false.) then ! Minimal printing
      write(6,*)"Initial coordinates"
      write(6,105)"    # symb","x (A)","y (A)","z (A)"
      do i=1,natom(im)
          ii=iatom(im)+i
          write(6,101)ii,symbol(ii),xxm(1,i)*autoang,
     &                xxm(2,i)*autoang,xxm(3,i)*autoang
      enddo
      write(6,*)
      endif

c WRITE INITIAL ENERGY AND ANGULAR MOMENTUM
      write(6,*)"Info for this AG only"
!      if (ldofrag) then
      call getpem(xxm,natom(im),pema,pemd,gpema,gpemd,dvec,symb)
        if (repflag.eq.1) then
          pe = pemd(nsurf,nsurf)
        else
          pe = pema(nsurf)
        endif
c       correct for this atom groups zero of energy
        pe = pe + ezero - ezeroim(im)
!      else
!        write(6,*)"Initial potential energy evalulation skipped"
!        pe = 0.d0
!      endif
      peagi(im)=pe
      call gettemp(ppm,mmm,natom(im),temp,ke)
!      call ange(xxm,ppm,mmm,natom(im),eig,bigj,bigjtot,erot,erottot)
      te = pe+ke
      write(6,106)"Initial potential energy   = ",
     &   pe*autoev," eV"
!      write(6,106)"Initial kinetic energy     = ",
!     &   ke*autoev," eV"
!      write(6,106)"Initial rotational energy  = ",
!     &   erottot*autoev," eV"
!      write(6,106)"Initial vibrational energy = ",
!     &   (ke-erottot)*autoev," eV"
!      write(6,106)"Initial total energy       = ",
!     &   te*autoev," eV"
!      write(6,106)"Initial temp               = ",
!     &   temp," K"
!      write(6,106)"Initial angular momentum   = ",
!     &   bigjtot," au"
      write(6,*)

c save AG info
!      erottotagi(im)=erottot
!      evibagi(im)=ke-erottot
!      eintagi(im)=ke
!      bigjtotagi(im)=bigjtot
!      do i=1,3
!      bigjagi(i,im)=bigj(i)
!      erotagi(i,im)=erot(i)
!      enddo

c write manipulated data back to xx and pp arrays
      do i=1,natom(im)
        ii=i+iatom(im)
        do j=1,3
          xx(j,ii) = xxm(j,i)
          pp(j,ii) = ppm(j,i)
        enddo
      enddo

      print *,"xx = "
      print *,((xx(i,j),i=1,natom(im)),j=1,3)!,((pp(i,j),i=1,3),j=1,3)

c temp
      rtrans=0.d0

      return

      end
