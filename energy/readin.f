      subroutine readin

c This subroutine reads from the standard input.  This is the only 
c subroutine that reads data.  Data are stored for later use in
c common blocks.

      implicit none

      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
      include 'c_traj.f'
      include 'c_output.f'

c local variables
      integer i,ii,im,k,ntmp,itmp(99),nrecl
      double precision rdum,murel,prel,dum
      character*2 chardum

c read potential interface flag
      write(6,*)"Potential interface"
      write(6,*)"-------------------"
      read(5,*)potflag
      write(6,*)"POTFLAG = ",potflag
      if (potflag.eq.0) then
      write(6,*)"Using the POTLIB MM-HO-1 interface"
      elseif (potflag.eq.2) then
      write(6,*)"Using the POTLIB MM-HE-1 interface"
      elseif (potflag.eq.1) then
      write(6,*)"Using the POTLIB 3-2V interface"
      elseif (potflag.eq.3) then
      write(6,*)"Using the multistate diabatic MM-HE-1 interface"
      elseif (potflag.eq.4) then
      write(6,*)"Using the multistate adiabatic MM-HE-1 interface"
c     ADD YOUR OWN INTERFACE HERE, USING POTFLAG = -1
c     YOU NEED ONLY MODIFY GETPEM
c      elseif (potflag.eq.-1) then
c      write(6,*)"Using a user-designed interface"
      else
      write(6,*)"POTFLAG = ",potflag," is not supported"
      stop
      endif
      write(6,*)

      call prepot

c LDOFRAG
!      ldofrag = .false.
      ldofrag = .true.

      nsurft = 1 ! TEMP DM
      nsurf0 = 1 ! TEMP DM
      nsurf = 1 ! TEMP DM

c initial surface, number of coupled surfaces, and electronic representation flag
!      IF (my_rank.eq.0) THEN
!      read(5,*)nsurf0,nsurft,methflag,repflag
!      ENDIF

c read number of molecules (atom groups)
      write(6,*)"Atom group (AG) specification"
      write(6,*)"-----------------------------"
      read(5,*)nmol,ezero
      write(6,*)"This is a ",nmol," AG calculation"
      write(6,*)"The zero of energy for the supermolecule is ",
     & ezero," eV"
      ezero=ezero/autoev
      write(6,*)
      if (nmol.gt.mnmol) then
        write(6,*)"ERROR:  nmol > mnmol, ",nmol,">",mnmol
        stop
      endif
      if (nmol.le.0) then
        write(6,*)"ERROR:  nmol must be > 0, nmol=",nmol
        stop
      endif

c read info for each molecule
      nat = 0
      DO im=1,nmol

      write(6,*)"AG ",im
      write(6,*)"-------"

      read(5,*)natom(im),initx(im),initp(im),initj(im),ezeroim(im)
      iatom(im)=nat
      nat=nat+natom(im)
      write(6,*)"This AG has ",natom(im)," atoms (",
     &  "atoms ",iatom(im)+1," to ",nat,")"
      if (natom(im).le.0) then
        write(6,*)"ERROR: cant have natom < 1, natom(",im,")=",natom(im)
        stop
      endif
      if (nat.gt.mnat) then
        write(6,*)"ERROR:  nat > mnat, ",nat,">",mnat
        stop
      endif
      write(6,*)"Zero of energy for this AG is ",ezeroim(im)," eV"
      ezeroim(im)=ezeroim(im)/autoev

      write(6,*)"Initial conditions prepared using INITx = ",initx(im)
!     &   ," INITp = ",initp(im)," and INITj = ",initj(im)

c INITx controls the initial coordinates
c read different info depending on the value of INITx(IM)
      if (initx(im).eq.0) then
c       crystal structure
        write(6,*)"INITx = 0:  Reading in initial coordinates:"
        write(6,102)"    # symb","mass (amu)","x (A)","y (A)","z (A)"
 102    format(a10,4a12)
        mmag(im)=0.d0
        do i=1,natom(im)
           ii=iatom(im)+i
c          read in atomic symbol, coords, and mass
           read(5,*)symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
           write(6,103)ii,symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
 103       format(i5,a5,4f12.5)
c          convert from amu to au
           mm(ii)=mm(ii)*amutoau
c          compute total mass for this AG
           mmag(im)=mmag(im)+mm(ii)
c          convert from A to bohr
           do k=1,3
              xx0(k,ii)=xx0(k,ii)/autoang
           enddo
        enddo
        write(6,*)
!      elseif (initx(im).eq.3) then
!c       atom-diatom initial conditions
!        write(6,*)"INITx = 3:  Atom-diatom initial conditions"
!        if (natom(im).ne.3) then
!          write(6,*)"INITx = 3 must be used with a triatomic molecule"
!          stop
!        endif
!        if (nmol.ne.1) then
!          write(6,*)"INITx = 3 may not be used with multiple AGs"
!          stop
!        endif
!        write(6,102)"    # symb","mass (amu)"
!        mmag(im)=0.d0
!        do i=1,natom(im)
!           ii=iatom(im)+i
!c          read in atomic symbol, and mass
!           read(5,*)symbol(ii),mm(ii)
!           write(6,103)ii,symbol(ii),mm(ii)
!c          convert from amu to au
!           mm(ii)=mm(ii)*amutoau
!c          compute total mass for this AG
!           mmag(im)=mmag(im)+mm(ii)
!        enddo
!        read(5,*)escatad,vvad,jjad,rrad,arrad
!        write(6,*)"Total energy = ",escatad," eV"
!!        write(6,*)"Initial v = ",vvad
!!        write(6,*)"Initial j = ",jjad
!        write(6,*)"Initial distance = ",rrad
!        write(6,*)"Initial arrangement = ",arrad
!!        write(6,*)"For INITx = 3, INITp is not used and is set to -1"
!        initp(im) = -1
!        escatad=escatad/autoev
!        rrad=rrad/autoang
!        write(6,*)
!      elseif (initx(im).eq.4) then
!c       NOT IMPLEMENTED. DO NOT USE.
!c       diatom initial conditions
!        write(6,*)"INITx = 4:  Diatom initial conditions"
!        if (natom(im).ne.2) then
!          write(6,*)"INITx = 4 must be used with a diatomic molecule"
!          stop
!        endif
!        read(5,*)vvdi(im),jjdi(im),rmindi(im)
!        if (vvdi(im).ge.0.d0) then
!        write(6,*)"Initial v = ",vvdi(im)
!        else
!        write(6,*)"Initial v selected from thermal distribution, T = "
!     &        ,-vvdi(im)," K"
!        endif
!        if (jjdi(im).ge.0.d0) then
!        write(6,*)"Initial j = ",jjdi(im)
!        else
!        write(6,*)"Initial j selected from thermal distribution, T = "
!     &        ,-jjdi(im)," K"
!        endif
!        write(6,*)"Estimated R = ",rmindi(im)
!        rmindi(im)=rmindi(im)/autoang
!        write(6,*)
!        write(6,102)"    # symb","mass (amu)"
!        mmag(im)=0.d0
!        do i=1,natom(im)
!           ii=iatom(im)+i
!c          read in atomic symbol, and mass
!           read(5,*)symbol(ii),mm(ii)
!           write(6,103)ii,symbol(ii),mm(ii)
!c          convert from amu to au
!           mm(ii)=mm(ii)*amutoau
!c          compute total mass for this AG
!           mmag(im)=mmag(im)+mm(ii)
!        enddo
!      elseif (initx(im).eq.6) then
!        write(6,*)"INITx = 6: Read in initial coordinates and momenta ",
!     &   "from a separate files"
!        read(5,*)samptot(im),lbinsamp(im),sampfilexx(im),sampfilepp(im),
!     &   lems(im)
!c     &   sampwell(im),lems(im)
!          sampwell(im)=1
!        write(6,*)"Sampling from ",samptot(im)," data contained in ",
!     &    "the XYZ file ",sampfilexx(im)," and the PxPyPz file ",
!     &    sampfilepp(im)
!        if (sampwell(im).eq.2) then 
!          read(5,*)samptot2(im),lbinsamp2(im),sampfilexx2(im),
!     &     sampfilepp2(im),relwell2(im)
!          relwell1(im)=1.d0-relwell2(im)
!        endif
!        if (lems(im)) then
!          read(5,*)emicr(im),ninc(im),nbrea(im),nemstot(im)
!          write(6,*)"Efficient Microcanonical Sampling (EMS) ",
!     &     "will be used to refine the resulting coordinates ",
!     &     "and momenta to generate a microcanonical ensemble ",
!     &     "with Etot = ",emicr(im)," eV"
!          emicr(im)=emicr(im)/autoev
!          emsw(im)=0.d0
!        endif
!        write(6,*)
!        write(6,102)"    # symb","mass (amu)"
!        mmag(im)=0.d0
!        do i=1,natom(im)
!           ii=iatom(im)+i
!c          read in atomic symbol, coords, and mass
!           read(5,*)symbol(ii),mm(ii)
!           write(6,103)ii,symbol(ii),mm(ii)
!c          convert from amu to au
!           mm(ii)=mm(ii)*amutoau
!c          compute total mass for this AG
!           mmag(im)=mmag(im)+mm(ii)
!c          convert from A to bohr
!        enddo
!!        write(6,*)
!!        write(6,*)"For INITx = 6, INITp is not used and is set to -1"
!!        write(6,*)
!        initp(im) = -1
      else
        write(6,*)"ERROR:  initial conditions flag INITx = ",initx(im),
     &      " is not supported"
        stop
      endif

c INITp controls the initial momenta
!      if (initp(im).eq.0) then
!c       initial target temperature
!        read(5,*)temp0im(im),escale0im(im)
!        write(6,101)"INITp = 0:  Target temp for random thermal ",
!     &    "dist is ",temp0im(im),"K"
!        if (escale0im(im).gt.0.d0) then
!          write(6,101)"Resulting momenta scaled such that the ",
!     &     "kinetic energy = ",escale0im(im)," eV"
!          escale0im(im)=escale0im(im)/autoev
!        else
!          write(6,*)"Resulting momenta are not scaled"
!        endif
!        write(6,*)
! 101    format(1x,a,a,f12.5,1x,a)
!      elseif (initp(im).eq.1) then
!        write(6,*)"INITp = 1:  Initial momenta set to zero"
!        write(6,*)
!      elseif (initp(im).eq.2) then
!        write(6,*)"INITp = 2:  Reading in initial momenta"
!        write(6,*)
!        write(6,105)"x (au)","y (au)","z (au)"
! 105    format(3a12)
!        do i=1,natom(im)
!           ii=iatom(im)+i
!c          read in atomic momenta
!           read(5,*)chardum,dum,(pp0(k,ii),k=1,3)
!           write(6,106)ii,(pp0(k,ii),k=1,3)
! 106       format(i5,4f12.5)
!        enddo
!      elseif (initp(im).eq.-1) then
!        write(6,*)"INITp = -1:  Initial momenta determined by INITx"
!        write(6,*)
!      else
!        write(6,*)"ERROR:  initial conditions flag = ",initp(im),
!     &      " is not supported"
!        stop
!      endif
!
!c INITj
!      if (initj(im).eq.0) then
!        write(6,*)"INITj = 0: Total angular momentum will ",
!     &  "not be adjusted"
!        write(6,*)
!      elseif (initj(im).eq.1) then
!        write(6,*)"INITj = 1: Thermal state selected J"
!        read(5,*)samptarg(im),letot(im),sampjmin(im),sampjmax(im),
!     &   sampjtemp1(im),sampjtemp2(im),sampjbrot1(im),sampjbrot2(im)
!        if (letot(im)) 
!     &      write(6,*)"Initial total energy is fixed"
!        if (.not.letot(im)) 
!     &      write(6,*)"Initial Evib+K2 energy is fixed"
!        if (samptarg(im).gt.0.d0) then
!          write(6,*)"Sampled vibrational momenta will be scaled to an ",
!     &     "energy of ", samptarg(im)," eV, independent of J"
!          samptarg(im)=samptarg(im)/autoev
!        elseif (samptarg(im).eq.0.d0) then
!          write(6,*)"Sampled momenta will not be scaled"
!        elseif (samptarg(im).lt.0.d0) then
!          write(6,*)"Sampled vibrational momenta will be scaled ",
!     &   "according to the rule: Evib(J)=(A+B*J+C*J^2)*D"
!          read(5,*)(ejsc(k,im),k=1,4)
!          write(6,*)"A=",ejsc(1,im)," eV"
!          write(6,*)"B=",ejsc(2,im)," eV"
!          write(6,*)"C=",ejsc(3,im)," eV"
!          write(6,*)"D=",ejsc(4,im)
!          if (sampjmin(im).lt.0.d0) then
!            write(6,*)"ERROR: Jmin > 0 is required with EJSCALE.f"
!            stop
!          endif
!        endif
!        write(6,*)
!        if (sampjmin(im).ge.0.d0.and.sampjtemp1(im).lt.0.d0) then
!          write(6,*)"Total angular momentum assigned evenly from ",
!     &      sampjmin(im)," to ",sampjmax(im)," hbar"
!       elseif (sampjmin(im).ge.0.d0.and.sampjtemp1(im).ge.0.d0) then
!         write(6,*)"Total angular momentum assigned from a ",
!     &     sampjtemp1(im)," to ",sampjtemp2(im),
!     &     " K thermal distribution "
!         write(6,*)" and from",
!     &     sampjmin(im)," to ",sampjmax(im)," hbar with a ",
!     &     "rotational constants of ",sampjbrot1(im)," cm-1",
!     &     " (x2) and ",sampjbrot2(im)," cm-1 (unique)"
!          sampjbrot1(im)=sampjbrot1(im)/autocmi
!          sampjbrot2(im)=sampjbrot2(im)/autocmi
!        else
!          write(6,*)"Total angular momentum will not be adjusted"
!        endif
!        write(6,*)
!      else
!        write(6,*)"ERROR:  initial conditions flag INITj= ",initj(im),
!     &      " is not supported"
!        stop
!      endif

      ENDDO

c AG orientational information
!      write(6,*)"AG orientational prescription"
!      write(6,*)"-----------------------------"
!      if (nmol.eq.1) then
!        write(6,*)"Single AG calculation, no orientational",
!     &                   " info required"
!        do k=1,3
!          comxx(k,1) = 0.d0
!          compp(k,1) = 0.d0
!        enddo
!      else if (nmol.gt.1) then
!        read(5,*)iorient,ldofrag
!        if (ldofrag) then
!        write(6,*)"Fragment energies will be computed."
!        else
!        write(6,*)"Fragment energies will not be computed."
!        endif
!        write(6,*)
!        if (iorient.eq.0) then
!c         read a list of x,y,z and Px,Py,Pz
!          do i=1,nmol
!          read(5,*)(comxx(k,i),k=1,3)
!          read(5,*)(compp(k,i),k=1,3)
!          write(6,*)"For AG ",i
!          write(6,*)"Center of mass coords (x,y,z) (A)  = ",
!     &        (comxx(k,i),k=1,3)
!          write(6,*)"Center of mass motion (x,y,z) (au) = ",
!     &        (compp(k,i),k=1,3)
!          do k=1,3
!            comxx(k,i) = comxx(k,i)/autoang
!          enddo
!          enddo
!        elseif (iorient.eq.1) then
!          if (nmol.ne.2) then
!            write(6,*)"IORIENT = 1 works with NMOL = 2 only. NMOL = ",
!     &       nmol
!            stop
!          endif
!          write(6,*)"Quasiclassical bimolecular collisions"
!c         quasiclassical state-selected orientation
!          read(5,*)rel0qc,tempqc,bminqc,bmaxqc
!          write(6,*)"Initial separation (A) = ",rel0qc
!          write(6,*)"Impact parameter range (A) = ",bminqc," to ",
!     & dabs(bmaxqc)
!          if(bmaxqc.gt.0.d0)write(6,*)"Sampled uniformly from 0-bmax"
!          if(bmaxqc.lt.0.d0)write(6,*)"Sampled uniformly from 0-bmax^2"
!          if (tempqc.gt.0.d0) then
!            write(6,*)"Relative temperature (K) = ",tempqc
!          else if (tempqc.eq.0.d0) then
!            write(6,*)"Relative temperature set to sampled ",
!     &      "J temperature for AG 1"
!          else if (tempqc.lt.0.d0) then
!            write(6,*)"Constant initial relative energy = ",
!     &                 dabs(tempqc)," eV"
!            tempqc=tempqc/autoev
!          endif
!          if (dabs(bmaxqc).gt.rel0qc) then
!             write(6,*)"Maximum impact parameter (",bmaxqc," A) cannot",
!     &          " be greater than initial separation (",rel0qc," A)"
!             stop
!          endif
!          bminqc=bminqc/autoang
!          bmaxqc=bmaxqc/autoang
!          rel0qc=rel0qc/autoang
!        else
!          write(6,*)"IORIENT = ",iorient," is not supported"
!          stop
!        endif
!      else
!        write(6,*)"Cant have AG < 0"
!        stop
!      endif
!      write(6,*)

c output flags
      write(6,*)"Output options"
      write(6,*)"--------------"
      read(5,*)ii,(itmp(k),k=1,ii)
      if (ii.eq.0) then 
      write(6,*)"Writing to all output files"
      do i=1,99
        lwrite(i) = .true.
      enddo
      else
      write(6,*)"Writing to units ",(itmp(k),k=1,ii)
      do i=1,99
        lwrite(i) = .false.
      enddo
      do i=1,ii
        lwrite(itmp(i)) = .true.
        if (itmp(i).eq.77) 
     & write(6,*)"Unit 77 chosen: Writing special named output files"
      enddo
      endif
      write(6,*)

c print headers
      if (lwrite(80)) then
        write(80,*)"[Molden Format]"
        write(80,*)"[GEOMETRIES] XYZ"
      endif
      if (lwrite(81)) then
        write(81,*)"title"
        write(81,*)"momenta XYZ au "
      endif
      nrecl=(3+5*nat)*10
      if (lwrite(82)) open(82,form='unformatted',
     &   status='unknown',access='direct',RECL=nrecl)
      if (lwrite(83)) open(83,form='unformatted',
     &   status='unknown',access='direct',RECL=nrecl)
      
      do im=1,nmol
      if (lems(im).or.lwell.gt.0) then
        if (lwell.gt.0) then
          nrecl=(3+5*nat)*10
        else
          nrecl=(3+5*natom(im))*10
        endif
        if (lwrite(86)) then
          write(86,*)"[Molden Format]"
          write(86,*)"[GEOMETRIES] XYZ"
        endif
        if (lwrite(87)) then
          write(87,*)"title"
          write(87,*)"momenta XYZ au "
        endif
        if (lwrite(88)) open(88,form='unformatted',
     &   status='unknown',access='direct',RECL=nrecl)
        if (lwrite(89)) open(89,form='unformatted',
     &     status='unknown',access='direct',RECL=nrecl)
      endif
      enddo

c check for bad input combinations
      if (potflag.eq.1.and.nat.ne.3) then
        write(6,*)"POTFLAG = 1 must be used with 3 atoms"
        stop 
      endif

      return
      end
