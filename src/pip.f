      program pipfit
c
      implicit double precision (a-h,p-z)
c
c MPI
      include 'mpif.h'
      integer my_id,nprocs,ierr
      integer status(MPI_STATUS_SIZE)

      include 'param.inc'
      dimension coef(maxcoef),sig(maxdata)
      dimension basis(maxterm)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension cut(100),nc(100),na(100),wc(100),wa(100),
     &          erra(100),errc(100)
      dimension x(maxatom),y(maxatom),z(maxatom)

      character*16 datafit,datatest
      character*2 dum

      integer batches,natom1,a,b
      logical linteronly

      common/foox/rrr,nncoef,natom1,linteronly

ccccc MPI
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

      IF (my_id.eq.0) THEN
        timestart=MPI_WTIME()
        write(6,*)"Began at ",timestart
      ENDIF
!      write(6,"(a, i3)") " MPI current proc: ", my_id

      IF (my_id.eq.0) THEN
!      write(6,"(a, i3)") " MPI num procs: ", nprocs

      write(6,'(100("*"))')
      write(6,*)
      write(6,*)"PIP: A Fortran code for fitting permutationally",
     &  " invariant polynomials"
      write(6,*)"     Daniel R. Moberg and Ahren W. Jasper, Argonne,",
     & " 2020"
      write(6,*)

      read(5,*)datafit,datatest
      read(5,*)epsilon,vvref
      read(5,*)ncut,(cut(j),j=1,ncut)
      ENDIF

      call MPI_BCAST(datafit, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(datatest,1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(epsilon, 1, MPI_DOUBLE_PRECISION, 0,
     &               MPI_COMM_WORLD, ierr)
      call MPI_BCAST(vvref, 1, MPI_DOUBLE_PRECISION, 0,
     &               MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ncut, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cut, ncut, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call prepot ! more parameters are read in prepot

!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!     IF (.false.) THEN ! Not fitting function for now
!     ncoef=nncoef
!
!     open(7,file=datafit)
!     read(7,*)ndat2
!     if (ndat2.gt.maxdata) then
!       write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata
!       stop
!     endif
!
!     ndat=0
!     vvmin=1.d50
!     do i=1,ndat2
!       read(7,*)natom
!       read(7,*)iz,dum,vvx
!       do j=1,natom
!         read(7,*)dum,x(j),y(j),z(j)
!       enddo
!
!       if (vvx.lt.cut(1)) then
!         if (vvx.lt.vvmin) vvmin=vvx
!         ndat=ndat+1
!         vv(ndat)=vvx
!         ii=0
!         a=natom
!         if (linteronly) a=natom1
!         do j=1,a
!           b=j
!           if (linteronly) b=natom1
!           do k=b+1,natom
!             ii=ii+1
!             rrr(ndat,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2
!    &                                         +(z(j)-z(k))**2)
!           enddo  
!         enddo  
!       sig(ndat)=1.d0/(epsilon/(dabs(vv(ndat)-vvref)+epsilon))
!       endif
!     enddo
!     close(7)
!
!     IF (my_id.eq.0) THEN
!     write(6,'(34("*  "))')
!     write(6,*)
!     write(6,*)"Fitting to the training data in ",datafit
!     write(6,*)"Using ",ndat," of the ",ndat2," provided data"
!     write(6,*)"Out of sample testing using data in ",datatest
!     write(6,*)"Weight function parameters: epsilon = ",
!    & epsilon," and vref = ",vvref
!     write(6,*)
!     ENDIF
!
!     call svdfit(vv,sig,ndat,coef,ncoef,ndat,ncoef)
!
!     do i=1,ncut
!     erra(i)=0.d0 ! RMSE for E < cut(i)
!     errc(i)=0.d0 ! RMSE for cut(i) < E < cut(i-1)
!     wc(i)=0.d0
!     wa(i)=0.d0
!     nc(i)=0
!     na(i)=0
!     enddo
!     vvimin=1d50
!     vvxmin=1d50
!     open(56,file="coef.dat")
!     do i=1,ndat
!       call funcs1(i,basis,ncoef) 
!      write(88,'(i10,10000e18.6)')i,1./sig(i),vv(i),(basis(j),j=1,ncoef)
!       vvx=0.d0
!       do j=1,ncoef
!         vvx=vvx+coef(j)*basis(j)
!         if (i.eq.1) write(56,'(i5,e20.10)')j,coef(j)
!       enddo
!       do k=1,ncut
!       if (vv(i).lt.cut(k)) then
!         erra(k)=erra(k)+(vvx-vv(i))**2/sig(i)**2
!         wa(k)=wa(k)+1.d0/sig(i)**2
!         na(k)=na(k)+1
!       endif
!       enddo
!       do k=1,ncut-1
!       if (vv(i).lt.cut(k).and.vv(i).ge.cut(k+1)) then
!         errc(k)=errc(k)+(vvx-vv(i))**2/sig(i)**2
!         wc(k)=wc(k)+1.d0/sig(i)**2
!         nc(k)=nc(k)+1
!       endif
!       enddo
!       if (vvx.lt.vvxmin) then
!         vvxmin=vvx
!         vvxa=vv(i)
!         vvxb=vvx
!       endif
!       if (vv(i).lt.vvimin) then
!         vvimin=vv(i)
!         vvia=vv(i)
!         vvib=vvx
!       endif
!     enddo
!     close(56)
!     errc(ncut)=erra(ncut)
!     wc(ncut)=wa(ncut)
!     nc(ncut)=na(ncut)
!     print *,'Fitting errors between (below) user-provided energies'
!     do k=1,ncut
!     erra(k)=dsqrt(erra(k)/wa(k))
!     errc(k)=dsqrt(errc(k)/wc(k))
!     write(6,'(f15.5,i10,f15.5," ( ",i10,f15.5," ) ")')
!    &  cut(k),nc(k),errc(k),na(k),erra(k)
!     enddo
!     print *
!     print *,"Comparision of low energy points found while fitting"
!     print *,"                    data           fit        ",
!    & "    difference "
!     write(6,'(a,3f15.5)')"   data minimum ",vvia,vvib,vvia-vvib
!     write(6,'(a,3f15.5)')" fitted minumum ",vvxa,vvxb,vvxa-vvxb
!
!     test set
!     open(7,file=datatest)
!     rewind(7)
!     read(7,*)ndat2
!     if (ndat2.gt.maxdata) then
!         write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata
!         stop
!     endif
! 
!       ndat=0
!       do i=1,ndat2
!         read(7,*)
!         read(7,*)iz,dum,vvx
!         do j=1,natom
!           read(7,*)dum,x(j),y(j),z(j)
!         enddo
!
!         xcut=cut(1) 
!         if (ncut.gt.1) xcut=cut(2) 
!         if (vvx.lt.xcut) then
!           ndat=ndat+1
!           vv(ndat)=vvx
!           ii=0
!           a=natom
!           if (linteronly) a=natom1
!           do j=1,a
!             b=j
!             if (linteronly) b=natom1
!             do k=b+1,natom
!               ii=ii+1
!       rrr(ndat,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
!             enddo  
!           enddo  
!        sig(ndat)=1.d0/(epsilon/(dabs(vv(ndat)-vvref)+epsilon))
!         endif
!       enddo  
!
!       err3=0.d0
!       ndat3=0
!       wn=0.d0
!       do i=1,ndat
!         call funcs1(i,basis,ncoef)
!         vvx=0.d0
!         do l=1,ncoef
!           vvx=vvx+coef(l)*basis(l)
!         enddo
!         err3=err3+(vvx-vv(i))**2/sig(i)**2
!         wn=wn+1.d0/sig(i)**2
!       enddo
!       err3=dsqrt(err3/wn)
!
!       write(6,*)
!       write(6,*)"Out of sample test set error: ",err3 
!       write(6,*)
!       write(6,*)"Optimized coefficients written to coef.dat"
!       write(6,*)
!       write(6,'(100("*"))')
!
!       ix=2
!       if (ncut.eq.1) ix=1 
!       print *
!       print *,"summary",ncoef,erra(ix),err3
!       print *
!
!     ENDIF ! false

      IF (my_id.eq.0) THEN
        timeend=MPI_WTIME()
!        write(6,*)"Ended at:"timeend
        write(6,*)"Total time:",(timeend-timestart)
      ENDIF

      call MPI_FINALIZE(ierr)

      end
