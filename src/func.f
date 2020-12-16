********************************************

      subroutine prepot

      implicit double precision(a-h,p-z)
c MPI
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      include 'param.inc'
      integer dgroup,a,b,ugg(maxatom,maxatom),g1,g2,ghi,glo,bad(maxterm)
      dimension iagroup(maxatom),itype(maxterm),
     &  ind(maxterm,maxpair),iatom(maxperm,maxatom),indord(maxterm),
     &  idum(maxatom),nngroup(maxatom),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),r(maxpair),
     &  dbasisdr(maxterm,maxpair),rrr(maxdata,maxpair),
     &  index(maxatom,maxatom),ix(maxperm,maxpair),
     &  dgroup(maxchan,maxatom)
      integer, allocatable :: ifailmpi(:),ibmpi(:),ifail(:),ibasis(:)
      integer, allocatable :: ib(:),term1(:),term2(:),ttlist(:),disps(:)
      character*2 symb(maxatom),dum
      logical lreadbasis,lnodisc,lnointra,linteronly,ldebug

      common/foox/rrr,nncoef,natom1,linteronly

      save npairs,nterm,ind,ibasis

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

c read
      ldebug=.true.
      ldebug=.false.

      IF (my_id.eq.0) THEN
      read(5,*)natom
      read(5,*)(symb(k),k=1,natom)
      read(5,*)(iagroup(k),k=1,natom)
      read(5,*)lreadbasis,ipow,ipowt,imode
      write(6,'(100("*"))')
      write(6,*)
      if (lreadbasis) then
      write(6,*)"Reading the PIP",ipow,ipowt," basis"
      else
      write(6,*)"Generating a PIP",ipow,ipowt,
     & " basis and writing it to basis.dat"
      endif
      if (imode.eq.0) then     ! use all terms
        lnodisc = .false.      ! remove disconnected terms?
        lnointra = .false.     ! remove intramolecular-only terms?
        linteronly = .false.   ! use only intermolecular terms?
        write(6,*)"Mode = 0: Using all bond distances"
      elseif (imode.eq.1) then ! remove disconnected and intra-only terms
        lnodisc = .true.
        lnointra = .true.
        linteronly = .false.
        write(6,*)"Mode = 1: Using all bond distances and removing",
     & " disconnected and intramolecular-only terms"
      elseif (imode.eq.2) then ! remove disconnected terms
        lnodisc = .true.
        lnointra = .false.
        linteronly = .false.
        write(6,*)"Mode = 2: Using all bond distances and removing",
     & " disconnected terms"
      elseif (imode.eq.3) then ! remove intramolecular terms
        lnodisc = .false.
        lnointra = .true.
        linteronly = .false.
        write(6,*)"Mode = 3: Using all bond distances and removing",
     & " intramolecular-only terms"
      elseif (imode.eq.-1) then ! use intermolecular distances only
        lnodisc = .false.
        lnointra = .false.
        linteronly = .true.
        write(6,*)"Mode = -1: Using only intermolecular bond distances"
      endif
      write(6,*)
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(natom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(symb,natom, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(iagroup,natom,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lreadbasis,1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ipow, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ipowt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(imode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_BCAST(lnodisc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lnointra, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(linteronly,1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (lnodisc.or.lnointra.or.linteronly) then
      IF (my_id.eq.0) THEN
        read(5,*)nchan
        do i=1,nchan
          read(5,*)(dgroup(i,k),k=1,natom)
        enddo
      ENDIF
      call MPI_BCAST(nchan, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dgroup, maxchan*maxatom,
     &               MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (linteronly) then
        if(nchan.ne.1) then
          IF (my_id.eq.0) write(6,*)
     &     "Mode = -1 requires 1 fragment channel with 2 groups"
           stop
        endif
        natom1=0        ! bimolecular only
        natom2=0
        do i=1,nchan
        do j=1,natom
          if (dgroup(i,j).eq.1) natom1=natom1+1
          if (dgroup(i,j).eq.2) natom2=natom2+1
          if (dgroup(i,j).ne.1.and.dgroup(i,j).ne.2) then
            IF (my_id.eq.0) write(6,*)
     &       "Mode = -1 requires 1 fragment channel with 2 groups"
            stop
          endif
        enddo
        enddo
        if (natom.ne.natom1+natom2) then
          IF (my_id.eq.0) print *,"Error specifying molecular groups"
          stop
        endif
      endif

      npairs=natom*(natom-1)/2                 ! Number of interatomic pairs
      if (linteronly) npairs = natom1*natom2

      IF (my_id.eq.0) THEN
      write(6,'(34("*  "))')
      write(6,*)
      write(6,*)"Atoms"
      write(6,'(1x,a8,100a5)')"Symbol",(symb(i),i=1,natom)
      write(6,'(a8,100i5)')"Index",(i,i=1,natom)
      write(6,'(a8,100i5)')"Group",(iagroup(i),i=1,natom)
      if (lnodisc.or.lnointra.or.linteronly) then
        write(6,'(3x,a)')"Fragment Channels"
        do i=1,nchan
          write(6,'(2x,i5,a,100i5)')i,":",(dgroup(i,j),j=1,natom)
        enddo
      endif
      ENDIF

      if (natom.gt.maxatom) then
        IF (my_id.eq.0) 
     &      print *,"natom (",natom,") > maxatom (",maxatom,")"
        stop
      endif

      if (npairs.gt.maxpair) then
        IF (my_id.eq.0) 
     &    print *,"npairs (",npairs,") > maxpair (",maxpair,")"
        stop
      endif

      do i=1,nchan
      do j=1,natom
        if (dgroup(i,j).ne.1.and.dgroup(i,j).ne.2) then
          IF (my_id.eq.0) write(6,*)
     &      "This version of the code is limited to ",
     &      "bimolecular fragments"
          stop
        endif
      enddo
      enddo

ccc GENERATE BASIS ccc
      IF (.not.lreadbasis) THEN
ccc GENERATE BASIS ccc

c generate atom permutation lists
      do i=1,natom
        nngroup(i)=0
      enddo
      ngroup=1
      do i=1,natom
        if (iagroup(i).gt.ngroup) ngroup=iagroup(i)
        nngroup(iagroup(i))=nngroup(iagroup(i))+1
      enddo

      nn=0

      do i=1,ngroup
        n=0
        do k=1,natom
          if (iagroup(k).eq.i) then
            n=n+1
            idum(n)=k
          endif
        enddo
      
        npermute=0
        call heapp(idum,n,n,idum2,npermute)
        nperm0(i)=nn+1
        nperm1(i)=nn+npermute

        do k=1,npermute
          nn=nn+1
          m=0
          do j=1,natom
            idum3(nn,j)=0
            if (iagroup(j).eq.i) then
              m=m+1
              idum3(nn,j)=idum2(k,m)
            endif
          enddo ! j=1,natom
        enddo ! k=1,npermute
      enddo ! i=1,ngroup

      ntmp=1
      do i=1,ngroup
        idum(i)=nperm0(i)
        ntmp=ntmp*(nperm1(i)-nperm0(i)+1)
      enddo

      npermute=0
      do while (.true.)
        npermute=npermute+1
        if (npermute.gt.maxperm) then
          IF (my_id.eq.0) 
     &      print *,"npermute (",npermute,") > maxperm (",maxperm,")"
          stop
        endif

        do i=1,natom
          iatom(npermute,i)=0
          do j=1,ngroup
            iatom(npermute,i)=iatom(npermute,i)+idum3(idum(j),i)
          enddo
        enddo

        idum(ngroup)=idum(ngroup)+1
 777    continue

        do i=1,ngroup
          if (idum(i).gt.nperm1(i)) then
            if (i.eq.1) go to 778
            idum(i)=nperm0(i)
            idum(i-1)=idum(i-1)+1
            go to 777
          endif
        enddo 

      enddo
 778  continue

      IF (my_id.eq.0) THEN
      print *
      print *,'Atom permutations = ',npermute
      do i=1,min(npermute,100) ! prints only first 100 permutations
        write(6,'(i7,a,1000i5)')i,":",(iatom(i,j),j=1,natom)
      enddo
      if (npermute.gt.100) write(6,*)" ** Truncating list **"
      write(6,*)
      ENDIF

      ii=0
      a=natom
      if (linteronly) a=natom1
      do i=1,a
        b=i
        if (linteronly) b=natom1
        do j=b+1,natom
          ii=ii+1
          index(i,j)=ii ! pair number
          index(j,i)=ii
        enddo
      enddo

      IF (my_id.eq.0) THEN
      if (ldebug) then ! DEBUG OUTPUT
      write(6,*)"Pair permutation list"
      if (.not.linteronly) then
      write(6,'(14x,1000(a3,"- ",a3,2x))')
     &   ((symb(i),symb(j),j=1+i,natom),i=1,natom) 
      write(6,'(13x,1000(i3," -",i3,2x))')((i,j,j=1+i,natom),
     &   i=1,natom)
      else
      write(6,'(14x,1000(a3,"- ",a3,2x))')
     &   ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1)
      write(6,'(13x,1000(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
     &   i=1,natom1)
      endif
      endif ! DEBUG OUTPUT
      ENDIF

      do ii=1,npermute
        iix=0
        a=natom
        if (linteronly) a=natom1
        do i=1,a
          b=i
          if (linteronly) b=natom1
          do j=b+1,natom
            iix=iix+1
            ix(ii,iix)=index(iatom(ii,i),iatom(ii,j))
          enddo
        enddo
        if (ii.le.100.and.ldebug) write(6,'(i7,a,1000i10)')
     &     ii,":",(ix(ii,iix),iix=1,npairs)
      enddo
      IF (my_id.eq.0) THEN
      if (npermute.gt.100.and.ldebug) write(6,*)" ** Truncating list **"
      if (ldebug) write(6,*)
      ENDIF

c generate terms using user-supplied power constraints
      ii=1
      indtot=0
      do i=1,npairs
        ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
          IF (my_id.eq.0) 
     &      print *,"number of terms (",ii,") > maxterm (",maxterm,")"
          stop
        endif

        indord(ii-1)=indtot
        do i=1,npairs
          ind(ii,i)=ind(ii-1,i)
        enddo
        ind(ii,npairs)=ind(ii,npairs)+1
 300    continue

        indtot=0
        do i=1,npairs
          indtot=indtot+ind(ii,i)
          if (ind(ii,i).gt.ipow.or.indtot.gt.ipowt) then ! ipow(i) would allow atom-atom-type-dependent limits
            if (i.eq.1) go to 400
            ind(ii,i)=0
            ind(ii,i-1)=ind(ii,i-1)+1
            go to 300
          endif
        enddo
      enddo
 400  continue
      nterm=ii-1

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (my_id.eq.0) t1=MPI_WTIME()
ccccc Symmetrize
      nbasis=0
      IF (my_id.eq.0) write(6,*)"Symmetrizing the expansion"

      allocate(term1(nproc))
      allocate(term2(nproc))
      allocate(ttlist(nproc))
      allocate(disps(nproc))

c calculate number of terms assigned to each process
      IF (my_id.eq.0) THEN
      do i=0,nproc-1
        term1(i+1)=(i*nterm)/nproc+1
        term2(i+1)=((i+1)*nterm)/nproc
        ttlist(i+1)=term2(i+1)-term1(i+1)+1
      enddo
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(term1, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(term2, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ttlist,nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

c displacements
      disps(1)=0
      do i=2,nproc
        disps(i)=disps(i-1)+ttlist(i-1)
      enddo

      it=my_id+1

      allocate(ifailmpi(ttlist(it)))
      allocate(ibmpi(ttlist(it)))

ccccc MPI DO
      if (my_id.eq.0) tl1=MPI_WTIME()
      DO ii=term1(it),term2(it)
!        IF (my_id.eq.0) THEN
          if (mod(ii,100).eq.0) print *,"  step ",ii," of ",nterm,
     &                " on proc ",my_id
!        ENDIF
        ifailmpi(ii-term1(it)+1)=0
        ! MOVE MPI LOOP TO HERE?
        do i=ii-1,1,-1
          if (indord(ii).eq.indord(i)) then
          do j=1,npermute
            ifailmpi(ii-term1(it)+1)=1
            do k=1,npairs
              if (ind(i,k).ne.ind(ii,ix(j,k))) then
                ifailmpi(ii-term1(it)+1)=0
                exit
              endif
            enddo
            if (ifailmpi(ii-term1(it)+1).eq.1) go to 1010
          enddo
          endif
        enddo
 1010 continue
        if (ifailmpi(ii-term1(it)+1).ne.0) then
          ibmpi(ii-term1(it)+1)=i
        endif
      ENDDO
ccccc MPI DO

      allocate(ifail(nterm))
      allocate(ibasis(nterm))
      allocate(ib(nterm))

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_id.eq.0) tl2=MPI_WTIME()
      call MPI_ALLGATHERV(ifailmpi, ttlist(it), MPI_INTEGER, ifail, 
     &                    ttlist, disps, MPI_INTEGER,
     &                    MPI_COMM_WORLD, ierr)
      call MPI_ALLGATHERV(ibmpi, ttlist(it), MPI_INTEGER, 
     &                    ib, ttlist, disps, MPI_INTEGER, 
     &                    MPI_COMM_WORLD, ierr)

      DO ii=1,nterm
        if (ifail(ii).eq.0) then
          nbasis=nbasis+1
          ibasis(ii)=nbasis
        else
          ibasis(ii)=ibasis(ib(ii))
        endif
      ENDDO

      IF (my_id.eq.0) write(6,*)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_id.eq.0) then
        t2=MPI_WTIME()
        write(6,'(" CPU time in MPI loop is ",f10.5," s",/)')(tl2-tl1)
        write(6,'(" CPU time in Symmetrizing is ",f10.5," s",/)')(t2-t1)
      endif

      nncoef=nbasis

c remove unconnected and intramolecular only terms if required

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      IF (my_id.eq.0) THEN
      IF (lnodisc.or.lnointra) THEN

        do k=1,nterm
          itype(k)=0
        enddo

        DO m=1,nchan     ! loop over fragment channels, m
        do k=1,nterm     ! loop over terms, k

          do i=1,natom
          do j=1,natom
            ugg(i,j)=0
          enddo
          enddo
          ghi=0
          glo=natom
          ii=0           ! loop over atom pairs i<j
          do i=1,natom
          do j=i+1,natom
            ii=ii+1
            g1=min(dgroup(m,i),dgroup(m,j))
            g2=max(dgroup(m,i),dgroup(m,j))
            ghi=max(g1,g2,ghi)
            glo=min(g1,g2,glo)
            ugg(g1,g2)=ugg(g1,g2)+ind(k,ii)
          enddo
          enddo

          ng=0
          do i=glo,ghi   ! loop over groups
            if (ugg(i,i).ne.0) ng=ng+1
          enddo

          ityp=0
          if (ng.eq.0.or.ugg(1,2).ne.0) then
            ityp=0 ! ok
          elseif (ng.eq.1.and.lnointra) then
            ityp=1 ! intra only
          elseif (ng.eq.2.and.lnodisc) then
            ityp=2 ! unconnected
          endif
          itype(ibasis(k)) = max(itype(ibasis(k)),ityp) ! If terms have different statuses in different channels, retain the larger assignment
        enddo

        nin=0
        nun=0
        do i=1,nbasis
          if (itype(i).eq.1) nin=nin+1
          if (itype(i).eq.2) nun=nun+1
          if (itype(i).eq.1.or.itype(i).eq.2) bad(nin+nun)=i
        enddo
        IF (my_id.eq.0) THEN
        if (m.eq.1) then
          print *,"Analyzing fragment channel ",m
        else
          print *,"Updating group types for fragment channel ",m
        endif
        if (lnointra) 
     &      print *,"Found ",nin," intramolecular-only groups"
        if (lnodisc) 
     &      print *,"Found ",nun," unconnected groups"
        print *
        ENDIF

        ENDDO

        nt=nterm
        nb=nbasis

        do ii=nterm,1,-1
          if (itype(ibasis(ii)).ne.0) then ! bad term
            nterm=nterm-1
            do jj=ii,nterm
              ibasis(jj)=ibasis(jj+1)
              do k=1,npairs
                ind(jj,k)=ind(jj+1,k)
              enddo
            enddo
          endif
        enddo ! ii=nterm,1,-1

        do ii=1,nterm
          nx=ibasis(ii)
          do i=1,nin+nun
            if (nx.gt.bad(i)) ibasis(ii)=ibasis(ii)-1
          enddo
        enddo ! ii=1,nterm

        nbasis=nbasis-nin-nun

        IF (my_id.eq.0) THEN
        print *,"Removing these groups results in the",
     & " following reductions:"
        print *,"Terms:  ",nt," --> ",nterm," ("
     &             ,(dble(nterm)/dble(nt)*100.)," % )"
        print *,"Groups: ",nb," --> ",nbasis," ("
     &             ,(dble(nbasis)/dble(nb)*100.)," % )"
        print *
        ENDIF
      ENDIF
      ENDIF ! my_id=0
      nncoef=nbasis

      IF (my_id.eq.0) THEN
      write(6,*)"Finished preparing the basis and writing it to",
     & " basis.dat"
      write(6,*)"This basis has ",nterm," (",nncoef,") terms (groups)"
      write(6,*)
      open(55,file="basis.dat")
      a=natom
      if (linteronly) a=natom1
      write(55,*)a,npairs,nncoef,nterm,
     & " ! atoms, atom pairs, groups/coefficients, terms"
      write(55,*)"   TERM   GROUP : EXPONENTS"
      do ii=1,nterm
        write(55,'(2i8," : ",1000i8)')
     &        ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)
      ENDIF

ccc READ BASIS ccc
      ELSE
        IF (my_id.eq.0) THEN
        open(55,file="basis.dat")
        a=natom
        if (linteronly) a=natom1
        read(55,*)a,npairs,nncoef,nterm
        read(55,*)
        do i=1,nterm
          read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
        enddo
        close(55)
        write(6,*)
        write(6,*)"Reading the basis from basis.dat"
        write(6,*)"This basis has ",nterm," (",nncoef,") terms (groups)"
        write(6,*)
        ENDIF
      ENDIF

      return

      entry funcs1(iii,basis,ncoef)

      do j=1,ncoef
        basis(j)=0.d0
!        do j=1,npairs
!          dbasisdr(i,j)=0.d0
!        enddo
      enddo

      do j=1,npairs
        r(j)=dexp(-rrr(iii,j))
!        r(j)=dexp(-rrr(iii,j)*autoang)
      enddo

      do i=1,nterm
        arg=1.d0
        do j=1,npairs
          arg=arg*(r(j)**ind(i,j))
        enddo
        basis(ibasis(i))=basis(ibasis(i))+arg
!        do j=1,npairs
!          dbasisdr(ibasis(i),j)=dbasisdr(ibasis(i),j)   ! dV/dy * dy/dr
!     &                         -arg*dble(ind(i,j))*autoang
!        enddo
      enddo

      return 
      end

***************************************************

      recursive subroutine heapp(ia,size,n,iia,ii)

      include 'param.inc'
      integer i,n,size,ii
      integer ia(maxatom)
      integer iia(maxperm,maxatom)
      integer iagroup(maxatom)

      if (size.eq.1) then
         ii=ii+1
         do i=1,n
           iia(ii,i)=ia(i)
         enddo
        return
      endif

      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
      endif

      enddo

      end subroutine

***************************************************

