********************************************

      subroutine prepot

      implicit double precision(a-h,p-z)
c MPI
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      include 'param.inc'
      integer changroup,a,b,ugg(maxatom,maxatom),
     &  g1,g2,ghi,glo,bad(maxterm),nchangroup(maxatom)
      integer, allocatable :: ifailmpi(:),ibmpi(:),ifail(:)
      integer, allocatable :: ib(:),term1(:),term2(:),ttlist(:),disps(:)
      dimension iagroup(maxatom),itype(maxterm),
     &  ind(maxterm,maxpair),iatom(maxperm,maxatom),indord(maxterm),
     &  idum(maxatom),nngroup(maxatom),idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),r(maxpair),ibasis(maxterm),
     &  dbasisdr(maxterm,maxpair),rrr(maxdata,maxpair),
     &  index(maxatom,maxatom),ix(maxperm,maxpair),
     &  changroup(maxchan,maxatom)
      character*2 symb(maxatom),dum
      logical lreadbasis,lnounc,lnointra,linter,ldebug,lwrite(-1:100)

      common/foox/rrr,nncoef,natom1,linter,lwrite

      save npairs,nterm,ind,ibasis

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

c read
      ldebug=.false.
      if (lwrite(1)) ldebug=.true.             ! print some extra information

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
        write(6,'(a,2i5,a)')" Generating a PIP",ipow,ipowt,
     &   " basis and writing it to basis.dat"
      endif

c     logical flags
c       lnounc = remove unconnected terms?
c       lnointra = remove intramolecular terms?
c       linter = use intermolecular distances only?
      if (imode.eq.0) then     ! use all terms
        lnounc = .false.
        lnointra = .false.
        linter = .false.
        write(6,*)"Mode = 0: Using all bond distances"
      elseif (imode.eq.1) then ! remove unconnected terms
        lnounc = .true.
        lnointra = .false.
        linter = .false.
        write(6,*)"Mode = 1: Using all bond distances and removing",
     & " unconnected terms"
      elseif (imode.eq.2) then ! remove unconnected and intra-only terms
        lnounc = .true.
        lnointra = .true.
        linter = .false.
        write(6,*)"Mode = 2: Using all bond distances and removing",
     & " unconnected and intramolecular-only terms"
      elseif (imode.eq.3) then ! remove intramolecular terms
        lnounc = .false.
        lnointra = .true.
        linter = .false.
        write(6,*)"Mode = 3: Using all bond distances and removing",
     & " intramolecular-only terms"
      elseif (imode.eq.-1) then ! use intermolecular distances only
        lnounc = .false.
        lnointra = .false.
        linter = .true.
        write(6,*)"Mode = -1: Using only intermolecular bond distances"
      endif
      write(6,*)
      ENDIF

      if (my_id.eq.0.and.nproc.gt.1) 
     & write(6,'(1x,a,i10,a,/)')"Distributing job over ",nproc," cores"

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(natom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(symb,natom, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(iagroup,natom,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lreadbasis,1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ipow, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ipowt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_BCAST(lnounc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lnointra, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(linter,1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (lnounc.or.lnointra.or.linter) then
      IF (my_id.eq.0) THEN
        read(5,*)nchan
        if (nchan.gt.maxchan) then
          write(6,*)" nchan (",nchan,") > maxchan (",maxchan,")"
          stop
        endif
        do i=1,nchan
          read(5,*)(changroup(i,k),k=1,natom)
        enddo
      ENDIF
      call MPI_BCAST(nchan, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(changroup, maxchan*maxatom,
     &               MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (linter) then
        if(nchan.ne.1) then
          IF (my_id.eq.0) write(6,*)
     &     "Mode = -1 requires 1 fragment channel with 2 groups"
          stop
        endif
        natom1=0        ! bimolecular only
        natom2=0
        do i=1,nchan
        do j=1,natom
          if (changroup(i,j).eq.1) natom1=natom1+1
          if (changroup(i,j).eq.2) natom2=natom2+1
          if (changroup(i,j).ne.1.and.changroup(i,j).ne.2) then
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

      if (natom.gt.maxatom) then
        IF (my_id.eq.0)
     &    print *,"natom (",natom,") > maxatom (",maxatom,")"
        stop
      endif

      if (npairs.gt.maxpair) then
        IF (my_id.eq.0)
     &    print *,"npairs (",npairs,") > maxpair (",maxpair,")"
        stop
      endif

      npairs=natom*(natom-1)/2                 ! Number of interatomic pairs
      if (linter) npairs = natom1*natom2

      termest=1.d0 ! number of terms for ipow=ipowt
      do i=max(npairs,ipowt)+1,npairs+ipowt
        termest=termest*dble(i)
      enddo
      do i=2,min(npairs,ipowt)
        termest=termest/dble(i)
      enddo

      if (termest.gt.dble(maxterm)) then
        IF (my_id.eq.0)
     &    print *,"Estimated number of terms (",
     &    termest,") > maxterm (",maxterm,")"
        stop
      endif

      IF (my_id.eq.0) THEN
      write(6,'(34("*  "))')
      write(6,*)
      write(6,*)"Atoms"
!      write(6,'(1x,a8,100a5)')"Symbol",(symb(i),i=1,natom)
!      write(6,'(a8,100i5)')"Index",(i,i=1,natom)
!      write(6,'(a8,100i5)')"Group",(iagroup(i),i=1,natom)
      write(6,'(1x,a8,*(a5))')"Symbol",(symb(i),i=1,natom)
      write(6,'(a8,*(i5))')"Index",(i,i=1,natom)
      write(6,'(a8,*(i5))')"Group",(iagroup(i),i=1,natom)
      ENDIF

      if (lnointra.or.lnounc.or.linter) then
      do i=1,nchan
      do j=1,natom
       if (changroup(i,j).gt.nchangroup(i)) nchangroup(i)=changroup(i,j)
      enddo
      enddo
      endif

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

      if (my_id.eq.0) 
     &    write(6,'(3x,a,i10,/)')"Number of permutations = ",ntmp
      if ((lnounc.or.lnointra.or.linter).and.my_id.eq.0) then
        write(6,'(3x,a)')"Fragment Channels"
        do i=1,nchan
!          write(6,'(1x,i5,a,100i5)')i," :",(changroup(i,j),j=1,natom)
          write(6,'(1x,i5,a,*(i5))')i," :",(changroup(i,j),j=1,natom)
        enddo
      endif

      if (ntmp.gt.maxperm) then
        IF (my_id.eq.0) 
     &    print *,"npermute (",ntmp,") > maxperm (",maxperm,")"
        stop
      endif

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
      if (ldebug) then ! DEBUG OUTPUT
      write(6,*)
      write(6,*)'Atom permutations = ',npermute
      do i=1,min(npermute,50) ! print up to 50 permutations
!        write(6,'(i6,a,1000i5)')i," :",(iatom(i,j),j=1,natom)
        write(6,'(i6,a,*(i5))')i," :",(iatom(i,j),j=1,natom)
      enddo
      if (npermute.gt.50) write(6,*)" ** Truncating list **"
      endif
      if (lwrite(10)) then
      open(10,file="basisinfo.dat")
      write(10,*)
      write(10,*)'Atom permutations = ',npermute
      do i=1,npermute
!        write(10,'(i7,a,1000i5)')i," :",(iatom(i,j),j=1,natom)
        write(10,'(i7,a,*(i5))')i," :",(iatom(i,j),j=1,natom)
      enddo
      endif
      ENDIF

      ii=0
      a=natom
      if (linter) a=natom1
      do i=1,a
        b=i
        if (linter) b=natom1
        do j=b+1,natom
          ii=ii+1
          index(i,j)=ii ! pair number
          index(j,i)=ii
        enddo
      enddo

      IF (my_id.eq.0) THEN
      if (ldebug) then ! DEBUG OUTPUT
        write(6,*)
        write(6,*)"Pair permutations"
        if (.not.linter) then
!          write(6,'(14x,1000(a3,"- ",a3,2x))')
          write(6,'(14x,*(a3,"- ",a3,2x))')
     &       ((symb(i),symb(j),j=1+i,natom),i=1,natom) 
!          write(6,'(13x,1000(i3," -",i3,2x))')((i,j,j=1+i,natom),
          write(6,'(13x,*(i3," -",i3,2x))')((i,j,j=1+i,natom),
     &       i=1,natom)
        else
!          write(6,'(14x,1000(a3,"- ",a3,2x))')
          write(6,'(14x,*(a3,"- ",a3,2x))')
     &       ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1)
!          write(6,'(13x,1000(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
          write(6,'(13x,*(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
     &       i=1,natom1)
        endif
      endif ! DEBUG OUTPUT
      if (lwrite(10)) then
        write(10,*)
        write(10,*)"Pair permutations"
        if (.not.linter) then
!          write(10,'(14x,1000(a3,"- ",a3,2x))')
          write(10,'(14x,*(a3,"- ",a3,2x))')
     &       ((symb(i),symb(j),j=1+i,natom),i=1,natom)
!          write(10,'(13x,1000(i3," -",i3,2x))')((i,j,j=1+i,natom),
          write(10,'(13x,*(i3," -",i3,2x))')((i,j,j=1+i,natom),
     &       i=1,natom)
        else
!          write(10,'(14x,1000(a3,"- ",a3,2x))')
          write(10,'(14x,*(a3,"- ",a3,2x))')
     &       ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1)
!          write(10,'(13x,1000(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
          write(10,'(13x,*(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
     &       i=1,natom1)
        endif
      endif
      ENDIF

      do ii=1,npermute
        iix=0
        a=natom
        if (linter) a=natom1
        do i=1,a
          b=i
          if (linter) b=natom1
          do j=b+1,natom
            iix=iix+1
            ix(ii,iix)=index(iatom(ii,i),iatom(ii,j))
          enddo
        enddo
        IF (my_id.eq.0) THEN
!        if (ii.le.50.and.ldebug) write(6,'(i6,a,1000i10)')
        if (ii.le.50.and.ldebug) write(6,'(i6,a,*(i10))')
     &     ii," :",(ix(ii,iix),iix=1,npairs)
!        if (lwrite(10)) write(10,'(i7,a,1000i10)')
        if (lwrite(10)) write(10,'(i7,a,*(i10))')
     &     ii," :",(ix(ii,iix),iix=1,npairs)
        ENDIF
      enddo
      IF (my_id.eq.0) THEN
      if (npermute.gt.50.and.ldebug) write(6,*)" ** Truncating list **"
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

c     Symmetrize the basis
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_id.eq.0) t1=MPI_WTIME()

      nbasis=0
      IF (my_id.eq.0) write(6,*)
      IF (my_id.eq.0) write(6,*)"Symmetrizing the expansion"

!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      it=my_id+1

      ndata=nterm/nproc+1*min(mod(nterm,nproc),1)

      allocate(ifailmpi(nterm))
      allocate(ibmpi(nterm))

      do i=1,nterm
        ifailmpi(i)=0
        ibmpi(i)=0
      enddo

      if (my_id.eq.0) tl1=MPI_WTIME()
ccccc MPI DO LOOP
      DO ii=my_id+1,nterm,nproc
        if (mod(ii,1000).eq.0) write(6,*)"  step ",ii," of ",nterm
c     &                ," on proc ",my_id
        ifailmpi(ii)=0
        do i=ii-1,1,-1
          if (indord(ii).eq.indord(i)) then
          do j=1,npermute
            ifailmpi(ii)=1
            do k=1,npairs
              if (ind(i,k).ne.ind(ii,ix(j,k))) then
                ifailmpi(ii)=0
                exit
              endif
            enddo
            if (ifailmpi(ii).eq.1) go to 1010
          enddo
          endif
        enddo
 1010 continue
        if (ifailmpi(ii).ne.0) then
          ibmpi(ii)=i
        endif
      ENDDO
ccccc MPI DO LOOP

      allocate(ifail(nterm))
      allocate(ib(nterm))

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_id.eq.0) tl2=MPI_WTIME()
c      call MPI_ALLREDUCE(ifailmpi, ifail, nterm, MPI_INTEGER,
c     &                   MPI_SUM, MPI_COMM_WORLD, ierr)
c      call MPI_ALLREDUCE(ibmpi, ib, nterm, MPI_INTEGER, 
c     &                   MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(ifailmpi, ifail, nterm, MPI_INTEGER,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(ibmpi, ib, nterm, MPI_INTEGER, 
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      DO ii=1,nterm
        if (ifail(ii).eq.0) then
          nbasis=nbasis+1
          ibasis(ii)=nbasis
        else
          ibasis(ii)=ibasis(ib(ii))
        endif
      ENDDO

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_id.eq.0) then
       t2=MPI_WTIME()
      write(6,*)
c       write(6,'(" CPU time in MPI loop is ",f20.10," s",/)')(tl2-tl1)
      write(6,'(" CPU time spent symmetrizing: ",f20.10," s")')(t2-t1)
      endif

      nncoef=nbasis

c     remove unconnected and intramolecular only terms if required
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      IF (my_id.eq.0) THEN
      IF (lnounc.or.lnointra) THEN

        do k=1,nterm
          itype(k)=0
        enddo

        DO m=1,nchan     ! loop over fragment channels, m
        nbadi=0
        nbadu=0
        ghi=0
        glo=natom
        do i=1,natom
        do j=i+1,natom
          g1=min(changroup(m,i),changroup(m,j))
          g2=max(changroup(m,i),changroup(m,j))
          ghi=max(g1,g2,ghi)
          glo=min(g1,g2,glo)
        enddo
        enddo
        do k=1,nterm     ! loop over terms, k
          do i=1,3
          do j=1,natom
            ugg(i,j)=0
          enddo
          enddo
          ii=0           ! loop over atom pairs i<j
          do i=1,natom
          do j=i+1,natom
            ii=ii+1
            g1=min(changroup(m,i),changroup(m,j))
            g2=max(changroup(m,i),changroup(m,j))
            if (g1.eq.g2) then
              ugg(1,g1)=ugg(1,g1)+ind(k,ii) ! u11
              do l=glo,ghi
              if (g1.ne.l) ugg(2,l)=ugg(2,l)+ind(k,ii) ! u22
              enddo
            endif
            if (g1.ne.g2) then
              ugg(3,g1)=ugg(3,g1)+ind(k,ii) ! u12
              ugg(3,g2)=ugg(3,g2)+ind(k,ii) ! u12
              do l=glo,ghi
              if (g1.ne.l.and.g2.ne.l) ugg(2,l)=ugg(2,l)+ind(k,ii) ! u22
              enddo
            endif
          enddo
          enddo

          ityp=0
          do i=glo,ghi   ! loop over groups
            if (ugg(1,i).eq.0.and.ugg(2,i).eq.0) ng=0
            if (ugg(1,i).gt.0.and.ugg(2,i).eq.0) ng=1
            if (ugg(1,i).eq.0.and.ugg(2,i).gt.0) ng=-1
            if (ugg(1,i).gt.0.and.ugg(2,i).gt.0) ng=2
            if (ng.eq.0.or.ugg(3,i).ne.0) then
              ityp=0 ! ok
            elseif (ng.eq.1.and.lnointra.and.ityp.eq.0) then
              ityp=1 ! intra only 
              nbadi=nbadi+1
            elseif (ng.eq.2.and.lnounc.and.ityp.eq.0) then
              ityp=2 ! unconnected
              nbadu=nbadu+1
            endif
          itype(ibasis(k)) = max(itype(ibasis(k)),ityp)
          enddo
        enddo

        nin=0
        nun=0
        do i=1,nbasis
          if (itype(i).eq.1) nin=nin+1
          if (itype(i).eq.2) nun=nun+1
          if (itype(i).eq.1.or.itype(i).eq.2) bad(nin+nun)=i
        enddo

        IF (my_id.eq.0) THEN
        write(6,*)
        if (m.eq.1) then
          write(6,*)"Analyzing fragment channel ",m
        if (lnointra)
     &    write(6,*)"Found ",nbadi," ( ",nin,
     &              ") intramolecular-only terms (groups)"
        if (lnounc)
     &    write(6,*)"Found ",nbadu," ( ",nun,
     &              ") unconnected terms (groups)"
        else
          write(6,*)"Updating group types for fragment channel ",m
        if (lnointra)
     &    write(6,*)"Found ",nbadi," intramolecular-only terms"
        if (lnointra)
     &    write(6,*)"The total number of intramolecular-only",
     &              "groups is now ",nin
        if (lnounc)
     &    write(6,*)"Found ",nbadu," unconnected terms"
        if (lnounc)
     &    write(6,*)"The total number of unconnected groups is now ",nun
        endif
        ENDIF

        ENDDO

        nt=nterm
        nb=nbasis

        ij=0
        do ii=1,nterm
          if (itype(ibasis(ii)).eq.0) then
          ij=ij+1
          ibasis(ij)=ibasis(ii)
          do k=1,npairs
            ind(ij,k)=ind(ii,k)
          enddo
          endif
        enddo
        nterm=ij

        do ii=1,nterm
          nx=ibasis(ii)
          do i=1,nin+nun
            if (nx.gt.bad(i)) ibasis(ii)=ibasis(ii)-1
          enddo
        enddo

        nbasis=nbasis-nin-nun

        IF (my_id.eq.0) THEN
        write(6,*)
        write(6,*)"Removing these groups results in the",
     & " following reductions:"
        write(6,'(a,i10,a,i10,a,i10,a,f15.3,a)')
     &             " Terms:  ",nt,"  - ",nt-nterm," = ",nterm," ("
     &             ,(dble(nterm)/dble(nt)*100.)," % )"
        write(6,'(a,i10,a,i10,a,i10,a,f15.3,a)')
     &             " Groups: ",nb,"  - ",nb-nbasis," = ",nbasis," ("
     &             ,(dble(nbasis)/dble(nb)*100.)," % )"
        ENDIF
      ENDIF

      nncoef=nbasis

      write(6,*)
      write(6,*)"Finished preparing the basis and writing it to",
     & " basis.dat"
      write(6,*)"This basis has ",nterm," (",nncoef,") terms (groups)"
      write(6,*)
      open(55,file="basis.dat")
      a=natom
      if (linter) a=natom1
      write(55,*)a,npairs,nncoef,nterm,
     & " ! atoms, atom pairs, groups/coefficients, terms"
      write(55,*)"   TERM   GROUP : EXPONENTS"
      if (lwrite(10)) then
        write(10,*)
        write(10,*)"Symmetrized Basis"
        if (.not.linter) then
!          write(10,'(23x,1000(a3,"- ",a3,2x))')
          write(10,'(23x,*(a3,"- ",a3,2x))')
     &       ((symb(i),symb(j),j=1+i,natom),i=1,natom)
!          write(10,'(22x,1000(i3," -",i3,2x))')((i,j,j=1+i,natom),
          write(10,'(22x,*(i3," -",i3,2x))')((i,j,j=1+i,natom),
     &       i=1,natom)
        else
!          write(10,'(23x,1000(a3,"- ",a3,2x))')
          write(10,'(23x,*(a3,"- ",a3,2x))')
     &       ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1)
!          write(10,'(22x,1000(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
          write(10,'(22x,*(i3," -",i3,2x))')((i,j,j=1+natom1,natom),
     &       i=1,natom1)
        endif
      endif
      do ii=1,nterm
!        write(55,'(2i8," : ",1000i4)')
        write(55,'(2i8," : ",*(i4))')
     &        ii,ibasis(ii),(ind(ii,j),j=1,npairs)
!        if (lwrite(10)) write(10,'(2i8,a,1000i4)')
        if (lwrite(10)) write(10,'(2i8,a,*(i4))')
     &     ii,ibasis(ii)," :",(ind(ii,j),j=1,npairs)
      enddo
      close(55)
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

ccc READ BASIS ccc
      ELSE
        IF (my_id.eq.0) THEN
        open(55,file="basis.dat")
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
      enddo

      do i=1,nterm
        arg=1.d0
        do j=1,npairs
          arg=arg*(r(j)**ind(i,j))
        enddo
        basis(ibasis(i))=basis(ibasis(i))+arg
!        do j=1,npairs
!          dbasisdr(ibasis(i),j)=dbasisdr(ibasis(i),j)   ! dV/dy * dy/dr
!     &                         -arg*dble(ind(i,j))
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

