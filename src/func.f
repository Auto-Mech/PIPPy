********************************************

      subroutine prepot

      implicit double precision(a-h,p-z)
      include 'param.inc'
      integer dgroup,a,b,ugg(maxatom,maxatom),g1,g2,ghi,glo,bad(maxterm)
      dimension iagroup(maxatom),itype(maxterm),
     &  ind(maxterm,maxpair),iatom(maxperm,maxatom),indord(maxterm),
     &  idum(maxatom),nngroup(maxatom),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  dbasisdr(maxterm,maxpair),rrr(maxdata,maxpair),
     &  index(maxatom,maxatom),ix(maxperm,maxpair),
     &  dgroup(maxchan,maxatom)
      character*2 symb(maxatom),dum
      logical lreadbasis,lnodisc,lnointra,linteronly,ldebug

      common/foox/rrr,nncoef,natom1,linteronly

      save npairs,nterm,ind,ibasis

c read
      ldebug=.true.
      ldebug=.false.
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

      if (lnodisc.or.lnointra.or.linteronly) then
      read(5,*)nchan
        do i=1,nchan
          read(5,*)(dgroup(i,k),k=1,natom)
        enddo
      endif
      if (linteronly) then
        if(nchan.ne.1) then
         write(6,*)"Mode = -1 requires 1 fragment channel with 2 groups"
           stop
        endif
        natom1=0        ! bimolecular only
        natom2=0
        do i=1,nchan
        do j=1,natom
          if (dgroup(i,j).eq.1) natom1=natom1+1
          if (dgroup(i,j).eq.2) natom2=natom2+1
          if (dgroup(i,j).ne.1.and.dgroup(i,j).ne.2) then
         write(6,*)"Mode = -1 requires 1 fragment channel with 2 groups"
            stop
          endif
        enddo
        enddo
        if (natom.ne.natom1+natom2) then
           print *,"Error specifying molecular groups"
           stop
        endif
      endif

      npairs=natom*(natom-1)/2                 ! Number of interatomic pairs
      if (linteronly) npairs = natom1*natom2

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

      if (natom.gt.maxatom) then
        print *,"natom (",natom,") > maxatom (",maxatom,")"
        stop
      endif

      if (npairs.gt.maxpair) then
        print *,"npairs (",npairs,") > maxpair (",maxpair,")"
        stop
      endif

      do i=1,nchan
      do j=1,natom
        if (dgroup(i,j).ne.1.and.dgroup(i,j).ne.2) then
          write(6,*)"This version of the code is limited to ",
     & "bimolecular fragments"
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
          print *,"npermute (",npermute,") > maxperm (",maxperm,")"
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

      print *
      print *,'Atom permutations = ',npermute
      do i=1,min(npermute,100) ! prints only first 100 permutations
       write(6,'(i7,a,1000i5)')i,":",(iatom(i,j),j=1,natom)
      enddo
      if (npermute.gt.100) write(6,*)" ** Truncating list **"
      write(6,*)

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
      if (npermute.gt.100.and.ldebug) write(6,*)" ** Truncating list **"
      if (ldebug) write(6,*)

c generate terms using user-supplied power constraints
      ii=1
      indtot=0
      do i=1,npairs
        ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
          print *,"number of terms (",ii,") > maxterm (",maxterm,")"
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

c symmetrize
      nbasis=0
      write(6,*)"Symmetrizing the expansion"
      DO ii=1,nterm
        if (mod(ii,100).eq.0) print *,"  step ",ii," of ",nterm
        ifail=0
        do i=ii-1,1,-1
          if (indord(ii).eq.indord(i)) then
          do j=1,npermute
            ifail=1
            do k=1,npairs
              if (ind(i,k).ne.ind(ii,ix(j,k))) then
                 ifail=0
                 exit
              endif
            enddo
            if (ifail.eq.1) go to 1010
          enddo
          endif
        enddo
 1010 continue
        if (ifail.eq.0) then
          nbasis=nbasis+1
          ibasis(ii)=nbasis
        else
          ibasis(ii)=ibasis(i)
        endif
      ENDDO
      write(6,*)

      nncoef=nbasis

c remove unconnected and intramolecular only terms if required
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
        print *,"Removing these groups results in the",
     & " following reductions:"
        print *,"Terms:  ",nt," --> ",nterm," ("
     &             ,(dble(nterm)/dble(nt)*100.)," % )"
        print *,"Groups: ",nb," --> ",nbasis," ("
     &             ,(dble(nbasis)/dble(nb)*100.)," % )"
        print *
      ENDIF

      nncoef=nbasis

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

ccc READ BASIS ccc
      ELSE
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

