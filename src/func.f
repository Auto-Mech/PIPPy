********************************************

      subroutine prepot

      implicit double precision(a-h,p-z)
      include 'param.inc'
      integer cc,comb,dtot,dpcheck,discind,pairpair,npow
      integer dgroup,dg1,dg2,ndg,dgtot,ntot,nintra,ndisc,a,b
      integer,allocatable :: itotal(:),itermtot(:)
      dimension iagroup(maxatom),!itermtot(maxterm),itotal(maxterm),
     &  ind(maxterm,maxpair),iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),idisc(maxatom,maxterm),
     &  idum2(maxperm,maxatom),idiscterm(maxatom,maxterm),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  dbasisdr(maxterm,maxpair),rrr(maxdata,maxpair),
     &  index(maxatom,maxatom),ix(maxperm,maxpair),
     &  dpcheck(maxpair),discind(maxatom,maxatom),
     &  pairpair(maxpair,maxpair),npow(maxpair),dgroup(maxatom,maxatom),
     &  ndg(maxatom),ndisc(maxatom),
     &  ndiscterm(maxatom),idisctotal(maxterm),idisctermtot(maxterm),
     &  iintratotal(maxterm),nintra(maxatom),nintraterm(maxatom),
     &  iintra(maxatom,maxterm),comb(2*maxterm),tcomb(2*maxterm),
     &  itmp(maxatom),itmp2(maxterm),itmp3(maxatom),itmp4(maxterm),
     &  iintraterm(maxatom,maxterm),iintratermtot(maxterm)
      character*2 symb(maxatom),dum
      logical lreadbasis,lreaddisc,lremintra,linteronly

      common/foox/rrr,nncoef,natom1,linteronly

      save npairs,nterms,ind,ibasis

      read(5,*)natom
      if (natom.gt.maxatom) then
        print *,"natom (",natom,") > maxatom (",maxatom,")"
        stop
      endif

      read(5,*)(symb(k),k=1,natom)
      read(5,*)(iagroup(k),k=1,natom)
      read(5,*)lreadbasis,ipow,ipowt,imode,nchans
      write(6,'(100("*"))')
      write(6,*)
      if (lreadbasis) then
      write(6,*)"Reading the PIP",ipow,ipowt," basis"
      else
      write(6,*)"Generating a PIP",ipow,ipowt,
     & " basis and writing it to basis.dat"
      endif
      if (imode.eq.0) then     ! use all terms
        lreaddisc = .false.  ! remove unconnected terms?
        lremintra = .false.  ! remove intramolecular-only terms?
        linteronly = .false. ! use only intermolecular terms?
        write(6,*)"Mode = 0: Using all bond distances"
      elseif (imode.eq.1) then ! remove unconnected and intra-only terms
        lreaddisc = .true.
        lremintra = .true.
        linteronly = .false.
        write(6,*)"Mode = 1: Using all bond distances and removing",
     & " unconnected and intramolecular-only terms"
      elseif (imode.eq.2) then ! remove unconnected terms
        lreaddisc = .true.
        lremintra = .false.
        linteronly = .false.
        write(6,*)"Mode = 2: Using all bond distances and removing",
     & " disconnected terms"
      elseif (imode.eq.3) then ! remove intramolecular terms
        lreaddisc = .false.
        lremintra = .true.
        linteronly = .false.
        write(6,*)"Mode = 3: Using all bond distances and removing",
     & " intramolecular-only terms"
      elseif (imode.eq.-1) then ! use intermolecular distances only
        lreaddisc = .false.
        lremintra = .false.
        linteronly = .true.
        write(6,*)"Mode = -1: Using only intermolecular bond distances"
      endif
      write(6,*)

      if (lreaddisc.or.lremintra.or.linteronly) then
        do i=1,nchans
          read(5,*)(dgroup(i,k),k=1,natom)
        enddo
      endif
      if (linteronly) then
        natom1=0        ! bimolecular only
        natom2=0
        do i=1,nchans
        do j=1,natom
          if (dgroup(i,j).eq.1) natom1=natom1+1
          if (dgroup(i,j).eq.2) natom2=natom2+1
        enddo
        enddo
        if (natom.ne.natom1+natom2) then
           print *,"Error specifying molecular groups"
           stop
        endif
      endif

      npairs=natom*(natom-1)/2                 ! Number of interatomic pairs
      if (linteronly) npairs = natom1*natom2

      if (npairs.gt.maxpair) then
        print *,"npairs (",npairs,") > maxpair (",maxpair,")"
        stop
      endif

      write(6,'(34("*  "))')
      write(6,*)
      write(6,*)"Atoms"
      write(6,'(1x,a8,100a5)')"Symbol",(symb(i),i=1,natom)
      write(6,'(a8,100i5)')"Index",(i,i=1,natom)
      write(6,'(a8,100i5)')"Group",(iagroup(i),i=1,natom)
      if (lreaddisc.or.lremintra.or.linteronly) then
      write(6,'(3x,a)')"Fragment Channels"
c      do loop over channels
      do i=1,nchans
        write(6,'(2x,i5,a,100i5)')i,":",(dgroup(i,j),j=1,natom)
      enddo
c      end do loop over channels
      endif

ccc GENERATE BASIS ccc
      IF (.not.lreadbasis) THEN
ccc GENERATE BASIS ccc

c     generate atom permutation lists
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
c      print *,"Group ",i," has ",(nperm1(i)-nperm0(i)+1)," permutations"
        ntmp=ntmp*(nperm1(i)-nperm0(i)+1)
      enddo
c      print *,"For a total of ",ntmp," permutations"

      npermute=0
      do while (.true.)
        npermute=npermute+1
        if (npermute.gt.maxperm) then
          print *,"npermute (",npermute,") > maxperm (",maxperm,")"
          print *,"NOTE: maxperm needs to be at least npermute + 1"
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
        if (ii.le.100) write(6,'(i7,a,1000i10)')
     &     ii,":",(ix(ii,iix),iix=1,npairs)
      enddo
      if (npermute.gt.100) write(6,*)" ** Truncating list **"
      write(6,*)

c generate terms using individual power constraints
      ii=1
      do i=1,npairs
        ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
          print *,"number of terms (",ii,") > maxterm (",maxterm,")"
          stop
        endif

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
      nterms=ii-1

c symmetrize
      nbasis=0
      write(6,*)"Symmetrizing the expansion"
      DO ii=1,nterms
        if (mod(ii,100).eq.0) print *,"  step ",ii," of ",nterms
        ifail=0
        do i=1,ii-1
          do j=1,npermute
            ifail=1
            do k=1,npairs
              if (ind(i,k).ne.ind(ii,ix(j,k))) ifail=0
            enddo
            if (ifail.eq.1) go to 1010
          enddo
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

      IF (lreaddisc) THEN

        DO m=1,nchans

        ndisc(m)=0
        ndiscterm(m)=0

        dgtot=0         ! number of DGs
        do i=1,natom
          ndg(i)=0      ! number of atoms in each DG
        enddo
        do i=1,npairs      
          dpcheck(i)=0   ! becomes 1 if pair matches DG conditions
          discind(i,1)=0
          discind(i,2)=0
        enddo

        ndp=0        ! number of unconnected pairs

        do i=1,natom
          ndg(dgroup(m,i))=ndg(dgroup(m,i))+1
          if (dgroup(m,i).gt.dgtot) dgtot=dgroup(m,i)
        enddo

ccc LOOP OVER MOLECULAR GROUPS
        DO dg1=1,dgtot-1
        if (ndg(dg1).le.1) then
          print *,"No unconnected groups possible"
        else
        DO dg2=dg1+1,dgtot
          if (ndg(dg2).le.1) then
            print *,"No unconnected groups possible"
          else
            print *,"Unconnected term list (channel ",m,")"
            write(6,'(26x,1000(a3,"- ",a3,2x))')
     &      ((symb(i),symb(j),j=1+i,natom),i=1,natom) 
            write(6,'(25x,1000(i3," -",i3,2x))')((i,j,j=1+i,natom),
     &            i=1,natom)
            do i=1,natom
              do j=i+1,natom
                ag=0
                dg=0
                if (iagroup(i).eq.iagroup(j)) ag=1
                if (dgroup(m,i).eq.dgroup(m,j)) dg=1
                if (dg.eq.1) then ! if dg = 1, all good
                  ndp=ndp+1
                  dpcheck(index(i,j))=1
                  discind(ndp,1)=i
                  discind(ndp,2)=j
                  do k=1,natom ! inner loop over atoms, find matching AGs
                    if (iagroup(k).eq.iagroup(i)) then
                      do l=k+1,natom
                        if (iagroup(l).eq.iagroup(j)) then
                          if (dpcheck(index(k,l)).eq.0) then
                            ndp=ndp+1
                            dpcheck(index(k,l))=1
                            discind(ndp,1)=k
                            discind(ndp,2)=l
                          endif
                        endif
                      enddo
                    elseif (iagroup(k).eq.iagroup(j)) then
                      do l=k+1,natom
                        if (iagroup(l).eq.iagroup(i)) then
                          if (dpcheck(index(k,l)).eq.0) then
                            ndp=ndp+1
                            dpcheck(index(k,l))=1
                            discind(ndp,1)=k
                            discind(ndp,2)=l
                          endif
                        endif
                      enddo
                    endif
                  enddo
                elseif (ag.eq.1) then ! dg = 0, ag = 1, also all good
                  ndp=ndp+1
                  dpcheck(index(i,j))=1
                  discind(ndp,1)=i
                  discind(ndp,2)=j
                endif
              enddo
            enddo
          endif
        ENDDO ! dg2
        endif
        ENDDO ! dg1

!        write(6,'(18x,100i10)')
!     &   (dpcheck(i),i=1,npairs)
      ! should give dpcheck(i)=1 for possible pairs
c
c       create pair-pair matrix
        do i=1,npairs
          npow(i)=0
          do j=1,npairs
            pairpair(i,j)=0
          enddo
        enddo

        do cc=1,ndp ! unconnected pairs
          i=discind(cc,1) ! pair 1, atom 1
          j=discind(cc,2) ! pair 1, atom 2
          do k=i+1,natom ! pair 2, atom 1
            do l=k+1,natom ! pair 2, atom 2
              if ((dpcheck(index(i,j)).eq.1).and.
     &            (dpcheck(index(k,l)).eq.1)) then
                if (j.ne.k.and.j.ne.l) then
                  pairpair(index(i,j),index(k,l))=1
                  pairpair(index(k,l),index(i,j))=1
                endif
              endif
            enddo
          enddo
        enddo

c TEST PAIRS
        do ii=1,nterms ! basis terms
          npx=0 ! number of pairs
          ipx=0 ! total power of pairs
          do i=1,npairs ! permutation pairs
            if (ind(ii,i).ne.0) then ! if there are powers of ii
              npx=npx+1
              ipx=ipx+i
              npow(npx)=i
            endif
          enddo

          if (npx.eq.2.and.pairpair(npow(1),npow(2)).eq.1) then ! pair is a monitored pair
            ndupe=0
            do j=1,ndisc(m)
              if (idisc(m,j).eq.ibasis(ii)) ndupe=1
            enddo
            ndiscterm(m)=ndiscterm(m)+1
            idiscterm(m,ndiscterm(m))=ii
!            if (ndupe.eq.0) then
              write(6,'(i8,"  (",i8,"):",i9,1000i10)')
     &            ii,ibasis(ii),(ind(ii,j),j=1,npairs)
            if (ndupe.eq.0) then
              ndisc(m)=ndisc(m)+1
              idisc(m,ndisc(m))=ibasis(ii)
            endif
          endif
        enddo ! ii=1,nterms
 
        print *
        print *,"Found ",ndiscterm(m)," (",ndisc(m),
     & ") unconnected terms (groups)"
        print *

        ENDDO ! nchans

c       Combine terms to remove from each channel
        if (nchans.ge.2) then
!          call combine(idisc(1,:),idisc(2,:),
!     &                 ndisc(1),ndisc(2),idisctotal,ndisctotal)
          call union(idisc(1,:),idisc(2,:),
     &                 ndisc(1),ndisc(2),itmp,ntmp)
          ndisctotal=ntmp
          do i=1,ndisctotal
            idisctotal(i)=itmp(i)
          enddo
          call union(idiscterm(1,:),idiscterm(2,:),
     &           ndiscterm(1),ndiscterm(2),itmp2,ntmp2)
          ndisctermtot=ntmp2
          do i=1,ndisctermtot
            idisctermtot(i)=itmp2(i)
          enddo
          if (nchans.ge.3) then
            DO m=1,nchans
!            call combine(idisctotal,idisc(m,:),
!     &                   ndisctotal,ndisc(m),itmp,ntmp)
              call union(idisctotal,idisc(m,:),
     &                   ndisctotal,ndisc(m),itmp,ntmp)
              ndisctotal=ntmp
              do i=1,ndisctotal
                idisctotal(i)=itmp(i)
              enddo
              call union(idisctermtot,idiscterm(m,:),
     &                   ndisctermtot,ndiscterm(m),itmp2,ntmp2)
              ndisctermtot=ntmp2
              do i=1,ndisctermtot
                idisctermtot(i)=itmp2(i)
              enddo
            ENDDO ! nchans
          endif
        else ! single fragment channel
          ndisctotal=ndisc(1)
          do i=1,ndisctotal
            idisctotal(i)=idisc(1,i)
          enddo
          ndisctermtot=ndiscterm(1)
          do i=1,ndisctermtot
            idisctermtot(i)=idiscterm(1,i)
          enddo
        endif

        print *
        write(6,'(a28,i8,a2,i8,a37)'),"Across all channels, found ",
     &         ndisctermtot,"(",ndisctotal,
     &         ") common unconnected terms (groups)"
!        do ii=1,ndisctermtot
!          write(6,'(i8,"  (",i8,"):",i9,1000i10)')
!     &         idisctermtot(ii),ibasis(idisctermtot(ii)),
!     &         (ind(idisctermtot(ii),k),k=1,npairs)
!        enddo
        print *

      ENDIF ! lreaddisc
c end remove unconnected terms

cccccccccccccc INTRAMOLECULAR TERMS
      IF (lremintra) THEN

        DO m=1,nchans

        nintra(m)=0
        nintraterm(m)=0

        print *,'Intramolecular-only term list'
        write(6,'(26x,1000(a3,"- ",a3,2x))')
     &   ((symb(i),symb(j),j=1+i,natom),i=1,natom) 
        write(6,'(25x,1000(i3," -",i3,2x))')((i,j,j=1+i,natom),
     &   i=1,natom)

        DO ii=1,nterms
          ifail=1
          g1=0
          if (ii.eq.1) then
            ifail=0 ! keep constant term
            goto 501
          endif
          do i=1,natom
          do j=i+1,natom
            if ((ind(ii,index(i,j)).ne.0).and.(ifail.eq.1)) then
              if (g1.eq.1) then
            if ((dgroup(m,i).ne.tmpgrp).or.(dgroup(m,j).ne.tmpgrp)) then
                  ifail=0
                  goto 501
                endif
              endif
              if (dgroup(m,i).ne.dgroup(m,j)) then
                ifail=0
                goto 501
              elseif (g1.eq.0) then
                tmpgrp=dgroup(m,i)
                g1=1
              endif
            endif
          enddo
          enddo
  501     continue

          if (ifail.eq.1) then ! intramolecular term
            ndupe=0
            do jj=1,nintra(m)
              if (iintra(m,jj).eq.ibasis(ii)) ndupe=1
            enddo
            nintraterm(m)=nintraterm(m)+1
            iintraterm(m,nintraterm(m))=ii
            if (ndupe.eq.0) then
              nintra(m)=nintra(m)+1
              iintra(m,nintra(m))=ibasis(ii)
              write(6,'(i8,"  (",i8,"):",i9,1000i10)')
     &         ii,ibasis(ii),(ind(ii,k),k=1,npairs)
            endif
          endif
        ENDDO

        print *
        print *,"Found ",nintraterm(m)," (",nintra(m),
     &   ") intramolecular-only terms (groups)"
        print *

        ENDDO ! nchans

        ! condense iintra into single array, no duplicate terms
        if (nchans.ge.2) then
!          call combine(iintra(1,:),iintra(2,:),
!     &                 nintra(1),nintra(2),iintratotal,nintratotal)
          call union(iintra(1,:),iintra(2,:),
     &                 nintra(1),nintra(2),itmp3,ntmp3)
          nintratotal=ntmp3
          do i=1,nintratotal
            iintratotal(i)=itmp3(i)
          enddo
          call union(iintraterm(1,:),iintraterm(2,:),
     &               nintraterm(1),nintraterm(2),itmp4,ntmp4)
          nintratermtot=ntmp4
          do i=1,nintratermtot
            iintratermtot(i)=itmp4(i)
          enddo
          if (nchans.ge.3) then
            DO m=1,nchans
!              call combine(iintratotal,iintra(m,:),
!     &                   nintratotal,nintra(m),itmp,ntmp)
              call union(iintratotal,iintra(m,:),
     &                   nintratotal,nintra(m),itmp3,ntmp3)
              nintratotal=ntmp3
              do i=1,nintratotal
                iintratotal(i)=itmp3(i)
              enddo
              call union(iintratermtot,iintraterm(m,:),
     &                   nintratermtot,nintraterm(m),itmp4,ntmp4)
              nintratermtot=ntmp4
              do i=1,nintratermtot
                iintratermtot(i)=itmp4(i)
              enddo
            ENDDO ! nchans
          endif
        else ! single channel
          nintratotal=nintra(1)
          do i=1,nintratotal
            iintratotal(i)=iintra(1,i)
          enddo
          nintratermtot=nintraterm(1)
          do i=1,nintratermtot
            iintratermtot(i)=iintraterm(1,i)
          enddo
        endif

        print *
        write(6,'(a28,i8,a2,i8,a44)'),"Across all channels, found ",
     &         nintratermtot,"(",nintratotal,
     &         ") common intramolecular-only terms (groups)"
!        do ii=1,nintratermtot
!          write(6,'(i8,"  (",i8,"):",i9,1000i10)')
!     &         iintratermtot(ii),ibasis(iintratermtot(ii)),
!     &         (ind(iintratermtot(ii),k),k=1,npairs)
!        enddo
        print *

      ENDIF

cccccccc REMOVE DISCONNNECTED AND/OR INTRAMOLECULAR-ONLY TERMS
      IF (lreaddisc.and.lremintra) THEN
        if (ndisctotal.eq.0) then
          ntot=nintratotal
          allocate(itotal(ntot))
          do i=1,ntot
            itotal(i)=iintratotal(i)
          enddo
          ntermtot=nintratermtot
          allocate(itermtot(ntermtot))
          do i=1,ntermtot
            itermtot(i)=iintratermtot(i)
          enddo
        elseif (nintratotal.eq.0) then
          ntot=ndisctotal
          allocate(itotal(ntot))
          do i=1,ntot
            itotal(i)=idisctotal(i)
          enddo
          ntermtot=ndisctermtot
          allocate(itermtot(ntermtot))
          do i=1,ntermtot
            itermtot(i)=idisctermtot(i)
          enddo
        else
          ! Combine Group Lists, Remove Duplicates
!        call combine(idisctotal,iintratotal,
!     &               ndisctotal,nintratotal,itotal,ntot)
!        call combine(idisctermtot,iintratermtot,
!     &               ndisctermtot,nintratermtot,itermtot,ntermtot)
          j=ndisctotal+1
          do i=1,ndisctotal
            comb(i)=idisctotal(i)
          enddo
          do i=1,nintratotal
            if (.not.any(iintratotal(i)==comb(1:ndisctotal))) then
              comb(j)=iintratotal(i)
              j=j+1
            endif
          enddo
          ntot=j-1
          allocate(itotal(ntot))
          do i=1,ntot
            itotal(i)=comb(i)
          enddo
          ! Combine Term Lists
          j=ndisctermtot+1
          do i=1,ndisctermtot
            tcomb(i)=idisctermtot(i)
          enddo
          do i=1,nintratermtot
            if (.not.any(iintratermtot(i)==tcomb(1:ndisctermtot))) then
              tcomb(j)=iintratermtot(i)
              j=j+1
            endif
          enddo
          ntermtot=j-1
          allocate(itermtot(ntermtot))
          do i=1,ntermtot
            itermtot(i)=tcomb(i)
          enddo
        endif
        !!!!!!!!!!!!!!!!!!!!!
      ELSEIF (lreaddisc) THEN
        ntot=ndisctotal
        allocate(itotal(ntot))
        do i=1,ntot
          itotal(i)=idisctotal(i)
        enddo
        ntermtot=ndisctermtot
        allocate(itermtot(ntermtot))
        do i=1,ntermtot
          itermtot(i)=idisctermtot(i)
        enddo
      ELSEIF (lremintra) THEN
        ntot=nintratotal
        allocate(itotal(ntot))
        do i=1,ntot
          itotal(i)=iintratotal(i)
        enddo
        ntermtot=nintratermtot
        allocate(itermtot(ntermtot))
        do i=1,ntermtot
          itermtot(i)=iintratermtot(i)
        enddo
      ENDIF

      call sort(itotal,ntot)
!      call sort(itermtot,ntermtot)

      print *
      write(6,'(a28,i8,a2,i8,a27)'),"Across all channels, found ",
     &         ntermtot,"(",ntot,
     &         ") terms (groups) to remove"
      print *
 
      IF (lreaddisc.or.lremintra) THEN
        do ii=nterms,1,-1
          ibad=0
          do j=1,ntot
            if (ibasis(ii).eq.itotal(j)) ibad=1
          enddo
          if (ibad.eq.1) then
            nterms=nterms-1
            do jj=ii,nterms
              ibasis(jj)=ibasis(jj+1)
              do k=1,npairs
                ind(jj,k)=ind(jj+1,k)
              enddo
            enddo
          endif
        enddo ! ii=nterms,1,-1

        do ii=1,nterms
          nx=ibasis(ii)
          do i=1,ntot
            if (nx.gt.itotal(i)) ibasis(ii)=ibasis(ii)-1
          enddo
        enddo ! ii=1,nterms
      ENDIF

c      print *
c      print *,"Basis # (Group):  Powers"

c      do ii=1,nterms
c        write(6,'(i5,"  (",i5,"):",100i8)')
c     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
c      enddo
c      write(6,*)

      nbasis=nbasis-ntot
      nncoef=nbasis

      write(6,'(34("*  "))')
      write(6,*)
      write(6,*)"Finished preparing the basis and writing it to",
     & " basis.dat"
      write(6,*)"This basis has ",nterms," (",nncoef,") terms (groups)"
      write(6,*)
      open(55,file="basis.dat")
      a=natom
      if (linteronly) a=natom1
      write(55,*)a,npairs,nncoef,nterms,
     & " ! atoms, atom pairs, groups/coefficients, terms"
      write(55,*)"   TERM   GROUP : EXPONENTS"
      do ii=1,nterms
        write(55,'(2i8," : ",1000i8)')
     &        ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)

ccc READ BASIS ccc
      ELSE
        open(55,file="basis.dat")
        a=natom
        if (linteronly) a=natom1
        read(55,*)a,npairs,nncoef,nterms
        read(55,*)
        do i=1,nterms
          read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
        enddo
        close(55)
      write(6,*)
      write(6,*)"Reading the basis from basis.dat"
      write(6,*)"This basis has ",nterms," (",nncoef,") terms (groups)"
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

      do i=1,nterms
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

      subroutine sort(a,n)

      integer n,a(n),temp,k,l

      do k=1,n-1
        do l=k+1,n
          if (a(k).gt.a(l)) then
            temp = a(k)
            a(k) = a(l)
            a(l) = temp
          endif
        enddo
      enddo
      return

      end subroutine

***************************************************

      subroutine combine(a,b,m,n,c,z)

      integer m,n,a(m),b(n),comb(m+n),z,c(m+n)

      j=m+1

      do i=1,m
        comb(i)=a(i)
      enddo

      do i=1,n
        if (.not.any(b(i).eq.comb(1:m))) then
          comb(j)=b(i)
          j=j+1
        endif
      enddo

      z=j-1

!      allocate(itotal(ntot))
      do i=1,z
        c(i)=comb(i)
      enddo

      return

      end subroutine

***************************************************

      subroutine union(a,b,m,n,comb,k)

      integer m,n,a(m),b(n),k,comb(m+n)

      k=0

      if (m.gt.n) then
        do i=1,n
          comb(i)=0
        enddo
        do i=1,n
          if (any(b(i).eq.a)) then
            k=k+1
            comb(k)=b(i)
          endif
        enddo
      else
        do i=1,m
          comb(i)=0
        enddo
        do i=1,m
          if (any(b(i).eq.a)) then
            k=k+1
            comb(k)=b(i)
          endif
        enddo
      endif

      return

      end subroutine

