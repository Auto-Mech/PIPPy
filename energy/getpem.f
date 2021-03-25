      subroutine getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)

c     This subroutine handles all potential calls.
c     The adiabatic matrix PEMA, its gradients GPEMA,
c     the diabatic matrix PEMD, its gradients GPEMD,
c     and the nonadiabatic coupling vector DVEC are returned
c     PEMA is a diagonal matrix, so it's a 1-D vector, with the
c     index running from 1-NSURFT.
c     PEMD is 2-D, with indices running from 1-NSURFT.
c     For GPEMA and GPEMD, the first index represents x,y,z, the
c     second represents the atom index, the third (and fourth for
c     PEMD) label the electronic states and run from 1-NSURFT.
c     NCLU is the number of atoms in this call, this number may not be
c     the same as NAT, but must be less than NAT

      implicit none
      include 'param.f'
      include 'c_sys.f'

      integer i,j,k,l,i1,i2,j1,j2
c TEMP AJ
      integer iskip
      double precision rad,gtmp(3,mnat),gtmp2(3,mnat)

c     input
      integer nclu
      double precision xx(3,mnat)
      character*2 symb(mnat)

c     output
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & pema(mnsurf),gpema(3,mnat,mnsurf),dvec(3,mnat,mnsurf,mnsurf)

c local (for 3V-2 interface)
      double precision r(3),du11(3),du12(3),du22(3),v1,v2,u11,u12,u22,
     &  cs1,sn1,cc22(2,2),tmp

c local (for HO-MM-1 and HE-MM-1 interfaces)
      double precision dx(nat),dy(nat),dz(nat),x(nat),y(nat),z(nat),v

c local (for multistate HO-MM-1)
      double precision ap(mnsurf*(mnsurf+1)/2),work(3*mnsurf),
     & cc(mnsurf,mnsurf),dum
      integer info

      nsurft = 1 ! TEMP DM

c ######################################################################
c     POTLIB HO-MM-1 interface
      if (potflag.eq.0) then
c ######################################################################
c     Single surface

      if (nsurft.ne.1) then
      write(6,*)"STOP in GETPEM"
      write(6,*)"NSURFT = ",nsurft," is not allowed for POTFLAG = 0"
      stop
      endif

      do i=1,nclu
      x(i) = xx(1,i)
      y(i) = xx(2,i)
      z(i) = xx(3,i)
      enddo
      call pot(x,y,z,v,dx,dy,dz,nclu,nclu)
      pema(1) = v
      pemd(1,1) = v
      do i=1,nclu
      gpema(1,i,1)=dx(i)
      gpema(2,i,1)=dy(i)
      gpema(3,i,1)=dz(i)
      gpemd(1,i,1,1)=dx(i)
      gpemd(2,i,1,1)=dy(i)
      gpemd(3,i,1,1)=dz(i)
      dvec(1,i,1,1)=0.d0
      dvec(2,i,1,1)=0.d0
      dvec(3,i,1,1)=0.d0
      enddo

c ######################################################################
c     POTLIB HE-MM-1 interface
      elseif (potflag.eq.2) then
c ######################################################################
c     Single surface

      if (nsurft.ne.1) then
      write(6,*)"STOP in GETPEM"
      write(6,*)"NSURFT = ",nsurft," is not allowed for POTFLAG = 2"
      stop
      endif

      do i=1,nclu
      x(i) = xx(1,i)
      y(i) = xx(2,i)
      z(i) = xx(3,i)
      enddo
      call pot(symb,x,y,z,v,dx,dy,dz,nclu,mnat)
      pema(1) = v
      pemd(1,1) = v
      do i=1,nclu
      gpema(1,i,1)=dx(i)
      gpema(2,i,1)=dy(i)
      gpema(3,i,1)=dz(i)
      gpemd(1,i,1,1)=dx(i)
      gpemd(2,i,1,1)=dy(i)
      gpemd(3,i,1,1)=dz(i)
      dvec(1,i,1,1)=0.d0
      dvec(2,i,1,1)=0.d0
      dvec(3,i,1,1)=0.d0
      enddo

c ######################################################################
c     POTLIB 3-2V interface
      elseif (potflag.eq.1) then
c ######################################################################
c     Triatomic, two-surface interface

      if (nclu.ne.3) then
      write(6,*)"STOP in GETPEM"
      write(6,*)"NAT = ",nclu," is not allowed for POTFLAG = 1"
      stop
      endif

      if (nsurft.ne.2) then
      write(6,*)"STOP in GETPEM"
      write(6,*)"NSURFT = ",nsurft," is not allowed for POTFLAG = 1"
      stop
      endif

c     r(1) = A-B
c     r(2) = B-C
c     r(3) = A-C
c     im = 1 = AB+C
c     im = 2 = BC+A
c     im = 3 = CA+B
c     transform to internals
        r(1)=0.d0
        r(2)=0.d0
        r(3)=0.d0
      do i=1,3
        r(1) = r(1) + (xx(i,1)-xx(i,2))**2
        r(2) = r(2) + (xx(i,2)-xx(i,3))**2
        r(3) = r(3) + (xx(i,1)-xx(i,3))**2
      enddo
      do i=1,3
        r(i)=dsqrt(r(i))
      enddo
      call pot(r,u11,du11,1,1)
      call pot(r,u12,du12,1,2)
      call pot(r,u22,du22,1,3)

c     save diabatic info
      pemd(1,1) = u11
      pemd(1,2) = u12
      pemd(2,1) = u12
      pemd(2,2) = u22
      do i=1,3
c       convert from internals to cartesian
        gpemd(i,1,1,1) = du11(1)/r(1)*(xx(i,1)-xx(i,2))
     *                  +du11(3)/r(3)*(xx(i,1)-xx(i,3))
        gpemd(i,2,1,1) = du11(1)/r(1)*(xx(i,2)-xx(i,1))
     *                  +du11(2)/r(2)*(xx(i,2)-xx(i,3))
        gpemd(i,3,1,1) = du11(2)/r(2)*(xx(i,3)-xx(i,2))
     *                  +du11(3)/r(3)*(xx(i,3)-xx(i,1))
        gpemd(i,1,1,2) = du12(1)/r(1)*(xx(i,1)-xx(i,2))
     *                  +du12(3)/r(3)*(xx(i,1)-xx(i,3))
        gpemd(i,2,1,2) = du12(1)/r(1)*(xx(i,2)-xx(i,1))
     *                  +du12(2)/r(2)*(xx(i,2)-xx(i,3))
        gpemd(i,3,1,2) = du12(2)/r(2)*(xx(i,3)-xx(i,2))
     *                  +du12(3)/r(3)*(xx(i,3)-xx(i,1))
        gpemd(i,1,2,1) = gpemd(i,1,1,2)
        gpemd(i,2,2,1) = gpemd(i,2,1,2)
        gpemd(i,3,2,1) = gpemd(i,3,1,2)
        gpemd(i,1,2,2) = du22(1)/r(1)*(xx(i,1)-xx(i,2))
     *                  +du22(3)/r(3)*(xx(i,1)-xx(i,3))
        gpemd(i,2,2,2) = du22(1)/r(1)*(xx(i,2)-xx(i,1))
     *                  +du22(2)/r(2)*(xx(i,2)-xx(i,3))
        gpemd(i,3,2,2) = du22(2)/r(2)*(xx(i,3)-xx(i,2))
     *                  +du22(3)/r(3)*(xx(i,3)-xx(i,1))
      enddo

c     compute adiabatic info
c     diagonalize the 2x2
      call dlaev2(u11,u12,u22,v2,v1,cs1,sn1)

c     phase convention (this part depends on how DLAEV2 works)
      if (v2.ge.v1) then
        cc22(1,1) = -sn1
        cc22(2,1) = cs1
        cc22(1,2) = cs1
        cc22(2,2) = sn1
      else
        tmp = v1
        v1 = v2
        v2 = tmp
        cc22(1,1) = cs1
        cc22(2,1) = sn1
        cc22(1,2) = -sn1
        cc22(2,2) = cs1
      endif

      pema(1) = v1
      pema(2) = v2

c     gradients of adiabatic surfaces v1 and v2
c     this bit of code can handle an arbitrary number of surfaces
c     if cc is defined properly
      do j1 = 1, 3
      do j2 = 1, nclu
      do i = 1, nsurft
        gpema(j1,j2,i)=(cc22(1,i)**2)*gpemd(j1,j2,1,1)
        do k = 2, nsurft
          gpema(j1,j2,i)=gpema(j1,j2,i)+(cc22(k,i)**2)*gpemd(j1,j2,k,k)
          do l = 1, k-1
            gpema(j1,j2,i)=gpema(j1,j2,i)
     &                     +2.d0*cc22(k,i)*cc22(l,i)*gpemd(j1,j2,k,l)
          enddo
        enddo
      enddo
      enddo
      enddo

c compute d
      do i1 = 1,3
      do i2 = 1,nclu
        do i = 2,nsurft
          do j = 1, i-1
            dvec(i1,i2,i,j) = 0.d0
            do k = 1, nsurft
            do l = 1, nsurft
              dvec(i1,i2,i,j) = dvec(i1,i2,i,j) 
     &                        + cc22(k,i)*cc22(l,j)*gpemd(i1,i2,k,l)
            enddo
            enddo
            if ((pema(j) - pema(i)) .ne. 0.0d0) then
              dvec(i1,i2,i,j) = dvec(i1,i2,i,j)/(pema(j)-pema(i))
            else
              dvec(i1,i2,i,j) = 0.0d0
            endif
            dvec(i1,i2,j,i) = -dvec(i1,i2,i,j)
          enddo
        enddo
      enddo
      enddo

c ######################################################################
c     HE-MM-1 interface with multiple diabatic surfaces
      elseif (potflag.eq.3) then
c ######################################################################

      do i=1,nclu
      x(i) = xx(1,i)
      y(i) = xx(2,i)
      z(i) = xx(3,i)
      enddo

      call pot(symb,x,y,z,pemd,gpemd,nclu,mnat,nsurft,mnsurf)

c     diagonalize the diabatic matrix
      do i=1,nsurft
      do j=i,nsurft
        ap(i+(j-1)*j/2)=pemd(i,j)
      enddo
      enddo
      call dspev( 'v','u',nsurft,ap,pema,cc,mnsurf,work,info )

c HACK HACK
c HACK HACK
c HACK HACK
      if (1.eq.2) then
c     compute adiabatic info
c     diagonalize the 2x2
      u11=pemd(1,1)
      u12=pemd(1,2)
      u22=pemd(2,2)
      call dlaev2(u11,u12,u22,v2,v1,cs1,sn1)

c     phase convention (this part depends on how DLAEV2 works)
      if (v2.ge.v1) then
        cc22(1,1) = -sn1
        cc22(2,1) = cs1
        cc22(1,2) = cs1
        cc22(2,2) = sn1
      else
        tmp = v1
        v1 = v2
        v2 = tmp
        cc22(1,1) = cs1
        cc22(2,1) = sn1
        cc22(1,2) = -sn1
        cc22(2,2) = cs1
      endif

c      print *,"a",cc(1,1),cc(1,2),cc(2,1),cc(2,2)
      rad = dsqrt(0.25d0*(u11-u22)**2+u12**2)
      v1 = 0.5d0*(u11+u22)-rad
      v2 = 0.5d0*(u11+u22)+rad

c      print *,"u",u11*autoev,u22*autoev
c      print *,"v",v1*autoev,v2*autoev
c      print *,"v",pema(1)*autoev,pema(2)*autoev
      v1 = 0.d0
      v2 = 0.d0
      do i=1,nsurft
      do j=1,nsurft
       v1 = v1 + cc(i,1)*pemd(i,j)*cc(j,1)
       v2 = v2 + cc(i,2)*pemd(i,j)*cc(j,2)
      enddo
      enddo
c      print *,"vv",v1*autoev,v2*autoev

c     gradients of adiabatic surfaces v1 and v2
      do j1 = 1, 3
      do j2 = 1, nclu
         gtmp(j1,j2) = 0.5d0*(gpemd(j1,j2,1,1)+gpemd(j1,j2,2,2))
     &    - 0.5d0*(0.5d0*(gpemd(j1,j2,1,1)-gpemd(j1,j2,2,2))
     *     *(pemd(1,1)-pemd(2,2))
     &     +2.d0*gpemd(j1,j2,1,2)*pemd(1,2))/rad
       gtmp2(j1,j2) = cc(1,1)**2*gpemd(j1,j2,1,1)
     &         +cc(1,1)*cc(2,1)*gpemd(j1,j2,1,2)
     &    +cc(1,1)*cc(2,1)*gpemd(j1,j2,1,2)
     &    +cc(2,1)**2*gpemd(j1,j2,2,2)
      enddo
      enddo

      endif
c HACK HACK
c HACK HACK
c HACK HACK

      do j1 = 1, 3
      do j2 = 1, nclu
      do i = 1, nsurft
         gpema(j1,j2,i)=0.d0
        do k = 1, nsurft
          do l = 1, nsurft
            gpema(j1,j2,i)=gpema(j1,j2,i)
     &                     +cc(k,i)*gpemd(j1,j2,k,l)*cc(l,i)
c     &                     +cc(k,i)*cc(i,l)*gpemd(j1,j2,k,l)
          enddo
        enddo
      enddo
      enddo
      enddo
c      write(6,99)((gpema(j1,j2,1),j1=1,3),j2=1,nclu)
c      write(6,99)((gtmp(j1,j2),j1=1,3),j2=1,nclu)
c      write(6,99)((gtmp2(j1,j2),j1=1,3),j2=1,nclu)
 99   format(100e11.3)

c compute d
      do i1 = 1,3
      do i2 = 1,nclu
        do i = 1,nsurft
          dvec(i1,i2,i,i) = 0.d0
        enddo
        do i = 2,nsurft
          do j = 1, i-1
            dvec(i1,i2,i,j) = 0.d0
            do k = 1, nsurft
            do l = 1, nsurft
              dvec(i1,i2,i,j) = dvec(i1,i2,i,j)
     &                        + cc(k,i)*cc(l,j)*gpemd(i1,i2,k,l)
c     &                        + cc(k,i)*cc(j,l)*gpemd(i1,i2,k,l)
            enddo
            enddo
            if ((pema(j) - pema(i)) .ne. 0.0d0) then
              dvec(i1,i2,i,j) = dvec(i1,i2,i,j)/(pema(j)-pema(i))
            else
              dvec(i1,i2,i,j) = 0.0d0
            endif
c AJ hack for CO2 PES
c            if (i.eq.3.and.j.eq.2) dvec(i1,i2,i,j)=0.d0
            dvec(i1,i2,j,i) = -dvec(i1,i2,i,j)
          enddo
        enddo
      enddo
      enddo

c ######################################################################
c     HE-MM-1 interface with multiple adiabatic surfaces
      elseif (potflag.eq.4) then
c ######################################################################

      do i=1,nclu
      x(i) = xx(1,i)
      y(i) = xx(2,i)
      z(i) = xx(3,i)
      enddo

      call pot(symb,x,y,z,pema,gpema,dvec,nclu,mnat,nsurft,mnsurf)

c     set diabatic energies to zero. these shouldn't be used
      do i=1,nsurft
      do j=i,nsurft
        pemd(i,j)=0.d0
      do j1 = 1, 3
      do j2 = 1, nclu
        gpemd(j1,j2,i,j)=0.d0
      enddo
      enddo
      enddo
      enddo

c ######################################################################
      else
      write(6,*)"Cant have POTFLAG = ",potflag," in GETPEM"
      stop
      endif
c ######################################################################

c     ADJUST ZERO OF ENERGY
      do i=1,nsurft
      pema(i)=pema(i)-ezero
      pemd(i,i)=pemd(i,i)-ezero
      enddo

      return

      end
