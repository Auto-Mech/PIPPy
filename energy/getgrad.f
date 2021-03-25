      subroutine getgrad(xx,pp,nsurf,v,cre,cim,gv,gcre,gcim,nclu,
     &  phop,dmag,dvec,pem,gpem,phase)

c Calls GETPEM once, and then computes several important variables:
c V = Potential energy
c GV = Gradient of the potential energy
c DMAG = Magnitude of the nonadiabatic coupling vector (for CSDM)
c GCRE and GCIM = Time derivatives of the real and imaginary parts
c   of the electronic variables, including the DM terms for the
c   SCDM and CSDM methods.
c PHOP = Hopping probability (for surface hoppping) or switching
c   probability (for SCDM and CSDM) divided by the stepsize.
c PEM = Array of adiabatic or diabatic potential energies.
c GPEM = Gradients of the adiabatic or diabatic potential energies.

      implicit none
      include 'param.f'
      include 'c_sys.f'

c temp
      double precision taup(3),taup2(3)
      common/tauprint/taup,taup2

      integer i,j,k,l,m,nsurf,ii,nclu,i2,j2
      double precision xx(3,mnat),pp(3,mnat),gv(3,mnat),v,
     & gcim(2*mnsurf),gcre(2*mnsurf),cre(2*mnsurf),cim(2*mnsurf)
      double precision pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),gpem(3,mnat,mnsurf),
     & rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf),vdotd(mnsurf,mnsurf),
     & phop(mnsurf),bb,pem(mnsurf),phase(mnsurf)

      double precision ppvib(3,mnat),svec(3,mnat,mnsurf),
     & Cparam,E0param,es,taui(mnsurf),ps(mnsurf),tmp,dum,dvecmag,pdotd,
     & u2b2(2,2),gu2b2(3,mnat,2,2),dvec2b2(3,mnat),gcred(mnsurf),
     & gcimd(mnsurf),vd(mnsurf),pdeco(3,mnat),smagms,
     & rhorc(mnsurf,mnsurf),rhoic(mnsurf,mnsurf),tmpr,
     & tmpi,dsum2,dmag,sntmp,cstmp,rhotmp,
     & rhorc2(mnsurf,mnsurf),rhoic2(mnsurf,mnsurf)

c msx search variables
       integer nupper
       double precision x1mag,x1(3,mnat),gupp(3,mnat),
     &   dot1,egap,gf(3,mnat),gg(3,mnat),msx_grad(3,mnat)

      print *,"xx = "
!      print *,xx
      print *,((xx(i,j),i=1,6),j=1,3)!,((pp(i,j),i=1,3),j=1,3)
c ZERO
      do i=1,nsurft
        phop(i) = 0.d0
      enddo

c GET ENERGIES
      call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symbol)

c POTENTIAL ENERGIES AND GRADIENTS
!      IF (METHFLAG.EQ.0.OR.METHFLAG.EQ.1.OR.METHFLAG.EQ.5) THEN
c     single surface propagation or surface hopping calculation

      if (repflag.eq.0) then
c       adiabatic
        v = pema(nsurf)
        do i=1,3
        do j=1,nclu
          gv(i,j)=gpema(i,j,nsurf)
        enddo
        enddo
      else if (repflag.eq.1) then
c       diabatic
        v = pemd(nsurf,nsurf)
        do i=1,3
        do j=1,nclu
          gv(i,j)=gpemd(i,j,nsurf,nsurf)
        enddo
        enddo
      else
        write(6,*)"REPFLAG = ",repflag," in GETPOT"
        stop
      endif

!      ELSE

!      write(6,*)"METHFLAG = ",methflag," is not allowed in GETPOT"
!      stop

!      ENDIF

c TIME-DERIVATIVES OF THE ELECTRONIC COORDINATES
c     All methods, DM methods add terms later
!      if (nsurft.eq.1) then
        gcre(1) = 0.d0
        gcim(1) = 0.d0
!      else
!        if (repflag.eq.1) then
!c         diabatic rep
!          do i=1,nsurft
!            gcim(i) = 0.d0
!            gcre(i) = 0.d0
!            do j=1,nsurft
!c integrate whole coefficient
!c              gcre(i) = gcre(i) - cim(j)*pemd(i,j)
!c              gcim(i) = gcim(i) + cre(j)*pemd(i,j)
!c end
!c integrate phase angle separately
!              if (i.ne.j) then
!                tmp = phase(j)-phase(i)
!                sntmp = dsin(tmp)
!                cstmp = dcos(tmp)
!                tmpr = -(sntmp*cre(j)-cstmp*cim(j))*pemd(i,j)
!                tmpi = -(sntmp*cim(j)+cstmp*cre(j))*pemd(i,j)
!                gcre(i) = gcre(i) + tmpr
!                gcim(i) = gcim(i) + tmpi
!                if (i.eq.nsurf) then
!                  phop(j) = -2.d0*(tmpr*cre(i)+tmpi*cim(i))
!     &               /(cre(i)**2+cim(i)**2)
!                endif 
!              endif 
!c end
!            enddo
!          enddo
!        elseif (repflag.eq.0) then
!c         adiabatic rep
!c         compute velocity-dot-dvec
!          do 10 k=1,nsurft
!          do 10 l=1,nsurft
!          vdotd(k,l) = 0.d0
!          do 10 i=1,3
!          do 10 j=1,nclu
!            vdotd(k,l) = vdotd(k,l) + dvec(i,j,k,l)*pp(i,j)/mm(j)
!  10      continue
!          do i=1,nsurft
!c integrate whole coefficient
!c            gcre(i) =  cim(i)*pema(i)
!c            gcim(i) = -cre(i)*pema(i)
!c end
!c integrate phase angle separately
!            gcre(i) = 0.d0
!            gcim(i) = 0.d0
!            do j=1,nsurft
!c end
!c integrate whole coefficient
!c              gcre(i) =  gcre(i) - cre(j)*vdotd(i,j)
!c              gcim(i) =  gcim(i) - cim(j)*vdotd(i,j)
!c end
!c integrate phase angle separately
!              if (i.ne.j) then
!                tmp = phase(j)-phase(i)
!                sntmp = dsin(tmp)
!                cstmp = dcos(tmp)
!                tmpr = -(cstmp*cre(j)+sntmp*cim(j))*vdotd(i,j)
!                tmpi = -(cstmp*cim(j)-sntmp*cre(j))*vdotd(i,j)
!                gcre(i) = gcre(i) + tmpr
!                gcim(i) = gcim(i) + tmpi
!                if (i.eq.nsurf) then
!                  phop(j) = -2.d0*(tmpr*cre(i)+tmpi*cim(i))
!     &               /(cre(i)**2+cim(i)**2)
!                endif
!              endif
!c end
!            enddo
!          enddo
!        else
!          write(6,*)"REPFLAG = ",repflag," in GETGRAD"
!          stop
!        endif
!      endif

   99 continue

c save PEM as the diagonal elements of whichever represenation we are using
      do i=1,nsurft
        if (repflag.eq.0) then
        pem(i) = pema(i)
        do j=1,3
        do k=1,nclu
          gpem(j,k,i) = gpema(j,k,i)
        enddo
        enddo
        else
        pem(i) = pemd(i,i)
        do j=1,3
        do k=1,nclu
          gpem(j,k,i) = gpemd(j,k,i,i)
        enddo
        enddo
        endif
      enddo

c overwrite gradients if we are doing a MSX search
      if (tflag(1).eq.4) then

c     get upper state index
      nupper=2
      if (pem(1).gt.pem(2)) nupper=1
      egap=pem(1)-pem(2)

c     get x1 and gupper

      x1mag = 0.d0
      do i=1,nclu
      do j=1,3
      x1(j,i)   = gpem(j,i,1)-gpem(j,i,2)    ! grad(E1-E2)
      x1mag     = x1mag+x1(j,i)**2
      gupp(j,i) = gpem(j,i,nupper)           ! grad of upper surface
      enddo
      enddo
      x1mag = dsqrt(x1mag)                   ! magnitude of grad(E1-E2)

c     for zero nonadiabatic coupling...
c     compute dot prod
      dot1=0.d0
      do i=1,nclu
      do j=1,3
        dot1=dot1+x1(j,i)*gupp(j,i)/x1mag    ! grad upper . x1/x1mag
      enddo
      enddo

c     project out componenets in the direction of x1
      do i=1,nclu
      do j=1,3
      gf(j,i) = 2.d0*egap*x1(j,i)/x1mag      ! f = 2*(E1-E2)*x1/x1mag
      gg(j,i) = gupp(j,i)-dot1*x1(j,i)/x1mag ! g = grad upper projected onto plane perp to f
      msx_grad(j,i) = gg(j,i)+gf(j,i)            ! total gradient
      enddo
      enddo

c     replace gradient with msx_grad
      do i=1,nclu
      do j=1,3
        gv(j,i)=msx_grad(j,i)
      enddo
      enddo

      endif

      return
      end
