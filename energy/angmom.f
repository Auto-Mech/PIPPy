      subroutine angmom(xx,pp,mm,natom,bigj,bigjtot)

c computes angular momentum around the origin using the position XX and momentum PP
c components of angmom are the cartesian components, not the components around the the principal axes

      implicit none
      include 'param.f'

      integer natom
      double precision xx(3,natom),pp(3,natom),mm(natom),
     &   bigj(3),bigjtot

c local
      integer i,k,l

c zero
      bigjtot=0.d0
      do i=1,3
       bigj(i)=0.d0
      enddo

c bigj = xx cross pp (J = R x P)
      do i=1,natom
!        print *,"xx",(xx(k,i),k=1,3),"pp",(pp(l,i),l=1,3)
        bigj(1)=bigj(1)+(xx(2,i)*pp(3,i)-xx(3,i)*pp(2,i))
        bigj(2)=bigj(2)-(xx(1,i)*pp(3,i)-xx(3,i)*pp(1,i))
        bigj(3)=bigj(3)+(xx(1,i)*pp(2,i)-xx(2,i)*pp(1,i))
      enddo

c bigj(i) is the ith component of the total angular momentum
c bigjtot is the magnitide of the total angular momentum
      do i=1,3
       bigjtot = bigjtot + bigj(i)**2
      enddo
      bigjtot = dsqrt(bigjtot)

      end
