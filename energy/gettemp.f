      subroutine gettemp(pp,mm,natom,temp,ke)

c Compute the temperature and kinetic energy from the momentum.

      implicit none
      include 'param.f'

      integer natom
      double precision pp(3,mnat),mm(mnat),temp,ke

c local
      integer i,j,k

      temp=0.d0
      ke=0.d0
      do i=1,natom
      do j=1,3
         ke = ke + 0.5d0*pp(j,i)*pp(j,i)/mm(i)
      enddo
      enddo
c <KE> = KE/Natom = 3/2 k T
      temp = 2.d0*ke/(kb*3.d0*dble(natom))
      return
      end
