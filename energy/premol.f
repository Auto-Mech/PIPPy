      subroutine premol(im)

c Precompute information for each atom group.

      implicit none
      include 'param.f'
      include 'c_sys.f'
 
      double precision xxm(3,mnat),ppm(3,mnat),mmm(mnat)
      double precision nmvecm(3,mnat,3*mnat),freqm(3*mnat),
     & rturnm(3*mnat),ewellm,nmqnm(3*mnat)
      integer im,i,j,ii,k
      character*2 symb(mnat)

      write(6,*)"Precomputing info for AG ",im
      write(6,*)"------------------------------"

c     store coords for this molecule in their own array for easier manipulation
      do i=1,natom(im)
       ii=i+iatom(im)
       mmm(i) = mm(ii)
       symb(i) = symbol(ii)
       do j=1,3
        xxm(j,i) = 0.d0
        ppm(j,i) = 0.d0
       enddo
      enddo

      return

      end
