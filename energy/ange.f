      subroutine ange(xx,pp,mm,natom,eig,bigj,bigjtot,erot,erottot)

c Compute rotational energies around principal axes

      implicit none
      include 'param.f'

      integer natom
      double precision xx(3,mnat),pp(3,mnat),mm(mnat),
     &   bigj(3),bigjtot,xxp(3,mnat),ppp(3,mnat)

c local
      integer i,j,info,ismall
      double precision erot(3),erottot,
     & temp1,temp2,
     & mom(3,3),ap(6),eig(3),rot(3,3),work(9),temp3,
     & mom2,rr

c single atom
      if (natom.eq.1) then
         bigj(1) = 0.d0
         bigj(2) = 0.d0
         bigj(3) = 0.d0
         bigjtot = 0.d0
         erot(1) = 0.d0
         erot(2) = 0.d0
         erot(3) = 0.d0
         erottot = 0.d0
         return
      endif
      
c I modified the polyatomic part to work for linears instead of this
c diatom
c      if (natom.eq.2) then
c      rr=0.d0
c      do j=1,3
c      rr=rr+(xx(j,1)-xx(j,2))**2
c      enddo
c      mom2=mm(1)*mm(2)/(mm(1)+mm(2))*rr   ! remember rr=R**2
c      print *,"mom2=",mom2
c      return
c      endif

c natom >= 2
c compute moment of intertia matrix mom
      do i=1,3
      do j=1,3
          mom(i,j) = 0.d0
      enddo
      enddo

      do i=1,natom
        mom(1,1)=mom(1,1)+mm(i)*(xx(2,i)**2+xx(3,i)**2)
        mom(2,2)=mom(2,2)+mm(i)*(xx(1,i)**2+xx(3,i)**2)
        mom(3,3)=mom(3,3)+mm(i)*(xx(1,i)**2+xx(2,i)**2)
        mom(1,2)=mom(1,2)-mm(i)*(xx(1,i)*xx(2,i))
        mom(1,3)=mom(1,3)-mm(i)*(xx(1,i)*xx(3,i))
        mom(2,3)=mom(2,3)-mm(i)*(xx(2,i)*xx(3,i))
      enddo
      mom(2,1)=mom(1,2)
      mom(3,1)=mom(1,3)
      mom(3,2)=mom(2,3)

c      print *,"mom"
c      print *,(mom(1,j),j=1,3)
c      print *,(mom(2,j),j=1,3)
c      print *,(mom(3,j),j=1,3)

c      call angmom(xx,pp,mm,natom,bigj,bigjtot)
c      print *,"j",bigj,bigjtot

c     diagonalize the mom matrix
      do i=1,3
      do j=i,3
        ap(i+(j-1)*j/2)=mom(i,j)
      enddo
      enddo
      call dspev( 'v','u',3,ap,eig,rot,3,work,info )

c      print *,"rot"
c      print *,(rot(1,j),j=1,3)
c      print *,(rot(2,j),j=1,3)
c      print *,(rot(3,j),j=1,3)
c      print *,"eig"
c      print *,(eig(j),j=1,3)

c rotate to diagonalize mom
      do i=1,natom
         temp1 = xx(1,i)
         temp2 = xx(2,i)
         temp3 = xx(3,i)
         xxp(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         xxp(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         xxp(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
         temp1 = pp(1,i)
         temp2 = pp(2,i)
         temp3 = pp(3,i)
         ppp(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         ppp(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         ppp(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
      enddo
      
      call angmom(xxp,ppp,mm,natom,bigj,bigjtot)

c rotational energies
c      write(6,*)"Ia = ",(eig(j),j=1,3)," au"
      erottot = 0.d0
      ismall=0
      do i=1,3
        if (eig(i).lt.1.d-10) then
          ismall=ismall+1
          erot(i) = 0.d0
        else 
          erot(i) = bigj(i)**2/(2.d0*eig(i))
        endif
        erottot=erottot+erot(i)
      enddo
c      write(6,*)"Erot = ",(erot(j)*autoev,j=1,3)," eV"
      write(6,*)"Ji = ",(bigj(j),j=1,3)
      if (ismall.eq.1) write(6,*)"Found & skipped 1 small Ia,",
     & " must be a linear molecule"
      if (ismall.gt.1) then
        write(6,*)"Found ",ismall," small Ia! Problem."
        write(6,*)"The values are:",(eig(j),j=1,3)," au"
        stop
      endif

c      print *,"ange"
c      print *,bigj(1),erot(1)*autoev,eig(1)
c      print *,bigj(2),erot(2)*autoev,eig(2)
c      print *,bigj(3),erot(3)*autoev,eig(3)

      return
 
      end
