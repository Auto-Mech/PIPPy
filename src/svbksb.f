      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION b(mp),u(mp,np),v(np,np),w(np),x(np)
      include 'param.inc'
c      PARAMETER (NMAX=2000)    ! Maximum anticipated value of n.
c      Solves A â X = B for a vector X, where A is specified by the arrays u, w, v as returned by
c      svdcmp. m and n are the logical dimensions of a, and will be equal for square matrices. mp
c      and np are the physical dimensions of a. b(1:m) is the input right-hand side. x(1:n) is
c      the output solution vector. No input quantities are destroyed, so the routine may be called
c      sequentially with different bÃ­s.
      INTEGER i,j,jj
      DOUBLE PRECISION s,tmp(maxdata)
      do 12 j=1,n             ! Calculate UTB.
      s=0.
      if(w(j).ne.0.)then      ! Nonzero result only if wj is nonzero.
      do 11 i=1,m
      s=s+u(i,j)*b(i)
 11   enddo
      s=s/w(j)                ! This is the divide by wj .
      endif
      tmp(j)=s
 12   enddo 
      do 14 j=1,n             ! Matrix multiply by V to get answer.
      s=0.
      do 13 jj=1,n
      s=s+v(j,jj)*tmp(jj)
 13   enddo 
      x(j)=s
 14   enddo 
      return
      END

