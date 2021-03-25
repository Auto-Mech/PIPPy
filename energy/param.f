c This file stores constants used for dimensioning arrays
c and for unit conversions.

c     max number of atom groups (aka molecules)
      integer mnmol
      parameter(mnmol=2)

c     max number of total atoms (total for all molecules)
      integer mnat
      parameter(mnat=100)

c     max number of electronic states
      integer mnsurf
      parameter(mnsurf=12)

c     max number of trajectories
      integer mntraj
      parameter(mntraj=10000000)

c     max number of termination condition outcomes
      integer mnoutcome
      parameter(mnoutcome=5)

c     maximum dimension of the Y-array that is integrated
      integer mnyarray
      parameter(mnyarray=6*mnat+5*mnsurf)

c     unit conversions
      double precision amutoau,kb,autoang,autoev,autofs,mu,autocmi,pi,
     & autokcal
      parameter(amutoau=1822.844987d0) ! best
c      parameter(amutoau=1822.888506d0) ! used by NAT
      parameter(kb=3.166829d-6)  ! Boltzmann in hartee/K
      parameter(autoang=0.52917706d0)
      parameter(autoev=27.2113961d0)
      parameter(autofs=0.024189d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)   ! IMLS
      parameter(mu=1.d0*amutoau)  ! mass-scaling mass
      parameter(pi=3.1415926536d0)
