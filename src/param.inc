      double precision autoang,autocmi
      integer maxatom,maxpair,maxdata,maxterm,maxperm,maxcoef,maxchan
      parameter(autoang=0.529177249d0)
      parameter(autocmi=219474.63067d0)
      parameter(maxatom=18)	! max total atoms
      parameter(maxchan=4)	! max fragment channels
      parameter(maxpair=maxatom*(maxatom-1)/2)	! max atom-atoms pairs
      parameter(maxdata=150000)	! max fitting data
      parameter(maxcoef=5000)	! max number of coefficients
      parameter(maxterm=100000000)	! max number of polynomial terms
      parameter(maxperm=1000000)	! max number of symmetry permutations

