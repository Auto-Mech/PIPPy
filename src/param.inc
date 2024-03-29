      integer maxatom,maxpair,maxdata,maxterm,maxperm,maxcoef,maxchan
      parameter(maxatom=16)	        ! max atoms+1
      parameter(maxchan=2)	        ! max fragment channels
      parameter(maxpair=maxatom*(maxatom-1)/2)	! max atom-atoms pairs
      parameter(maxdata=2000000)	        ! max fitting data
      parameter(maxterm=6000000)	! max number of polynomial terms
      parameter(maxcoef=4000)	        ! max number of coefficients
      parameter(maxperm=70000)		! max number of symmetry permutations

