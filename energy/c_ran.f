c COMMON BLOCK FOR FREQUENTLY USED VARIABLES THAT RELATE TO SPRNG

c NOTE: This bit of code has memory problems. 
c       When trajlist is used, SPRNG will crash when called.
c       It seems to be fixed once I reordered the elements in the block.
c       If problems occur, maybe don't call sprng using a variable in 
c       the common block; use a temporary variable instead?

c ranseed = randum number seed
c rng_stream = randum number stream
c trajlist = list of restart trajectory indices
c maxtraj = maximum trajectory index for this run

      integer ranseed,rng_stream,maxtraj,trajlist(mntraj)

      common /c_ran/ranseed,maxtraj,trajlist
      common /c_ran0/rng_stream
