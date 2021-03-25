c COMMON BLOCKS FOR FREQUENTLY USED VARIABLES THAT ARE THE SAME FOR ALL
c TRAJECTORIES

c Note: When "i" is an index, the array runs over atom groups,
c   and when "j" is an index, the array runs over atoms,
c   and when "k" is an index, the array runs over surfaces.


c **********************************************************************
c CONTROL VARIABLES
c **********************************************************************
c potflag    Flag that controls which potential interface to use
c nprint     Print info to unit 6 every nprint steps
c ntraj      Number of trajectories
c methflag   Flag for which single surface or non-BO method is used
c tflag      Special trajectory options
c repflag    Representation flag, =0 for adiabatic, =1 for diabatic
c ldofrag    If true, compute fragment energies. This option is useful 
c            for direct dynamics calculations where the fragments 
c            require different input options than the supermolecule.
c            Some options cannot be performed with this turned off.
c stodecoflag If true, use the stochastic decoherence.

      integer potflag,nprint,ntraj,methflag,tflag(4),repflag
      logical ldofrag,stodecoflag
      common/c_control/potflag,nprint,ntraj,methflag,tflag,repflag,
     &        ldofrag,stodecoflag
c **********************************************************************


c **********************************************************************
c SYSTEM
c **********************************************************************
c mm(j)       mass of atom j (input in amu, converted immediately to au)
c mmag(i)     total mass of atom group i (au, computed in readin)
c nsurft      number of electronic surfaces
c nmol        # of atom groups
c nat         # of total atoms (for all atom groups)
c natom(i)    # of atoms in AG i
c iatom(i)    index of first atom in AG i
c symbol(j)   atomic symbol of atom j
c ezero       zero of energy for the full potential (all atom groups), 
c             this number gets subtracted from the energies of every 
c             surface for all potential calls in GETPEM (input in eV, 
c             converted to au)
c ezeroim(i)  zero of energy for each AG, this value is used to correct
c             the zero of energy for each atom group (input in eV,
c             converted to au)
c emin2(i)    energy of the classical minimum for each well (input in eV)

      double precision ezero
      double precision mm(mnat),mmag(mnmol),emin2(mnmol),ezeroim(mnmol)
      integer nsurft,nmol,nat,natom(mnmol),iatom(mnmol)
      character*2 symbol(mnat)
      common/c_sys/ezero,ezeroim,mm,mmag,emin2,nsurft,nmol,nat,
     &             natom,iatom,symbol
c **********************************************************************


c **********************************************************************
c INITIAL CONDITIONS
c **********************************************************************
c nsurf0     Initial electronic surf for all trajectories
c initx(i)   Initial conditions flag for AG i (coords)
c initp(i)   Initial conditions flag for AG i (momenta)
c initj(i)   Initial conditions flag for AG i (ang mom)
c xx0(3,j)   Initial coordinates, read by readin (input in A, converted 
c            immediately to bohr)
c pp0(3,j)   Initial momenta, read by readin (input in au)
c xx02(3,j)   Initial coordinates for min of PE weel 2, optional, for INITx=8 
c rran(3)    Numbers used to generate randum clusters if INITx = 1 (input in 
c            angstroms, converted immediately to bohr)
c temp0im(i) Target temperature of AG i when INITp = 0 (Kelvin)
c lescale0(i) If TRUE -> scale momenta to a target kinetic energy
c escale0im(i) Target kinetic energy of AG i when INITp = 0 and lescale0=TRUE 
c            (input in eV, immediately converted to hartree)
c scale0im(i) Target total energy of AG i when INITx = 5 and scale0im > 0

c The following are used for atom-diatom scattering conditions (INITx = 3)
c vvad        Initial vibrational quantum number of the diatom
c jjad        Initial rotational quantum number of the diatom
c arrad      Initial molecular arrangement
c escatad    Fixed total energy (input in eV, converted immediately to au)
c rrad        Initial atom-diatom distance (input in A, converted 
c            immediately to bohr)
c rinad      Inner turning point of initial diatom
c routad     Outer turning point of initial diatom
c tauad      Period in au of initial diatom
c bmaxad     Maximum impact parameter in au
c bminad     Minimum impact parameter in au
c pprelad    Initial relative momenta in au
c ecolad     Initial energy in au
c easym(k,m) Minimum energy for the diatom in arragement m
c rasym(k,m) Minimum energy distance for the diatom in arragement m

c The following are used for the methods based on normal modes
c nmqn(n,i)    Initial normal mode quanutm number for mode n
c nmtype(i)    Initial geometry is a well (0) or a saddle point (1).
c nmvec(3,j,m,i) Normal mode vectors for the m^th normal mode
c freq(m,i)    Normal mode frequencies for the m^th normal mode
c rturn(m,i)   Normal mode turning points for the m^th normal mode
c ewell(i)    Energy of the initial structure when using normal mode method for initial
c            conditions.  Should correspond to the energy of the bottom of some well.
c lreadhess  T=read Hessian from fort.70. F=compute Hessian numerically from gradients
c nmvec2(3,j,m,i)
c freq2(m,i)
c rturn2(m,i)
c ewell2(i)  similar variables for EMS sampling from 2 PE wells    

c The following are used for diatom initial conditions (INITx = 4)
c vvdi(i)        Initial vibrational quantum number of the diatom for AG i
c jjdi(i)        Initial rotational quantum number of the diatom for AG i
c rmindi(i)      Guess at min diatomic separation for AG i, will be refined by code

c The following are used for reading initial conditions from a separate file (INITx=6)
c samptot(i)     Total number of samples contained in the separate file
c lbinsamp(i)    Samples files binary (direct access) or not
c lbinsamp2(i)   Samples files binary or not for second PE well
c sampfile(i)    Filename containing sampled data
c sampfile2(i)   Filename containing sampled data for the second PE well
c sampwell(im)   Number of wells to sample from
c relwell(im)    Weights for sampling from each well
c lems(im)       If true use EMS
c emicr(im)      Total energy of the microcanonical ensemble 
c nemstot(im)    Total number of points to output using EMS (efficient microcan. sampl.)
c ninc(im)       Incubation number of steps for EMS
c nbrea(im)      Number of steps in between points(trajectories) sampled in EMS
c emsw(im)       last calculated weight in EMS 

c The following are used for photoexcited trajectories (TFLAG(3) = 1)
c ntarget     The target electronic state
c ephoton     The photon energy in au
c wphoton     The photon width / 2 in au (i.e., trajectories are excited 
c             if the electronic energy gap is ephoton +/- wphoton)

      integer nsurf0,initx(mnmol),initp(mnmol),initj(mnmol)
      double precision xx0(3,mnat),pp0(3,mnat),rran(3),xx02(3,mnat),
     &    temp0im(mnmol),escale0im(mnmol),scale0im(mnmol),emicr(mnmol)

      integer vvad,jjad,arrad
      double precision escatad,rrad,rinad,routad,
     & tauad,bmaxad,bminad,pprelad,ecolad,
     & easym(mnsurf,3),rasym(mnsurf,3)

      integer nmtype(mnmol)
      double precision nmvec(3,mnat,3*mnat,mnmol),freq(3*mnat,mnmol),
     & rturn(3*mnat,mnmol),ewell(mnmol),nmqn(3*mnat,mnmol),
     & nmvec2(3,mnat,3*mnat,mnmol),freq2(3*mnat,mnmol),
     & rturn2(3*mnat,mnmol),ewell2(mnmol)
      logical lreadhess,lbinsamp(mnmol),lbinsamp2(mnmol),lems(mnmol),
     & letot(mnmol)

      double precision rmindi(mnmol),vvdi(mnmol),jjdi(mnmol),
     &    rindi(mnmol),routdi(mnmol),taudi(mnmol),xtaudi(mnmol)

      integer samptot(mnmol),samptot2(mnmol),nemstot(mnmol),
     & sampwell(mnmol),ninc(mnmol),nbrea(mnmol)
      character*10 sampfilexx(mnmol),sampfilepp(mnmol)
      character*10 sampfilexx2(mnmol),sampfilepp2(mnmol)
      double precision samptarg(mnmol),sampjmin(mnmol),sampjmax(mnmol),
     & sampjtemp1(mnmol),sampjtemp2(mnmol),
     & sampjbrot1(mnmol),sampjbrot2(mnmol),
     & relwell1(mnmol),relwell2(mnmol),ejsc(4,mnmol),emsw(mnmol)

      integer ntarget
      double precision ephoton,wphoton

      common/c_initial/xx0,pp0,rran,xx02,
     & temp0im,escale0im,scale0im,emicr,escatad,rrad,rinad,routad,
     & tauad,bmaxad,bminad,pprelad,ecolad,easym,rasym,
     & nmvec,nmqn,freq,rturn,ewell,
     & nmvec2,freq2,rturn2,ewell2,
     & rmindi,vvdi,jjdi,rindi,routdi,taudi,xtaudi,
     & ephoton,wphoton,
     & samptarg,sampjmin,sampjmax,sampjtemp1,sampjtemp2,emsw,
     & sampjbrot1,sampjbrot2,ejsc,relwell1,relwell2,
     & sampfilexx,sampfilepp,
     & sampfilexx2,sampfilepp2,
     & ntarget,
     & nsurf0,initx,initp,initj,
     & vvad,jjad,arrad,
     & nmtype,
     & samptot,samptot2,nemstot,sampwell,ninc,nbrea,
     & lreadhess,lbinsamp,lbinsamp2,lems,letot
c **********************************************************************


c **********************************************************************
c INTEGRATOR VARIABLES
c **********************************************************************
c hstep0 = initial stepsize (input in fs, converted immediately to au)
c eps = converge energy for each step to this tolerance when using
c       the BS integrator (au)
c intflag = flag that controls which integrator to use

      double precision hstep0,eps
      integer intflag
      common/c_integrator/hstep0,eps,
     & intflag
c **********************************************************************

c **********************************************************************
c ORIENTATIONAL VARIABLES
c **********************************************************************
c iorient = determines how to orient AGs
c bmaxqc = maximium impact parameter
c bminqc = minimium impact parameter
c erelqc = relative energy
c tempqc = temperature of the relative energy distribution
c rel0qc = initial separation
c The following pertain to orienting the AGs with respect to each other
c comxx(i)   Initial CoM coords for AG i (input in A, converted
c            immediately to bohr)
c compp(i)   Initial CoM momenta for AG i (input in au)
      double precision
     & comxx(3,mnmol),compp(3,mnmol)
      integer iorient
      double precision bmaxqc,erelqc,rel0qc,bminqc,tempqc
      common/c_orient/comxx,compp,bminqc,bmaxqc,erelqc,tempqc,
     & rel0qc,iorient

c **********************************************************************
c TERMINATION CONDITIONS
c **********************************************************************
c t_stime = For TERMFLAG =1, Terminate after T_STIME fs.
c t_gradmag = For TERMFLAG =2, Terminate when the magnitude of the
c             gradient is less than T_GRADMAG in eV/A.
c t_r(o) = For TERMFLAG = 3, Termination distance, input in A, 
c          immediately converted to bohr.  
c t_nstep = Maximum number of steps for each trajectory.
c lwell = 0 geoms in one PE well are not selected, 1-acetylene, 2-vinyledene

      double precision t_stime,t_gradmag,t_r(mnoutcome)
      integer termflag,t_nstep,t_noutcome,lwell
      character*2 t_symb(mnoutcome,2)
      common/c_term/
     & t_r,t_stime,t_gradmag,termflag,t_nstep,t_noutcome,t_symb,lwell
c **********************************************************************


c **********************************************************************
c SPECIAL TRAJECTORY OPTIONS
c **********************************************************************
c ramptime = ramp trajectory momentum every RAMPTIME time units 
c            (for TFLAG(1) = 2).  Read in fs, converted to au immediately.
c rampfact = ramp trajectory momentum by RAMPFACT (for TFLAG(1) = 2).
c nramp = ramp trajectory momentum NRAMP times (for TFLAG(1) = 2).
c andersen_temp = used by the Andersen thermostat, temperature
c andersen_freq = used by the Andersen thermostat, collision frequency
c scandth = Total energy to scale momenta after Andersen resampling (if >0)
      double precision ramptime,rampfact,andersen_temp,andersen_freq,
     & scandth
      integer nramp
      common/c_special/ramptime,rampfact,
     & andersen_temp,andersen_freq,scandth,
     & nramp
c **********************************************************************

c **********************************************************************
c MONTE CARLO INTEGRATION VARIABLES
c **********************************************************************
c mccurv = number of normal modes to treat as something other than a 
c          cartesian displacement (e.g., as a torsion or bend)
c mcmode = list of modes to be treated as a torsion, bend, etc.
c mctype = type of mode (1 = torsion)
c mcpar = additional parameters
      integer mccurv,mcmode(3*mnat),mctype(3*mnat),mcpar(2,3*mnat)
      common/c_montecarlo/mccurv,mcmode,mctype,mcpar
c **********************************************************************
