c COMMON BLOCK WITH VARIBLES THAT EVOLVE ALONG EACH SINGLE TRAJECTORY

c Note: When "i" is an index, the array runs over atom groups,
c   and when "j" is an index, the array runs over atoms,
c   and when "k" is an index, the array runs over surfaces,
c   and when "a" is an index, the array runs over the x,y,z Cartesian
c   components.

c **********************************************************************
c INITIAL VALUES
c **********************************************************************
c XXI(a,j) : Initial values of XX for current trajectory, assigned in DRIVER
c PPI(a,j) : Initial values of PP for current trajectory, assigned in DRIVER
c CREI(k)  : Initial value of the real part of the electronic coefficients, assigned in INITELEC
c CIMI(k)  : Initial value of the imaginary part of the electronic coefficients
c BIGJI(a) : Initial values of components of total angular momentum
c BIGJTOTI : Initial values of total angular momentum
c EROTI(a) : Initial values of components of the total rotational energy
c EROTTOTI : Initial values of total rotational energy
c PEI      : Initial values of potential energy
c PEAGI(i) : Initial values of potential energy for AG i
c KEI      : Initial values of kinetic energy
c TEI      : Initial values of total energy
c TEMPI    : Initial temperature
c BIGJTOTAGI(i): Initial rotational energy, assigned in INITMOL
c BIGJAGI(a,i) : Principal components of the initial rotational energy, assigned in INITMOL
c EROTTOTAGI(i): Initial total ang mom, assigned in INITMOL
c EROTAGI(a,i) : Principal components of the initial ang mom, assigned in INITMOL
c EINTAGI(i)   : Initial internal energy (= the kinetic energy of the isolated AG), assigned in INITMOL
c EVIBAGI(i)   : Initial vibrational energy (= EINTAGI - EROTTOTAGI), assigned in INITMOL
c BQCI     : Initial value of the impact parameter (used for for IORIENT = 1)
c LQCI     : Initial value of the orbital angular momentum (used for IORIENT = 1)
c EORBQCI  : Initial value of the orbital energy (used for IORIENT = 1)
c ERELQCI  : Initial value of the relative energy (used for IORIENT = 1)
c LSAMPT   : If TRUE -> thermal sampling of J
c SAMPJTEMP: Temperature for sampling J for htis trajectory
c **********************************************************************


c **********************************************************************
c NUCLEAR AND ELECTRONIC VARIABLES
c **********************************************************************
c XX(a,j)   : Unscaled cartesian coordinates.
c PP(a,j)   : Nuclear momenta.
c GV(a,j)   : Nuclear gradient.
c CRE(2*k)  : Real part of the electronic coefficients
c CIM(2*k)  : Imaginary part of the electronic coefficients
c GCRE(2*k) : Time derivative of real part of electronic coefficients
c GCIM(2*k) : Time derivative of imaginary part of electronic coefficients
c Note      : For the CSDM method, CRE, CIM, GCRE, and GCIM also contain
c             the coherent parts of these quantities stored as 
c             elements NSURFT+1 to 2*NSURFT.
c DVEC(a,j,k,k) : Components of the nonadiabatic coupling vector for
c                 every set of potential energy surfaces.
c PHASE(k)  : Integration of the potnential energies.
c PEM(k)    : Potential energies for whichever representation we are using.
c **********************************************************************


c **********************************************************************
c OTHER VARIABLES
c **********************************************************************
c HSTEP   : Stepsize
c BIGJ(a) : Components of total angular momentum
c BIGJTOT : Total angular momentum
c EROT(a) : Components of rotational energy
c EROTTOT : Total rotational energy
c PE      : Current potential energy
c KE      : Current kinetic energy
c TE      : Current total energy
c TEMP    : Current temperature
c ATEMP   : <temperature>*time over entire trajectory (i.e., ATEMP/TIME = <temperature>)
c ATEMP2  : <temperature^2>*time over entire trajectory
c STEMP   : Standard deviation for the temperature over the entire trajectory
c TIME    : Current time
c NSURF   : Current electronic surface
c ISTEP   : Number of integrator steps
c ISTEPW  : Number of integrator steps with geom in PE well LWELL
c ITRAJ   : Trajectory number
c NISTEP  : Total number of integrator steps for all trajectories up to current step
c NISTEPW : Same for geom in PE well LWELL
c OUTCOME : Outcome index for TERMFLAG = 3
c AIND    : Used for TERMFLAG =3, atom labels of broken bond
c LDRIVER : if 1 reject traj due to pei>sample0im for vibr state selected coords
c WELLI   : tells in which well the current geom is (1-acetylene,2-vinyledene)
c **********************************************************************

      double precision xx(3,mnat),pp(3,mnat),crei(mnsurf),cimi(mnsurf),
     &    gv(3,mnat),bqci,lqci,eorbqci,erelqci,peagi(mnmol),
     &    xxi(3,mnat),ppi(3,mnat),bigji(3),bigjtoti,pei,kei,tei,tempi,
     &    cre(2*mnsurf),cim(2*mnsurf),hstep,pe,ke,te,bigj(3),
     &    bigjtot,temp,erot(3),erottot,erottotagi(mnmol),
     &    erotagi(3,mnmol),bigjtotagi(mnmol),bigjagi(3,mnmol),
     &    evibagi(mnmol),eintagi(mnmol),erottoti,eroti(3),
     &    gcre(2*mnsurf),gcim(2*mnsurf),time,phase(mnsurf),pem(mnsurf),
     &    atemp,stemp,atemp2,sampjtemp(mnmol),rtrans
      integer nsurf,istep,itraj,outcome,aind(mnoutcome,2)
      integer ldriver,welli,istepw,nistep,nistepw

      common /c_traj/ xx,pp,crei,cimi,gv,cre,cim,xxi,ppi,hstep,phase,
     &   bigji,bigjtoti,pei,kei,tei,tempi,bigj,bigjtot,erot,erottot,
     &   erottotagi,erotagi,bigjtotagi,bigjagi,evibagi,eintagi,peagi,
     &   pe,ke,te,pem,temp,atemp,stemp,atemp2,sampjtemp,rtrans,
     &   gcre,gcim,time,bqci,lqci,
     &   eorbqci,erelqci,nsurf,istep,itraj,outcome,aind,
     &   istepw,nistep,nistepw
