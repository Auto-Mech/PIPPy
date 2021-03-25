#Makefile for dint

#set POT here or use 'gmake POT=potname'
#POT=exp6
TINKER=0
POT=h3o
#POT=tbplusexp6all
#POT=tbplusexp6-2
#POT=ch4hALL
#POT=ch4ohALL
#POT=ch4oALL
#POT=h2oh2o

#TINKER PES
# atom types: ~/tinker/tinker/params/mm3.prm 
# NOTE: obj_tinker contains a link to the mm3 PES
# NOTE: obj_tinker2 contains a link to the amoeba09 PES
# NOTE: obj_tinker3 contains a link to the oplsaa PES
#TINKER=1
#POT=rohMP
#POT=roohT
#POT=alkT
#POT=cih2o2m
#POT=cih2o2mPIP
#TINKER PES

EXEDIR = ../exe
TARG = $(EXEDIR)/pese-$(POT).x$(EXESUFFIX)
OBJDIR = obj
OBJDIR_OPT = obj_opt
SPRNGLIB = /home/mobergdr/codes/sprng

#  FC = f77
  FC = gfortran
  LD = gfortran
  OPT_FFLAGS = -O3
  FTPPFLAGS = -I$(SPRNGLIB)/SRC -I.. -c
ifeq ($(TINKER),0)
  FFLAGS = -c
else
  TINKLIB = /home/jasper/KTP/dint/tinker/tinker/source
  FFLAGS = -c -I$(TINKLIB)
endif
  COMPLIB = -llapack -lblas
  RNGLIB = -lmlfg
  LIBS = $(COMPLIB) $(RNGLIB)
#  LDFLAGS = -O3 -g -L../../sprng/lib
#  LDFLAGS = -O3 -g -L../../sprng/lib -L/soft/apps/packages/lapack-3.0/lib
  LDFLAGS = -O3 -g -L$(SPRNGLIB)/lib -L/usr/lib64 -L/home/jasper/lib

VPATH = ..:../../autopot

SRC =	$(POT).f \
        angmom.f \
        ange.f \
	driver.f \
        energy.F \
	getgrad.f \
	getpem.f \
        gettemp.f \
        header.f \
        initmol.f \
        premol.f \
        readin.f

OBJ = $(patsubst %.f,%.o,$(SRC))
OBJ := $(patsubst %.F,%.o,$(OBJ))

EXEDIR:=../$(EXEDIR)
MAKEOPT = -f ../Makefile
all: yes_opt

no_opt:
	$(MAKE) $(MAKEOPT) -C obj targ
yes_opt:
	$(MAKE) $(MAKEOPT) -C obj_opt EXESUFFIX=.opt OPT="$(OPT_FFLAGS)" targ

targ: $(TARG)

# dependencies
angmom.o:     param.f
driver.o:     param.f c_sys.f c_traj.f c_ran.f c_output.f
energy.o:     param.f c_sys.f c_traj.f c_ran.f
getgrad.o:    param.f c_sys.f
getpem.o:     param.f c_sys.f
gettemp.o:    param.f
header.o:
initmol.o:    param.f c_sys.f c_traj.f
premol.o:     param.f c_sys.f
readin.o:     param.f c_sys.f c_traj.f c_ran.f

COMPILE = $(FC) $(FFLAGS) $(OPT)

$(TARG): $(OBJ)
ifeq ($(TINKER),0)
	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) -o $(TARG)
else
#	$(LD) $(LDFLAGS) -ffast-math -fopenmp -static-libgcc $(OBJ) $(LIBS) -o $(TARG) ../obj_tinker/libtinker.a ../obj_tinker/libfftw3_threads.a ../obj_tinker/libfftw3.a
	$(LD) $(LDFLAGS) -ffast-math -fopenmp -static-libgcc $(OBJ) $(LIBS) -o $(TARG) $(TINKLIB)/libtinker.a $(TINKLIB)/libfftw3_threads.a $(TINKLIB)/libfftw3.a
endif

%.o : %.f
	$(COMPILE) $< -o $@

%.o : %.F
	$(COMPILE) $(FTPPFLAGS) $< -o $@

clean:
	-rm -f obj/* obj_opt/*

realclean:
	-rm -f obj/* obj_opt/* ../exe/*

