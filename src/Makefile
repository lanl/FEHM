#***********************************************************************
#
# Makefile for FEHM
#
#***********************************************************************
# Usage: make [VARIABLE] TYPE
#
# Types:
#    help                      view the help screen
#    all                       build release FEHM into ./src/     
#    debug                     build debug FEHM into ./src/   
#    install                   build release FEHM into PREFIX
#    uninstall                 remove release FEHM from PREFIX        
#    clean                     remove all object and mod files
#    test                      run FEHM unit tests
#
# Variables:
#    FC                        FORTRAN 77/90 compiler
#    EXE                       Compiled FEHM executable
#    SRCDIR                    Path to FEHM source (default: cwd)
#    WITHSILO 				   Build with SILO HDF5 TPL (DEFUALT)
#***********************************************************************
# GNU Make version 3.79.1 or later should be used
# make should be executed from the directory where the objects 
# will reside.
# Edit SRCDIR to point to the directory where the source code is located
# Edit FC,FFLAGS,EXE for your environment and compiler 
#***********************************************************************

# Compiletime variables
FC = gfortran
EXE = xfehm_v3.3.1
SRCDIR = ../src/
PREFIX = /usr/local/bin
WITHSILO := 1

# Build system specific variables
DEPEND = ${SRCDIR}Makefile.depends
OPSYS = $(shell uname -s )
DATETAG = $(shell date '+%y-%m-%d')
DATE = $(shell date '+.%d%b%y')
FFLAGS = -frecord-marker=4 -ffixed-line-length-132
DEBUGFLAGS = -g -O0 -frecord-marker=4 -fbounds-check -Wall


# Define the written OS tag
ifeq (${OPSYS}, Linux)
	OSTAG = lbUbuntu16
else ifeq (${OPSYS}, Darwin)
	OSTAG = macOS
else ifeq (${OPSYS}, Windows)
	OSTAG = Windows
else ifeq (${OPSYS}, Sun)
	OSTAG = Sun
else
	OSTAG = OPSYS
endif


# Define special permissions for gfortran
# Namely, accessing two Fortran files from an
# alternate directory due to compatibility issues
ifeq (${FC}, gfortran)
	ALTDIR = ${SRCDIR}GFORTRAN/
else
	ALTDIR = ${SRCDIR}
endif

SILOFLAGS := -cpp
ifeq (${WITHSILO}, 0)
	SILOFLAGS += -DNOSILO
else
	FFLAGS += -L/scratch/ymp/smckinney/silo/silo_install/lib -L/scratch/ymp/smckinney/hdf5/hdf5-1.8.11-linux-x86_64-shared/lib/ -I/scratch/ymp/smckinney/silo/silo_install/include -L/scratch/ymp/smckinney/silo/silo_install/lib -L/scratch/ymp/smckinney/hdf5/hdf5-1.8.11-linux-x86_64-shared/lib/ -lhdf5 -lstdc++
endif

# Define the helpscreen that will show on `make help`
define help
#====================================================================================
#
#   FEHM V3.3.1
#   (c) 2018 LOS ALAMOS NATIONAL LABORATORY
#
#====================================================================================

    make .or. make all          build release FEHM into ./src/     
    make debug                  build debug FEHM into ./src/   
    make install                build release FEHM into PREFIX (default: /usr/local/bin)  
    make uninstall              remove release FEHM from PREFIX         
    make clean                  remove all object and mod files
    make test                   run FEHM unit tests
    make help                   view this help screen

MAKEFILE CONFIGURATION

    EXE:                        ${EXE}
    SRCDIR:                     ${SRCDIR}
    OS:                         ${OSTAG}
    FC:                         ${FC}
    FFLAGS:                     ${FFLAGS}
    DEBUGFLAGS:                 ${DEBUGFLAGS}
    DATETAG:                    ${DATETAG}
    DATE:                       ${DATE}
    PREFIX:                     ${PREFIX}
endef
export help

# Define FORTRAN and FORTRAN90 objects
OBJECTS     := $(patsubst ${SRCDIR}%.f,%.o,$(wildcard ${SRCDIR}*.f))
OBJECTS_F90 := $(patsubst ${SRCDIR}%.f90,%.o,$(wildcard ${SRCDIR}*.f90))

# Define targets
.PHONY: all
all: FFLAGS += -O2
all: dated xfehm

.PHONY: install
install: FFLAGS += -O2
install: dated xfehm
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp ${EXE} $(DESTDIR)$(PREFIX)/bin/fehm
	@echo "Installed to $(DESTDIR)$(PREFIX)/fehm"

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/fehm      
	@echo "$(DESTDIR)$(PREFIX)/fehm removed."

.PHONY: clean
clean:
	rm *.o
	rm *.mod

.PHONY: help
help:
	@echo "$$help"

.PHONY: test
test:
	@echo "Running tests..."
	$(eval FEHM_CWD:=$(shell pwd))
	cd ${SRCDIR}/../fehmpytests/; python fehmpytests.py ${FEHM_CWD}/${EXE}
	cd ${FEHM_CWD}

.PHONY: dated
dated:
	rm -f ${SRCDIR}dated.f
	rm -f dated.o
	sed s/OS\ DATE/${OSTAG}\ ${DATETAG}/ ${SRCDIR}dated.template > ${SRCDIR}dated.f

xfehm: ${OBJECTS} ${OBJECTS_F90}
	${FC} ${FFLAGS} ${OBJECTS} ${OBJECTS_F90} -o ${EXE} -lsiloh5 -lhdf5 -lstdc++

debug: FFLAGS += -g	-fbounds-check
debug: dated ${OBJECTS} ${OBJECTS_F90}
	${FC} ${DEBUGFLAGS} ${OBJECTS} ${OBJECTS_F90} -o ${EXE}


# Define rules
%.o : ${SRCDIR}%.f
	${FC} ${FFLAGS} ${SILOFLAGS} $< -c

%.o : ${SRCDIR}%.f90
	${FC} ${FFLAGS} $< -c

dated.o : ${SRCDIR}dated.f
	${FC} ${FFLAGS} $< -c


# Special handling of two files
inrestart.o : ${ALTDIR}inrestart.f
	${FC} ${FFLAGS} $< -c
insptr.o : ${ALTDIR}insptr.f
	${FC} ${FFLAGS} $< -c


# Include `Makefile.depends`
include ${DEPEND}
