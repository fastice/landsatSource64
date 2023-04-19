#
# Modify this section to point to where stuff is.
# Current names are for a specific file system.
# ROOTDIR: root directory for code (e.g. /Users/username directory). Likely should change for linux
# PROGDIR: location for top source code directory (default ROOTDIR/progs/GIT64)
# BINHOME: root directory for binaries
# BINNAME: archetecture dependent basename for bin dir
# BINDIR: directory for binaries (default BINHOME/bin/MACHTYPE) (will create if doesn't exist)
# INCLUDEPATH: include path (default PROGDIR anything else could cause a problem)
# Various directors can be overridden with environment variable or from make command
# make BINHOME=/a/path/to/binhome
#
# Base directory for code
USER =	$(shell id -u -n)
#
# Default rootdir
ifneq ($(ROOTDIR)),)
	ROOTDIR =	/Users/$(USER)
endif
$(info ROOTDIR="$(ROOTDIR)")
#
# Default root for source code
ifneq ($(PROGDIR)),)
	PROGDIR =       $(ROOTDIR)/progs/GIT64
endif
$(info PROGDIR ="$(PROGDIR)")
#
# Default location root for compiled programs
ifneq ($(BINHOME)),)
	BINHOME =		~$(USER)
endif
$(info BINHOME="$(BINHOME)")
#
# Default binary directory
ifneq ($(BINDIR)),)
	BINDIR =	$(BINHOME)/bin/$(MACHTYPE)
endif
$(info BINDIR="$(BINDIR)")
#
# Create bin dir if it doesn't exist
$(shell mkdir -p $(BINDIR))
#
# Default include path
ifneq ($(INCLUDEPATH)),)
	INCLUDEPATH =	$(PROGDIR)
endif
$(info INCLUDEPATH ="$(INCLUDEPATH)")
#
C =		gcc
CFLAGS =	'-O3  -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3  -D$(MACHTYPE) $(COMPILEFLAGS)'
#-Wunused-variable'

CCFLAGS1= '-O3'
# uncomment to debug
#CFLAGS =	'-g -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'
#
# ******** SHOULD NOT NEED TO MODIFY BELOW HERE *********
#

C =		gcc
#ROOTDIR =	/Users/ian
#PROGDIR =       $(ROOTDIR)/progs/GIT64
#INCLUDEPATH =	$(ROOTDIR)/progs/GIT64
###BINDIR =	$(IHOME)/bin/$(MACHTYPE)
CFLAGS =	-O3 -I$(INCLUDEPATH) $(COMPILEFLAGS) -m3
CCFLAGS =  '-O3  -D$(MACHTYPE) $(COMPILEFLAGS)'

#CFLAGS =	'-gbdb -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-ggdb -D$(MACHTYPE) $(COMPILEFLAGS)'

FFLAGS  =	 -g
ERS1CODEDIR =	ers1Code

RECIPES  =	cRecipes/$(MACHTYPE)-$(OSTYPE)/polint.o cRecipes/$(MACHTYPE)-$(OSTYPE)/nrutil.o \
            cRecipes/$(MACHTYPE)-$(OSTYPE)/ratint.o cRecipes/$(MACHTYPE)-$(OSTYPE)/four1.o \
            cRecipes/$(MACHTYPE)-$(OSTYPE)/svdfit.o cRecipes/$(MACHTYPE)-$(OSTYPE)/svdcmp.o \
            cRecipes/$(MACHTYPE)-$(OSTYPE)/svbksb.o cRecipes/$(MACHTYPE)-$(OSTYPE)/pythag.o \
            cRecipes/$(MACHTYPE)-$(OSTYPE)/hunt.o cRecipes/$(MACHTYPE)-$(OSTYPE)/svdvar.o

LSTRACK  =	Lstrack/$(MACHTYPE)-$(OSTYPE)/runLStrack.o
LSFIT  =	Lsfit/$(MACHTYPE)-$(OSTYPE)/readLSTiePoints.o \
			Lsfit/$(MACHTYPE)-$(OSTYPE)/computeLSFit.o \
			Lsfit/$(MACHTYPE)-$(OSTYPE)/interpLSTies.o \
			Lsfit/$(MACHTYPE)-$(OSTYPE)/readLSOffsets.o \
			Lsfit/$(MACHTYPE)-$(OSTYPE)/xyscale.o

STANDARD =	clib/$(MACHTYPE)-$(OSTYPE)/standard.o

LSTRACKDIRS =	Lstrack clib

TARGETS =	lstrack lsfit

all: $(TARGETS)

lstrack:
	@for i in ${LSTRACKDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
		make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0; \
		cd $(PROGDIR); \
		); done
		gcc    -O3 \
                Lstrack/$(MACHTYPE)-$(OSTYPE)/lstrack.o $(LSTRACK) $(STANDARD)   -lm -lgeotiff -ltiff -lfftw3f -o $(BINDIR)/lstrack


LSFITDIRS =	Lsfit Lstrack clib cRecipes

lsfit:
	@for i in ${LSFITDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0; \
			cd $(PROGDIR); \
		); done
		gcc    -O3 \
                Lsfit/$(MACHTYPE)-$(OSTYPE)/lsfit.o $(LSFIT) $(LSTRACK) $(STANDARD) $(RECIPES)   -lm -lgeotiff -ltiff -lfftw3f -o $(BINDIR)/lsfit



