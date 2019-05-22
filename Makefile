C =		gcc
ROOTDIR =	/Users/ian
PROGDIR =       $(ROOTDIR)/progs/GIT
INCLUDEPATH =	$(ROOTDIR)/progs/GIT
BINDIR =	$(IHOME)/bin/$(MACHTYPE)
CFLAGS =	'-O3 -I$(INCLUDEPATH) $(COMPILEFLAGS) -m3'
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



