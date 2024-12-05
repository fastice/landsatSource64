CC =		gcc
CFLAGS  =	$(FLAGS) -c -I$(INCLUDEPATH)
FFLAGS  =	 -g

OBJS=	gammln.o poidev.o  sort.o polint.o nrutil.o ratint.o four1.o \
        svdfit.o svbksb.o  svdcmp.o pythag.o hunt.o svdvar.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)





