CC =		gcc
CFLAGS  =	$(FLAGS) -c  -I$(INCLUDEPATH) -D_LARGEFILE64_SOURCE -D_FILE64 -D_LARGEFILE_SOURCE -DHPLARGEFILES

OBJS= lsfit.o readLSTiePoints.o readLSOffsets.o interpLSTies.o computeLSFit.o xyscale.o
USER =	$(shell id -u -n)
MACHTYPE = $(shell uname -m)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)

