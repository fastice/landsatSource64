CC =		gcc
CFLAGS =	$(FLAGS) -D$(MACHTYPE) -c -I$(INCLUDEPATH)

OBJS =	standard.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)
.KEEP_STATE:
