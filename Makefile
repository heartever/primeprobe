PROGRAMS=newattack
SPYOBJS=spy.o
NEWATTACKOBJS=newattack.o
OBJS=cachemap.o probe.o pageset.o cpuid.o comm.o clock.o analyse.o 
ALLOBJS=${OBJS} ${SPYOBJS} ${NEWATTACKOBJS}

#INCLUDES=-Ilibcpuid
#LDFLAGS=-Llibcpuid
#LDLIBS=-lm -lcpuid
INCLUDES=
LDFLAGS=
LDLIBS=-lm 
#CFLAGS=-O2 -g -D_FILE_OFFSET_BITS=64 -std=gnu99 -Wall
CFLAGS=${INCLUDES} -g -std=gnu99





all: ${PROGRAMS}

newattack: ${NEWATTACKOBJS} analyse.o 
	g++ -O2 ${NEWATTACKOBJS} analyse.o -o newattack
	
${NEWATTACKOBJS}: newattack.c
	g++ -c newattack.c -o newattack.o

clean:
	rm -f ${ALLOBJS} ${PROGRAMS} *~ core.*
