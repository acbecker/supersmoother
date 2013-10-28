# written by John Doug Reynolds, July 1996
# modified by Mark Pratt Oct 1996

# ACB
# to get debugging in the fortran code...
#
# f77 -ff2c-library -g -c -o subs.o subs.f
# make lib
# g++ -ggdb -o fazBasic fazBasic.o -L../../clif/kernel -L../../network/site -L/usr/openwin/lib -L. -L../fit -L/usr/X11R6/lib -lm -lefence -lclif -lsite -lcfit -lcphase -lf2c

CCC             = g++
#CCC             = g++ -ggdb 

_ostype_	= _${OSTYPE}_
__ostype__	= _${_ostype_:__=_solaris_}_

_solaris_	= ${__ostype__:__solaris__=1}
solaris		= ${_solaris_:__%__=}

_sunos4_	= ${__ostype__:__sunos4__=1}
sunos4		= ${_sunos4_:__%__=}

#------------------------------------------------------------------------------
APPS		= fazBasic
phaseOBJ	= rphase.o dphase.o subs.o

SOURCE		= ../..
BINDIR		= ${SOURCE}/bin
NETWORK		= ${SOURCE}/network
CLIF		= ${SOURCE}/clif

DEFINE		= -D${__ostype__} ${MACHOTYPE:%=-D__production__}

INCLUDE		= -Isite -Iclif

LDPATH		= -Lclif -Lsite -L/usr/openwin/lib -L. -Lfit -L/usr/X11R6/lib


XLIBS		= -lX11 ${sunos4:%1=-lsuntool -lsunwindow -lpixrect}
#SMLIBS		= -lplotsub -ldevices -lutils ${XLIBS}
SMLIBS		= ${XLIBS}
SYSLIBS         = ${solaris:1=-lnsl} -lm  #-lefence

CCFLAGS		+= ${DEFINE} ${INCLUDE} -Wno-deprecated
LDFLAGS		+= ${LDPATH}

#------------------------------------------------------------------------------

all : lib ${APPS}

lib : libcphase.a

libcphase.a : ${phaseOBJ}
	ar ruv libcphase.a $?
	@ if [ -x /usr/bin/ranlib ]; then set -x; ranlib libcphase.a; fi

#${APPS} : $$@.o lib
#	${CCC} -o $@ $@.o ${LDFLAGS} ${SYSLIBS} -lclif -lsite -lcfit -lcphase 

fazBasic : fazBasic.o lib
	${CCC} -o fazBasic fazBasic.o ${LDFLAGS} ${SYSLIBS} -lclif -lsite -lcfit -lcphase 

installme: all
	cp  ${APPS} ${HOME}/bin

tidy : 
	rm -f ${phaseOBJ} *~
	rm -fr .sb

clean : tidy
	rm -f libcphase.a ${APPS}
	rm -f *.o

tar :
	tar cf - Makefile parameters *.[CHf] | compress > lcphase.tar.Z


#------------------------------------------------------------------------------

.KEEP_STATE:

.SUFFIXES: .C .o

.C.o:
	${CCC} -c ${CCFLAGS} $<
.f.o: 
	f77 -ff2c-library -g -c $<
