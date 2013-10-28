# written by John Doug Reynolds, July 1996
# modified by Mark Pratt Oct 1996

_ostype_	= _${OSTYPE}_
__ostype__	= _${_ostype_:__=_solaris_}_

_solaris_	= ${__ostype__:__solaris__=1}
solaris		= ${_solaris_:__%__=}

_sunos4_	= ${__ostype__:__sunos4__=1}
sunos4		= ${_sunos4_:__%__=}

#------------------------------------------------------------------------------

CCC             = g++

APPS		= fitml fsamp
fitOBJ		= lcfit.o FiniteSourceAmp.o

SOURCE		= ../..
BINDIR		= ${SOURCE}/bin
NETWORK		= ${SOURCE}/network
CLIF		= ${SOURCE}/clif

INCLUDE		= -I../site -I../clif

LDPATH		= -L../clif -L../site -L/usr/openwin/lib -L. 

XLIBS		= -lX11 ${sunos4:%1=-lsuntool -lsunwindow -lpixrect}
SMLIBS		= -lplotsub -ldevices -lutils ${XLIBS}
SYSLIBS         = ${solaris:1=-lnsl} -lm


CCFLAGS		+= ${INCLUDE} -Wno-deprecated
LDFLAGS		+= ${LDPATH}

#------------------------------------------------------------------------------

#all : lib ${APPS}
all : lib

lib : libcfit.a

libcfit.a : ${fitOBJ}
	ar ruv libcfit.a $?
	@ if [ -x /usr/bin/ranlib ]; then set -x; ranlib libcfit.a; fi

${APPS} : $$@.o lib
	${CCC} -o $@ $@.o ${LDFLAGS} -lm -lsite -lclif -lcfit ${SMLIBS} -lcurses

installme: all
	cp ${APPS} ${HOME}/bin

tidy :
	rm -f ${fitOBJ} *~
	rm -fr .sb

clean : tidy
	rm -f libcfit.a ${APPS}

tar :
	tar cf - Makefile parameters *.[CHf] | compress > lcfit.tar.Z


#------------------------------------------------------------------------------

.KEEP_STATE:

.SUFFIXES: .C .o

.C.o:
	${CCC} -c ${CCFLAGS} $<
