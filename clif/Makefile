# written by John Doug Reynolds, March 1998

CCC             = g++ 

_ostype_	= _${OSTYPE}_
__ostype__	= _${_ostype_:__=_solaris_}_

#------------------------------------------------------------------------------

#OBJ		= clif.o cuts.o reimann.o utils.o almanac.o lcplot.o subs.o
OBJ		= clif.o cuts.o utils.o almanac.o subs.o

SOURCE		=  ../..
BINDIR		= ${SOURCE}/bin

DEFINE		= -D${__ostype__}

INCLUDE		= -I../site

CCFLAGS		+= ${DEFINE} ${INCLUDE} -Wno-deprecated

#------------------------------------------------------------------------------

all lib : libclif.a

libclif.a : ${OBJ}
	ar ruv libclif.a $?
	@ if [ -x /usr/bin/ranlib ]; then set -x; ranlib libclif.a; fi

tidy :
	rm -f ${OBJ} *~
	rm -fr .sb

clean : tidy
	rm -f libclif.a

tar :
	tar cf - Makefile parameters *.[CHf] | compress > kernel.tar.Z

#------------------------------------------------------------------------------

INSTALLANY = installinclude installlib install
${INSTALLANY} := INCLUDE = -I${MACHO}/include -I/usr/macholocal/include
${INSTALLANY} := LIBS_DIR = -L${MACHO}/lib

installinclude :
	cvscheck
	@ if cmp -s libclif.H ${MACHO}/include/libclif.H; then :; else\
		echo install -m 660 libclif.H ${MACHO}/include;\
		install -m 660 libclif.H ${MACHO}/include;\
	fi

installlib : installinclude lib
	install -m 660 libclif.a ${MACHO}/lib

install : installlib

offinst : lib

#------------------------------------------------------------------------------

.KEEP_STATE:

.SUFFIXES: .C .o

.C.o:
	${CCC} -c ${CCFLAGS} $<
