# written by John Doug Reynolds, February 1998

_ostype_	= _${OSTYPE}_
__ostype__	= _${_ostype_:__=_solaris_}_

#------------------------------------------------------------------------------

CCC             = g++

HDR		= libsite.H prototypes.h

OBJ		= site.o etime.o SodophotSet.o SodCalib.o badobs.o StarNote.o\
		  DaophotLc.o FitsHeader.o SodLcDb.o LcSummary.o fits_header.o

SOURCE		= ../..
BINDIR		= ${SOURCE}/bin

INCLUDE		= 
DEFINE		= 
CCFLAGS		+= ${DEFINE} ${INCLUDE} -Wno-deprecated

#------------------------------------------------------------------------------

all lib : libsite.a
libsite.a : ${OBJ}
	ar ruv libsite.a $?
	@ if [ -x ranlib ]; then set -x; ranlib libsite.a; fi

client install_client : lib

server install_server : lib

tidy :
	rm -f ${OBJ}
	rm -fr .sb

clean : tidy
	rm -f libsite.a

tar :
	tar cf - Makefile *.[chCH1] | compress > site.tar.Z

#------------------------------------------------------------------------------

INSTALLANY = installinclude installlib install
${INSTALLANY} := INCLUDE = -I${MACHO}/include
${INSTALLANY} := LIBS_DIR = -L${MACHO}/lib

installinclude :
	cvscheck
	@ for h in ${HDR}; do\
	   if cmp -s $$h ${MACHO}/include/$$h; then :; else\
		echo install -m 660 $$h ${MACHO}/include;\
		install -m 660 $$h ${MACHO}/include;\
	   fi;\
	done

installlib : installinclude libsite.a
	install -m 660 libsite.a ${MACHO}/lib

install : installlib

#------------------------------------------------------------------------------

.KEEP_STATE:

.SUFFIXES: .C .o

.C.o:
	${CCC} -c ${CCFLAGS} $<

.c.o:
	${CC} -c ${CFLAGS} $<
