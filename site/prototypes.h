/* written by John Doug Reynolds, July 1997 */

#ifndef _prototypes_h
#define _prototypes_h

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUG__
#if defined __sunos4__
   void		bcopy( char *b1, char *b2, int length );
   int		munmap( caddr_t addr, int len );
#elif defined __solaris__
   int		ftruncate( int fd, off_t length );
   int		gethostname( char *name, int namelen );
#elif defined __linux__
   char*	cuserid( char *buf );
#endif
   int		select( int width, fd_set *readfds, fd_set *writefds
			, fd_set *exceptfds, struct timeval *timeout );
   char*	inet_ntoa( struct in_addr in );
#endif /* __GNUG__ */

#ifdef __sunos4__
   char*	strerror( int err );
#endif

#ifdef __cplusplus
}
#endif

#endif /* _prototypes_h */
