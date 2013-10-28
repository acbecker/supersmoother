/*
FITS header utility package

version 1.0 by Doug Reynolds, November 1993	<jdr@lensing.physics.ucsb.edu>

This package is intended to provide a reasonably robust
and convenient toolkit for manipulating FITS format headers.

General Notes:
o  card numbers start with card number zero
o  [card_num, comment, header_size, err] are optional parameters
o  NULL may be passed instead of a valid address for any optional parameter
o  card_num, when supplied, MUST BE INITIALIZED with the search starting card
o  always remember to deallocate headers using free_header (never use free)

documentation for each routine is provided in the source file: fits_header.c
*/

/* -------------------------------------------------------------------------*/

#ifndef _fits_header_h
#define _fits_header_h

#include <sys/types.h>		/* to get size_t and off_t */
#include <stdio.h>		/* to get FILE */

#ifdef __cplusplus
extern "C" {
#endif

enum { h_no_err = 0, h_logic_err, h_format_err, h_mem_err, h_file_err };

typedef char CARD[80];
struct header_control {
    int		end;
    CARD*	H;
};
typedef struct header_control *HEADER;

HEADER	new_header();
HEADER	copy_header( HEADER srce_hdr );
HEADER	read_header( int fd, size_t* header_size, int* err );
HEADER	fread_header( FILE* fi, size_t* header_size, int* err );
int	write_header( int fd, HEADER hdr );
int	fwrite_header( FILE* fi, HEADER hdr );
void	free_header( HEADER hdr );
int	enlarge_header( HEADER hdr, int cards );
off_t	header_offset( HEADER hdr );

char*	h_strerror( int error_code );

char*	h_get_card( HEADER hdr, int card );
char*	h_find_card( HEADER hdr, char* key, int* card_num );
char*	h_card_search( HEADER hdr, char* str, int* card_num );

int	h_read_card( HEADER hdr, char* key, char* type, void* value );
int	h_write_card( HEADER hdr, int card, char* key,
		      char* type, void* value, char* comment );
int	h_change_card( HEADER hdr, char* key,
		       char* type, void* value, char* comment );
int	h_append_card( HEADER hdr, char* key,
		       char* type, void* value, char* comment );

#ifdef __cplusplus
}
#endif

#endif /* _fits_header_h */
