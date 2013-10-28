/* some prototypes for SuperMongo library functions */

#ifdef __cplusplus
extern "C" {
#endif

extern int sm_device( const char *name );
extern void sm_alpha( void );
extern void sm_box( int, int, int, int );
extern void sm_conn( const float *x, const float *y, int length );
extern void sm_ctype( const char *color );
extern void sm_curs( float *x, float *y, int *c );
extern void sm_draw( double x, double y );
extern void sm_erase( void );
extern void sm_errorbar( const float *x
			, const float *y, const float *e, int k, int n );
extern void sm_gflush( void );
extern void sm_graphics( void );
extern void sm_hardcopy( void );
extern void sm_label( const char *label );
extern void sm_limits( double x1, double x2, double y1, double y2 );
extern void sm_points( const float *x, const float *y, int length );
extern void sm_ptype( const float *pp, int length );
extern void sm_putlabel( int i, const char *label );
extern void sm_relocate( double x, double y );
extern void sm_toscreen( double ux, double uy, int *sx, int *sy );
extern void sm_window( int nx, int ny, int x1, int y1, int x2, int y2 );
extern void sm_xlabel( const char *label );
extern void sm_ylabel( const char *label );

#ifdef __cplusplus
}
#endif
