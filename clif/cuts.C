/* written by John Doug Reynolds, March 1998 */

#include <Macho.h>
#include <libsite.H>
#include "libclif.H"

#include <sys/param.h>
#include <sys/stat.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>

class RunningMean {
public:
   RunningMean( int window );			// set buffer length
   ~RunningMean( void );			// delete buffer
   operator float() const;			// return running mean
   int valid( void ) const;			// true after first add
   RunningMean&	operator += ( float );		// add to running mean
private:
   float	*buf;
   float	total;
   int		count;
   int		top;
   int		max;
};

/*--- BadObsList methods ----------------------------------------------------*/
// written by Mark Pratt

BadObsList::BadObsList( void ) {
   nrange = maxn = 0;
   min = max = 0;
   *filename = 0;
   filter = Passband::Invalid;
}

BadObsList::~BadObsList( void ) {
   delete [] max;
   delete [] min;
}

int BadObsList::reset( const int N ) {
   filter = Passband::Invalid;		// invalidate list
   lastr = nrange = 0;
   *filename = 0;
   if( N < maxn )
      return SiteErr::errNone;
   maxn = 0;
   delete [] max;
   delete [] min;
   if( !(max = new int[N]) || !(min = new int[N]) )
      return SiteErr::errMemory;
   maxn = N;
   return SiteErr::errNone;
}

int BadObsList::reorder( void ) {
   if( nrange <= 0 )
      return SiteErr::errNone;

   int *omin = min;
   int *omax = max;
   if( !( min = new int[maxn] ) || !( max = new int[maxn] ) )
      return SiteErr::errMemory;

   int *idx = new int[nrange];
   int r, m, n;
   isort( idx, omin, nrange );
   for( m = n = 0; n < nrange; ++n ) {
      r = idx[n];

      /*--- consistency check ranges ---*/
      if( omin[r] > omax[r] ||			// improper range
	 (m && omax[r] < max[m-1] ) ) {		// redundant range
	 continue;
      } else if( m && omin[r] <= max[m-1] ) {	// merge overlapping ranges
	 max[m-1] = omax[r];
	 continue;
      }
      min[m] = omin[r];
      max[m] = omax[r];
      ++m;
   }
   delete [] omin;
   delete [] omax;
   delete [] idx;
   return nrange = m;
}
   
int	BadObsList::match( const short fltr, const char *str )
{
   char *fstr;
   switch ( fltr ) {
   case Passband::R_MACHO:	fstr = "macho.r";	break;
   case Passband::V_MACHO:	fstr = "macho.v";	break;
   case Passband::I_CTIO:	fstr = "ctio.i";	break;
   case Passband::R_CTIO:	fstr = "ctio.r";	break;
   case Passband::V_CTIO:	fstr = "ctio.v";	break;
   case Passband::B_CTIO:	fstr = "ctio.b";	break;
   case Passband::U_CTIO:	fstr = "ctio.u";	break;
   case Passband::I_MJUO:	fstr = "mjuo.i";	break;
   case Passband::R_MJUO:	fstr = "mjuo.r";	break;
   case Passband::V_MJUO:	fstr = "mjuo.v";	break;
   case Passband::B_MJUO:	fstr = "mjuo.b";	break;
   case Passband::U_MJUO:	fstr = "mjuo.u";	break;
   case Passband::I_UTSO:	fstr = "utso.i";	break;
   case Passband::R_UTSO:	fstr = "utso.r";	break;
   case Passband::V_UTSO:	fstr = "utso.v";	break;
   case Passband::B_UTSO:	fstr = "utso.b";	break;
   case Passband::U_UTSO:	fstr = "utso.u";	break;
   case Passband::I_WISE:	fstr = "wise.i";	break;
   case Passband::R_WISE:	fstr = "wise.r";	break;
   case Passband::V_WISE:	fstr = "wise.v";	break;
   case Passband::B_WISE:	fstr = "wise.b";	break;
   case Passband::U_WISE:	fstr = "wise.u";	break;
   default: return 0;
   }
   int flen = strlen( fstr );

   int c, n;
   char tok[64], *p;
   for ( c = 0; sscanf( str+c, "%s%n", tok, &n ) != EOF; c += n ) {
      int tlen = strlen( tok );
      for ( p = tok; *p = tolower(*p); ++p );
      if ( tok[0] == '.' && tlen == 1 )
	 return 1;		// '.' matches anything
      else if ( tok[0] == '.' && strncmp( tok, fstr+flen-tlen, tlen ) == 0 )
	 return 1;		// '.YYY' matches 'XXX.YYY'
      else if ( tok[tlen-1] == '.' && strncmp( tok, fstr, tlen ) == 0 )
	 return 1;		// 'XXX.' matches 'XXX.YYY'
      else if ( strcmp( tok, fstr ) == 0 )
	 return 1;		// 'XYZ' matches 'XYZ'
   }
   return !c;			// empty str matches any filter
}

int BadObsList::load( const short fltr, const char *fname, const char *lcdir,
		      const char *object ) {

   if ( fltr == Passband::Invalid || !fname || ( object && !lcdir ) ) {
      reset( 0 );
      return SiteErr::errParams;
   }

   FILE *fi = 0;
   char buf[1024];

   /*--- create filename and open badobs file ---*/
   if( object ) {		// use only object-specific badobs 
      if( lcdir_star_dir( lcdir, object, 'r', buf ) ) {	// object dir exists
	 sprintf( buf+strlen(buf)  , "/%s", fname );
	 if( strcmp( buf, filename ) == 0 && filter == fltr)
	    return SiteErr::errNone; 	 		// file already loaded
	 else if(fi = fopen( buf, "r" ) ) {
	    strcpy( filename, buf );
	 }
      }	// else fail
   } else if( strcmp( filename, fname ) == 0 && filter == fltr ) {
      return SiteErr::errNone;				// file already loaded
   } else if( fi = fopen( fname, "r" ) ) {		// check locally
      strcpy( filename, fname );
   } else if( lcdir && *fname != '/' ) {		// check in lcdir
      sprintf( buf, "%s/%s", lcdir, fname );
      if( strcmp( buf, filename ) == 0 && filter == fltr)
	 return SiteErr::errNone; 	 		// file already loaded
      if( fi = fopen( buf, "r" ) )
	 strcpy( filename, buf );
   } // else fail
   if( !fi ) {
      filter = Passband::Invalid;	// invalidate list
      //fprintf( stderr, "BadObsList cannot read %s\n", filename );
      return SiteErr::errMissing;
   } 

   /*--- find size by reading file - this should be improved! ---*/
   int N = 0;
   while( fgets( buf, 1024, fi ) )
      ++N;
   rewind( fi ); 
   reset( N );
   filter = fltr;
   //fprintf( stderr, "BadObsList reading %s\t\t%d\n", filename, N );
   

   /*--- read matching ranges from badobs file ---*/
   int disordered = 0;
   N = 0;
   while( fgets( buf, 1024, fi ) ) {
      int n, len;
      char tok[64];
      if( *buf == '#' ) continue;		// comment line

      sscanf( buf, "%s%n", tok, &len );		// scan obsid range
      if( !match( filter, buf+len ) ) {		// filter match remaining buf
	 continue;
      }
      
      /*--- parse obsid range ---*/ 
      if ( sscanf( tok, "%ld-%ld", min+N, max+N ) != 2 )
	 if ( tok[0] == '-' && strlen(tok) == 1 )
	    min[N] = -(max[N] = INT_MAX);	// '-' matches all obs
	 else if ( tok[0] == '-' && sscanf( tok+1, "%ld", max+N ) == 1 )
	    min[N] = -INT_MAX;			// '-XXX' matches -infty - XXX 
	 else if ( sscanf( tok, "%ld%n", min+N, &n ) == 1 )
	    max[N] =
	       tok[n] == '-' ? INT_MAX : min[N]; // 'XXX-' matches XXX - infty
	 else
	    continue;
      
      /*--- consistency check range ---*/
      if( min[N] > max[N] ) {			// improper range
	 continue;					// skip
      } 
      if( N )
	 if( min[N] < min[N-1] ) {		// min out of order
	    ++disordered;				// retain, sort later
	 } else if( max[N] < max[N-1] ) {	// redundant range
	    continue;					// skip
	 } else if( min[N] <= max[N-1] ) {	// merge overlapping ranges
	    max[N-1] = max[N];				// merge w. last
	    continue;
	 }
      ++N;
   }
   nrange = N;
   if( disordered ) nrange = reorder( );
   return SiteErr::errNone;
}

int BadObsList::ok( const int obsid, const int ordered ) {
   /*------------------------------------------------------------------------
     ok assumes properly ordered and distinct list of ranges 
     ordered ==> look up sequentially from last call
     !ordered ==> look up obsid with binary search
    *----------------------------------------------------------------------*/
     
   if( nrange <= 0 )
      return 1;

   if( ordered ) {	// sequential search
      int r = lastr;
      while( min[r] > obsid && --r >= 0 ) ;
      if( r == -1 ) return 1;		// out of range
      while( max[r] < obsid && ++r < nrange ) ;
      if( r == nrange )	return 1;	// out of range
      lastr = r;
      if( min[r] > obsid ) return 1;
      return 0;
   } else {		// binary search
      int a = 0, m, z = nrange - 1;
      while( a <= z ) {
	 m = (a + z) / 2;
	 if( max[m] < obsid ) a = m + 1;
	 else if( min[m] > obsid ) z = m - 1;
	 else break;
      }
      lastr = m;
      return a > z ? 1 : 0;
   } 
}

/*--- SodFlags methods ------------------------------------------------------*/

SodFlags::SodFlags( void )
{
   maxnobs = nobs = 0;
   crwd = chi2 = calib = 0;
   airmass = seeing = 0;
}

SodFlags::~SodFlags( void )
{
   delete [] crwd;
   delete [] chi2;
   delete [] calib;
   delete [] airmass;
   delete [] seeing;
}

int SodFlags::reset( int size )
{
   nobs = 0;
   if( maxnobs >= size ) return 0;
      
   maxnobs = 0;
   delete [] crwd;
   delete [] chi2;
   delete [] calib;
   delete [] airmass;
   delete [] seeing;

   if ( !(crwd = new short[size])	||
	!(chi2 = new short[size])	||
	!(calib = new short[size])	||
	!(airmass = new float[size])	||
	!(seeing = new float[size])
	)
      return SiteErr::errMemory;

   maxnobs = size;
   return 0;
}

/* --- SodophotCuts methods ------------------------------------------------ */

SodophotCuts::SodophotCuts( const char *lcdir_, const char *active_ )
{
   badobsFile = 0;
   verbose = doCuts = publicCuts = heliocentric = calibrate = 0;
   typeMod20 = type9ok = 0;
   typeMax = chi2Max = crwdMax = cosmMax = missMax = 0;
   errMin = errMax = magMin = chi2magMax = seeingMax = 0;
   errScale = errExtra2 = 0;
   refJD = NaN;

   badobsloaded = 0;
   active = duplicate(active_);
   lcdir = duplicate(lcdir_);
}

SodophotCuts::SodophotCuts( const ClifParameters &params
			    , const char *lcdir_, const char *active_ )
{
   load( params );
   badobsloaded = 0;
   active = duplicate(active_);
   lcdir = duplicate(lcdir_);
}

SodophotCuts::~SodophotCuts( void )
{
   if ( badobsFile ) delete [] badobsFile;
   if ( active ) delete [] active;
   if ( lcdir ) delete [] lcdir;
}

void	SodophotCuts::set( const char *lcdir_, const char *active_ )
{
   if ( active ) delete [] active;
   if ( lcdir ) delete [] lcdir;
   active = duplicate(active_);
   lcdir = duplicate(lcdir_);
}

void	SodophotCuts::load( const ClifParameters &params )
{
   refJD = params.get( "SodophotCuts.refJD", 2448623.5 );
   doCuts = params.get( "SodophotCuts.doCuts", 1 );
   verbose = params.get( "SodophotCuts.verbose", 1 );
   publicCuts = params.get( "SodophotCuts.publicCuts", 0 );
   heliocentric = params.get( "SodophotCuts.heliocentric", 0 );
   calibrate = params.get( "SodophotCuts.calibrate", 0 );
   typeMod20 = params.get( "SodophotCuts.type.mod20", 1 );
   type9ok = params.get( "SodophotCuts.type.9ok", 1 );
   typeMax = params.get( "SodophotCuts.type.max", 6 );
   chi2Max = params.get( "SodophotCuts.chi2.max", 140 );
   crwdMax = params.get( "SodophotCuts.crowding.max", 169 );
   cosmMax = params.get( "SodophotCuts.cosmic.max", 1 );
   missMax = params.get( "SodophotCuts.missing.max", 50 );
   errMin = params.get( "SodophotCuts.err.min", 0 );
   errMax = params.get( "SodophotCuts.err.max", 3 );
   magMin = params.get( "SodophotCuts.mag.min", -15 );
   chi2magMax = params.get( "SodophotCuts.chi2mag.max", -7 );
   seeingMax = params.get( "SodophotCuts.seeing.max", 8 );
   errScale = params.get( "SodophotCuts.error.scale", 1 );
   errExtra2 = pow(params.get("SodophotCuts.error.extra",0.014)*2.5/log(10),2);
   SoP = params.get( "SodophotCuts.SideOfPier", 0 );	// 1=W, -1=E, 0=both

   if ( badobsFile ) delete [] badobsFile;
   const char *str = params["SodophotCuts.badobs.file"];
   badobsFile = str && strcmp(str,"0") ? duplicate(str) : 0;
}

int	SodophotCuts::loadbadobs( const char *object )
{
   badobsloaded = 0;
   if ( ! badobsFile ) return 0;	// badobs disabled

   char fname[MAXPATHLEN+1];
   struct stat sbuf;

   /* badobs.file set to '+':	look for badobs.macho in lcdir, but let
				object-specific version override if existing */

   if ( strcmp(badobsFile,"+") == 0 ) {
      if ( ! lcdir_star_dir(lcdir,object,'r',fname) ) {
	 if ( verbose )
	    if ( ! lcdir ) fprintf( stderr, "loadbadobs: no lcdir\n" );
	    else if ( ! object ) fprintf( stderr, "loadbadobs: no object\n" );
	    else fprintf( stderr, "loadbadobs: filename too long\n" );
	 return SiteErr::errParams;
      }
      strcat(fname,"/badobs.macho");
      if ( stat(fname,&sbuf) && errno == ENOENT ) {
	 sprintf(fname,"%s/badobs.macho",lcdir);
	 if ( stat(fname,&sbuf) && errno == ENOENT ) {
	    if ( verbose )
	       fprintf( stderr, "badobs.macho not found in %s\n", lcdir );
	    return SiteErr::errMissing;
	 }
      }
   }

   /* otherwise:	look for badobsFile in cwd first, then in lcdir */

   else {
      strcpy(fname,badobsFile);
      if ( stat(fname,&sbuf) && errno == ENOENT )
	 if ( lcdir && *badobsFile != '/' ) {
	    sprintf(fname,"%s/%s",lcdir,badobsFile);
	    if ( stat(fname,&sbuf) && errno == ENOENT ) {
	       if ( verbose )
		  fprintf( stderr, "%s not found in %s\n", badobsFile, lcdir );
	       return SiteErr::errMissing;
	    }
	 }
	 else {
	    if ( verbose )
	       if ( lcdir ) fprintf( stderr, "%s not found\n", badobsFile );
	       else fprintf(stderr,"%s not found and no lcdir\n",badobsFile);
	    return SiteErr::errMissing;
	 }
   }

   SiteErr err = badobs.load(fname,verbose);
   if ( ! err ) badobsloaded = 1;
   return err;
}

static int	initlc( CommonLc *lc, const char *title, int nobs )
{
   if ( !lc  ||  !title ) return 1;
   if ( lc->resize( nobs ) ) return 1;
   lc->title = new char [ strlen(title) + 1 ];
   if ( !lc->title ) return 1;
   strcpy( lc->title, title );
   return 0;
}

static void	unpack( SodophotSet::StarObs &s
			, short m, short e, short f1, short f2, short f3 )
{
   s.mag = (m==INVALID_MAG) ? NaN : 0.001 * m;
   s.err = (e < 0) ? (e==INVALID_SHORT) ? NaN : e : 0.001 * e;

   if ( f1 != INVALID_SHORT ) {
      s.dsky = ((int) f1 + 32200) / 100 - 322;
      s.type = f1 - 100 * s.dsky;
   }
   else s.dsky = s.type = 999;
   
   if ( f2 != INVALID_SHORT ) {
      int tmp = f2 + 32768;
      s.crwd = tmp >> 8;
      s.chi2 = tmp & 255;
   }
   else s.crwd = s.chi2 = 999;
   
   if ( f3 != INVALID_SHORT ) {
      int tmp = f3 + 32768;
      s.mpix = tmp >> 8;
      s.cosm = tmp & 255;
   }   
   else s.mpix = s.cosm = 999;
}

int	SodophotCuts::ok( const SodophotSet::StarObs &s, double seeing ) const
{
   if ( !doCuts ) return !isnan(s.err);

   if ( isnan(s.mag+s.err+seeing) || s.type == 999 || s.crwd == 999
	|| s.chi2 == 999 || s.mpix == 999  || s.cosm == 999
      )
      return 0;

   int type = typeMod20 ? s.type % 20 : s.type;

   if ( s.mag < magMin || s.err < errMin || s.err > errMax || s.crwd > crwdMax
	|| s.mpix > missMax || s.cosm > cosmMax || seeing > seeingMax
	|| ( type > typeMax && !(type9ok && type == 9) )
	|| ( s.mag > chi2magMax && s.chi2 > chi2Max )
      )
      return 0;

   return 1;
}

int	SodophotCuts::clean( LcSet &set, SodLcMetaDb &db ) {
   if ( !db.star || heliocentric || calibrate ) return 1;

   char object[64], title[64];
   if ( db.hdr->iwjsver < 8 ) {
      int lot = 1000 * db.field + db.hdr->chunk;
      int item = (db.hdr->tileno[0] - db.star->tile << 15) + db.star->seq;
      sprintf( object, "%d.%d", lot, item );
   }
   else
      sprintf( object, "%d.%d.%d", db.field, db.star->tile, db.star->seq );
   sprintf( title, "MACHO object %s", object );

   if ( badobsFile &&
	( ! badobsloaded || strcmp(badobsFile,"+") == 0 ) &&
	loadbadobs(object)
	)
      return 1;
   
   int nobs = db.obscount();
   CommonLc *rlc = set.newlc();
   CommonLc *blc = set.newlc();
   if ( initlc(rlc,title,nobs) || initlc(blc,title,nobs) ) {
      set.recycle( rlc );
      set.recycle( blc );
      return 1;
   }

   if ( db.hdr->iwjsver < 8 ) {
      const double rad_mas = M_PI / 648e6;
      rlc->ra = blc->ra = db.star->ra * rad_mas;
      rlc->dec = blc->dec = db.star->dec * rad_mas;
   }

   rlc->refJD = blc->refJD = refJD;
   rlc->mag0 = blc->mag0 = 0;
   rlc->norm1 = blc->norm1 = 1;
   rlc->norm0 = blc->norm0 = 0;

   rlc->units = blc->units = CommonLc::Magnitude;
   rlc->filter = Passband::R_MACHO;
   blc->filter = Passband::V_MACHO;
   rlc->phtpkg = blc->phtpkg = CommonLc::oSodophot;

   int two_years_ago = int( time(0) - 63072000 );
   int badobsIndex = 0;

   rlc->nobs = blc->nobs = 0;
   for ( db.reset(db.star->tile,db.star->seq); db.obsindex != -1; db.next() ) {
      SodLcDb::Observation &obs = db.obs[ db.obsindex ];
      SodLcDb::Photometry &phot = db.star->phot[ db.obsindex ];
      SodophotSet::StarObs uphot;

      if ( publicCuts && two_years_ago < obs.time ) continue;
      if( SoP && obs.ipier != SoP ) continue;

      double date = juliantime( timemacho(obs.time) ) - refJD;

      int badamps = badobsFile ? badobs.amps(obs.obsno,badobsIndex) : 0;

      if ( ! (badamps & 0x00FF) ) {
	 short rerr = ( phot.rerr == -99 ) ? INVALID_SHORT : phot.rerr;
	 unpack(uphot,phot.rmag,rerr,phot.rflag1,phot.rflag2,phot.rflag3);
	 if ( ok( uphot, obs.seeing[0] ) ) {
	    int i = rlc->nobs;
	    rlc->date[i] = date;
	    rlc->obsid[i] = obs.obsno;
	    rlc->value[i] = uphot.mag;
	    double err = errScale * uphot.err;
	    rlc->error[i] = sqrt( err * err + errExtra2 );
	    ++rlc->nobs;
	 }
      }

      if ( ! (badamps & 0xFF00) ) {
	 short berr = ( phot.berr == -99 ) ? INVALID_SHORT : phot.berr;
	 unpack(uphot,phot.bmag,berr,phot.bflag1,phot.bflag2,phot.bflag3);
	 if ( ok( uphot, obs.seeing[1] ) ) {
	    int i = blc->nobs;
	    blc->date[i] = date;
	    blc->obsid[i] = obs.obsno;
	    blc->value[i] = uphot.mag;
	    double err = errScale * uphot.err;
	    blc->error[i] = sqrt( err * err + errExtra2 );
	    ++blc->nobs;
	 }
      }
   }
   if ( set.append(rlc) < 0 || set.append(blc) < 0 ) return 1;
   else return 0;
}

int	SodophotCuts::clean( const SodophotSet &sod, int index
			     , class CommonLc *rlc, class CommonLc *blc
			     , class SodFlags *rsfp, class SodFlags *bsfp )
			     
{
   char object[32], title[64];
   sod.starid( index, object );
   sprintf( title, "MACHO object %s", object );
   class SodCalibrator calibrator;
   class RunningMean meancolor(20);
   const double rad_mas = M_PI / 648e6;
   double ra = sod.star[index].ra * rad_mas;
   double dec = sod.star[index].dec * rad_mas;
   int rwcid = sod.star[index].wcid[Red];
   int bwcid = sod.star[index].wcid[Blue];
   int two_years_ago = int( time(0) - 63072000 );
   int badobsIndex = 0;
   int nobs = sod.hdr->count.o;

   short trmag = sod.star[index].mag[Red];
   short tbmag = sod.star[index].mag[Blue];
   float tbmr = trmag != INVALID_MAG && tbmag != INVALID_MAG
      ? 0.001 * (tbmag - trmag) : NaN;

   if( rsfp ) rsfp->reset( nobs );
   if( bsfp ) bsfp->reset( nobs );

   if ( calibrate > 0 ) {
      if ( rwcid >= 128 && bwcid >= 128 ) {
	 if ( verbose )
	    fprintf( stderr, "%s has no template information\n", object );
	 return 1; // probably a corrupted sodset
      }
      if ( abs(calibdb.field) != sod.hdr->field ) {
	 SiteErr err = calibdb.read(sod.hdr->field,lcdir,active);
	 if ( err && verbose ) err.perror( "calibdb.read" );
	 if ( err ) return 1;
      }
      int redwestchunk = rwcid < 128 ? rwcid : (~ChunkToken(West,bwcid)).id();
      if ( calibdb.get(calibrator,redwestchunk) ) {
	 if ( verbose )
	    fprintf( stderr, "%s is probably corrupted\n", object );
	 return 1;
      }
   }

   if ( badobsFile &&
	( ! badobsloaded || strcmp(badobsFile,"+") == 0 ) &&
	loadbadobs(object)
	)
      return 1;

   if ( publicCuts && two_years_ago < sod.hdr->incept ) return 1;
   if ( heliocentric && isnan(ra+dec) ) return 1;

   short phtpkg = (sod.hdr->source == 'p')
      ? CommonLc::pSodophot : CommonLc::oSodophot;

   if ( rlc ) {
      if ( initlc(rlc,title,nobs) ) return 1;
      rlc->ra = ra;
      rlc->dec = dec;
      rlc->refJD = refJD;
      rlc->mag0 = 0;
      rlc->norm1 = 1;
      rlc->norm0 = 0;
      rlc->units = CommonLc::Magnitude;
      rlc->filter = calibrate < 2 ? Passband::R_MACHO : Passband::R_Std;
      rlc->phtpkg = phtpkg;
      rlc->nobs = 0;
   }

   if ( blc ) {
      if ( initlc(blc,title,nobs) ) return 1;
      blc->ra = ra;
      blc->dec = dec;
      blc->refJD = refJD;
      blc->mag0 = 0;
      blc->norm1 = 1;
      blc->norm0 = 0;
      blc->units = CommonLc::Magnitude;
      blc->filter = calibrate < 2 ? Passband::V_MACHO : Passband::V_Std;
      blc->phtpkg = phtpkg;
      blc->nobs = 0;
   }

   for ( int o = 0; o < nobs; ++o ) {
      SodophotSet::Observation &obs = sod.obs[o];
      SodophotSet::StarObs rphot, bphot;
      double seeing, date;

      if ( publicCuts && two_years_ago < obs.date ) continue;
      if( SoP && obs.pierside != SoP ) continue;

      if ( heliocentric ) {
	 if ( obs.exposure == INVALID_SHORT ) continue;
	 double jd = juliantime( obs.date - obs.exposure / 2 );
	 date = ::heliocentric(jd,ra,dec) - refJD;
      }
      else date = juliantime(obs.date) - refJD;

      int coords = 0, badamps = 0;
      if ( badobsFile &&
	   (badamps = badobs.amps(obs.obsid,badobsIndex)) &&
	   badamps != 0xFFFF
	   )
	 coords = 'a'; // compute amp coordinates only if necessary

      sod.unpack( rphot, index, o, Red, coords );
      sod.unpack( bphot, index, o, Blue, coords );

      int rflag = isnan(rphot.mag) ? 0 : 1;
      int bflag = isnan(bphot.mag) ? 0 : 1;
      SideOfPier sop = (SideOfPier) obs.pierside;
      short calibFlag = 0;
      float bmr;

      if ( calibrate > 0 )
	 if ( rflag && bflag ) {
	    bmr = bphot.mag - rphot.mag;
	    calibFlag = 4;
	 }
	 else if ( meancolor.valid() ) {
	    bmr = meancolor;
	    calibFlag = 3;
	 }
	 else if ( ! isnan(tbmr) ) {
	    bmr = tbmr;
	    calibFlag = 2;
	 }
	 else {
	    bmr = 0;
	    calibFlag = 1;
	 }

      if (rlc && !(badamps & (rphot.ampid == -1 ? 0x00FF : 1<<rphot.ampid))) {
	 seeing = rwcid < 128 ? sod.chnk[rwcid][o].seeing : NaN;
	 if ( ok( rphot, seeing ) ) {
	    int i = rlc->nobs;
	    rlc->date[i] = date;
	    rlc->obsid[i] = obs.obsid;
	    if ( calibrate < 2 || !rflag )
	       rlc->value[i] = rphot.mag;
	    else
	       rlc->value[i] = calibrator.R(rphot.mag,bmr);
	    double err = errScale * rphot.err;
	    rlc->error[i] = sqrt( err * err + errExtra2 );
	    if( rsfp ) {
	       rsfp->crwd[i] = rphot.crwd;
	       rsfp->chi2[i] = rphot.chi2;
	       rsfp->calib[i] = calibFlag;
	       rsfp->seeing[i] = seeing;
	       rsfp->airmass[i] = obs.airmass;
	    }
	    ++rlc->nobs;
	    ++rflag;
	 }
      }

      if (blc && !(badamps & (bphot.ampid == -1 ? 0xFF00 : 1<<bphot.ampid))) {
	 seeing = bwcid < 128 ? sod.chnk[bwcid][o].seeing : NaN;
	 if ( ok( bphot, seeing ) ) {
	    int i = blc->nobs;
	    blc->date[i] = date;
	    blc->obsid[i] = obs.obsid;
	    if ( calibrate < 1 || !bflag )
	       blc->value[i] = bphot.mag;
	    else if ( calibrate == 1 )
	       blc->value[i] = calibrator.B(sop,bphot.mag,bmr);
	    else
	       blc->value[i] = calibrator.V(sop,bphot.mag,bmr);
	    double err = errScale * bphot.err;
	    blc->error[i] = sqrt( err * err + errExtra2 );
	    if( bsfp ) {
	       bsfp->crwd[i] = bphot.crwd;
	       bsfp->chi2[i] = bphot.chi2;
	       bsfp->calib[i] = calibFlag;
	       bsfp->seeing[i] = seeing;
	       bsfp->airmass[i] = obs.airmass;
	    }
	    ++blc->nobs;
	    ++bflag;
	 }
      }

      if ( calibrate > 0 && rflag == 2 && bflag == 2 ) meancolor += bmr;
   }
   if ( rlc && rsfp ) rsfp->nobs = rlc->nobs;
   if ( blc && bsfp ) bsfp->nobs = blc->nobs;
   return 0;
}

int	SodophotCuts::clean( LcSet &set, const SodophotSet &sod, int index
			     , class SodFlags *rsfp, class SodFlags *bsfp )
{
   CommonLc *rlc = set.newlc();
   CommonLc *blc = set.newlc();
   if ( clean(sod,index,rlc,blc,rsfp,bsfp) ) {
      set.recycle( rlc );
      set.recycle( blc );
      return 1;
   }
   return (set.append(rlc) < 0 || set.append(blc) < 0);
}

/* --- DaophotCuts methods ------------------------------------------------- */

DaophotCuts::DaophotCuts( const char *path )
{
   refJD = NaN;
   doCuts = normMin = 0;
   sharpMin = sharpMax = aspectMax = quadErrMax = dMagChi2Max = 0;
   errScale = errExtra2 = 0;
   lcdir = path ? strdup(path) : 0;
   badobsfile = 0;
}

DaophotCuts::DaophotCuts( const ClifParameters &params, const char *path )
{
   lcdir = path ? strdup(path) : 0;
   load( params );
}

DaophotCuts::~DaophotCuts( void )
{
   free( lcdir );
   free( badobsfile );
}

void	DaophotCuts::load( const ClifParameters &params )
{
   refJD = params.get( "DaophotCuts.refJD", 2448623.5 );
   doCuts = params.get( "DaophotCuts.doCuts", 1 );
   normMin = params.get( "DaophotCuts.min.norm", 3 );
   sharpMin = params.get( "DaophotCuts.sharpness.min", -0.5 );
   sharpMax = params.get( "DaophotCuts.sharpness.max", 1 );
   aspectMax = params.get( "DaophotCuts.psfAspect.max", 0.4 );
   quadErrMax = params.get( "DaophotCuts.quadErr.max", 3 );
   dMagChi2Max = params.get( "DaophotCuts.dMagChi2.max", 4 );
   errScale = params.get( "DaophotCuts.error.scale", 1 );
   errExtra2 = pow(params.get("DaophotCuts.error.extra",0)*2.5/log(10),2);
   const char *p = params["DaophotCuts.badobs.file"];
   badobsfile = p && strcmp(p,"0") ? strdup(p) : 0;
   badobs.reset( 0 );
}

void	DaophotCuts::set( const char *lcdir_ )
{
   if ( lcdir ) delete [] lcdir;
   lcdir = duplicate(lcdir_);
}

int	DaophotCuts::clean( CommonLc &lc, const DaophotLc &ilc )
{
   char title[32];
   int field, tile, seqn;
   if ( sscanf( ilc.object, "%d.%d.%d", &field, &tile, &seqn ) == 3 )
      sprintf( title, "MACHO object %s", ilc.object );
   else
      strcpy( title, ilc.object );
   if ( initlc( &lc, title, ilc.nobs ) ) return 1;

   lc.refJD = refJD;
   lc.mag0 = 25;
   lc.norm1 = 1;
   lc.norm0 = 0;

   lc.units = CommonLc::Magnitude;
   lc.phtpkg = CommonLc::Daophot;
   lc.filter = Passband::info(ilc.filter)->id;
   lc.nobs = 0;

   if( badobsfile && *badobsfile == '+' )
      badobs.load( lc.filter, "badobs.gman", lcdir, ilc.object );
   else
      badobs.load( lc.filter, badobsfile, lcdir );
   
   for ( int o = 0; o < ilc.nobs; ++o ) {
      DaophotLc::Photometry &phot = ilc.phot[o];

      int ok = !isnan( phot.mag + phot.err );
      if ( ok  &&  doCuts ) {
	 double err, aspect;
	 err = sqrt( phot.dMagErr * phot.dMagErr + phot.err * phot.err );
	 aspect = fabs( (phot.xFWHM-phot.yFWHM) / (phot.xFWHM+phot.yFWHM) );
	 if ( aspect > aspectMax
	      || err > quadErrMax
	      || phot.dMagChi2 > dMagChi2Max
	      || phot.sharpness < sharpMin
	      || phot.sharpness > sharpMax
	      || phot.normstars < normMin
	      || !badobs.ok( phot.obsid, 1 ) )
	    ok = 0;
      }

      if ( ok ) {
	 int i = lc.nobs;
	 lc.date[i] = juliantime(phot.date) - refJD;
	 lc.obsid[i] = phot.obsid;
	 lc.value[i] = phot.mag - lc.mag0;
	 double err = phot.err * errScale;
	 lc.error[i] = sqrt( err * err + errExtra2 );
	 ++lc.nobs;
      }
   }
   return 0;
}

int	DaophotCuts::clean( LcSet &set, const DaophotSet &dao )
{
   if ( !dao.count() ) return 0;

   for ( int b = 0; b < dao.count(); ++b ) {
      CommonLc *lc = set.newlc();
      if ( clean(*lc,*dao[b]) ) {
	 set.recycle( lc );
	 return 1;
      }
      else if ( set.append(lc) < 0 ) return 1;
   }
   return 0;
}

/* --- RunningMean methods ------------------------------------------------- */

RunningMean::RunningMean( int window )
{
   max = window > 0 ? window : 0;
   buf = max > 0 ? new float [max] : 0;
   total = 0;
   count = 0;
   top = -1;
}

RunningMean::~RunningMean( void )
{
   if ( buf ) delete [] buf;
}

RunningMean::operator float() const
{
   return count ? total / count : NaN;
}

int	RunningMean::valid( void ) const
{
   return count > 0;
}

class RunningMean&	RunningMean::operator +=( float f )
{
   if ( max > 0 )
      if ( count == max ) {
	 if ( ++top >= max ) top = 0;
	 total += f - buf[top];
	 buf[top] = f;
      }
      else {
	 top = count++;
	 total += buf[top] = f;
      }
   return *this;
}
