#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <libsite.H>
#include <libclif.H>
#include "lcfit.H"

int field, tile, seq;
int lcs_residual( LcSet &lcs, GeneralLcFit *lcfit, int undo );
int read_ascii_set( LcSet &set,  const char *lcdir, const char *starid );

int main( int argc, char *argv[] )
{
   if( argc < 2 ) {
      printf("\nUsage:\n\t%s Field.Tile.Seq\n", argv[0] );;
      exit(1);
   }

   char *starid = argv[1];
   int m, n, p, l, st, ct, narg, dp;
   int pt, pd, dlc;
   int npb;
   int done = 0;
   int nfits = 0;
   int tamper = 0;
   int rflag = 0;
   int field;
   double t, dt;
   char c, lcdir[256], cbuf[256], outfi[256];
   FILE *outfp;

   clifpar.load( );
   
   LcSet lcs, hlcs;
   PeakFinder pf( clifpar );

   SodophotCuts sodcuts( clifpar, lcdir );
   DaophotCuts daocuts( clifpar, lcdir );
   SodophotSet ss;

   lcs.recycle( );
   /*--- Find star in local sodset cache if possible ---*/
   //   if( sscanf( argv[1], "%d.%d.%d", &field, &tile, &seq) == 3 &&
   //       getenv( "SSPATH" ) ) {
   //      sprintf( cbuf, "%s/F_%d/sodset_%d.%d",
   //	       getenv("SSPATH"), field, field, tile );
   //      if ( ss.open( cbuf ) ) {
   //	 fprintf( stderr, "%s cannot open %s\n", argv[0], cbuf );
   //      }
   //      sodcuts.clean( lcs, ss, ss.find( seq ) );
   //   } else {
   //      printf("environment variable SSPATH not set\n" );
   //   }

   /*--- Find star in local clc cache if possible ---*/
   if( getenv( "LCDIR" ) ) {
      sprintf( lcdir, "%s", getenv("LCDIR") );

      lcdir_sodlc_file( lcdir, starid, 'p', cbuf );
      if( m = ss.open( cbuf ) )
	 printf("Cannot open %s  %d\n", cbuf, m );
      else 
	  sodcuts.clean( lcs, ss, 0 );
      
      lcdir_sodlc_file( lcdir, starid, 'p', cbuf );
      if( !ss.open( cbuf ) ) sodcuts.clean( lcs, ss, 0 );

      lcdir_sodlc_file( lcdir, starid, 'o', cbuf );
      if( !ss.open( cbuf ) ) sodcuts.clean( lcs, ss, 0 );

      lcdir_sodlc_file( lcdir, starid, 'n', cbuf );
      if( !ss.open( cbuf ) ) sodcuts.clean( lcs, ss, 0 );

      lcdir_sodlc_file( lcdir, starid, 'g', cbuf );
      if( !ss.open( cbuf ) ) sodcuts.clean( lcs, ss, 0 );

      lcdir_sodlc_file( lcdir, starid, 'e', cbuf );
      if( !ss.open( cbuf ) ) sodcuts.clean( lcs, ss, 0 );

   } else	{
      printf("environment variable LCDIR not set\n" );
   }
   
   /*--- get local cached daophot lc  ---*/
   DaophotSet daoset;
   daoset.read( lcdir, starid );
   daocuts.clean( lcs, daoset );

   /*--- get local cached ascii data ---*/
   if ( read_ascii_set( lcs, lcdir, starid ) ) {
      fprintf( stderr, "problem with ascii data\n" );
      //return 1;
   }

   p = 0;
   while( p != lcs.count( ) ) {
      if( lcs[p]->nobs == 0 ) {
	 lcs[p]->erase( );
	 lcs.unhook( lcs[p] );
      } else {
	 ++p;
      }
   }

   if( lcs.count( ) == 0 ) {
      printf("Exiting: no data\n" );
      exit( SiteErr::errMissing );
   }

   for( p = 0; p < lcs.count( ); ++p ) {
      printf("Loaded %d points from %s in %s\n",
	     lcs[p]->nobs,lcs[p]->title,Passband::info(lcs[p]->filter)->name);
   }


   LcMerge( lcs );				// combine like passbands
   for( p = 0; p < lcs.count( ); ++p )		// and convert to linear units
      (*lcs[p]).convert( CommonLc::Linear );

   p = 0;
   while( p != lcs.count( ) ) {
      if( lcs[p]->nobs <= 1 ) {
	 printf("Removing %s: too few observations\n", Passband::info(lcs[p]->filter)->name);
	 lcs[p]->erase( );
	 lcs.unhook( lcs[p] );
      } else {
	 ++p;
      }
   }

   /*--- run PointFilter on LcSet ---*/
   pf.eval( lcs );
   
   /*--- Append theory Lc's for later use ---*/
   CommonLc *tlc;
   npb = lcs.count();
   for( n = 0; n < npb; ++n ) {
      lcs.append( tlc = new CommonLc );
      lc_hdr_cpy( *lcs[n], *tlc );
      tlc->units = CommonLc::Linear;
      tlc->phtpkg = CommonLc::Theory;
   }
   
   double tlo, thi, vlo, vhi;
   char instr[256], *ibuf;

   LcMoment M;
   GeneralLcFit *lcfit = 0;
   GaussianFit gfit( clifpar );
   MicroFit micro( clifpar );
   MicroBlendFit fmicro( clifpar );
   MicroBinSrcFit bsmicro( clifpar );
   MicroFinSrcFit fsmicro( clifpar );
   MicroParFit parmicro( clifpar );
   MicroBinSrcFit2 bsmicro2( clifpar );
   MicroYukFit yukmicro ( clifpar );
   gfit.init( lcs, pf );
   micro.init( lcs, pf );
   fmicro.init( lcs, pf );
   bsmicro.init( lcs, pf );
   fsmicro.init( lcs, pf );
   parmicro.init( lcs, pf );
   bsmicro2.init( lcs, pf );
   yukmicro.init( lcs, pf );
   lcfit = &micro;

   char *title = lcs[0]->title;

   /*--- Show raw data for parameter guessing ---*/
   LcPlot plot;
   plot.config( lcs );
   plot.xnot = 0.2;
   if( plot.xlimhi - plot.xlimlo > 200 )
      tlo = plot.xlimlo = plot.xlimhi - 400;
   thi = plot.xlimhi += 20;
   plot.draw( lcs );

   do {
      printf("%s\n\n", title );
      printf("f\t\tFit\n");
      printf("l\t\tLimits <[tlo thi [vlo vhi]]>\n");
      printf("d\t\tShow data <Lc # [start count]>\n");
      printf("r\t\tToggle residual plot\n");
      printf("s\t\tStatistics\n");
      printf("p\t\tPlot <device>\n" );
      printf("a\t\tRestore all data\n");
      printf("x\t\tRemove <Lc #>\n");
      printf("q\t\tQuit\n");
      printf("<space>\t\tRefresh\n" );
      fflush( stdin );
      do {
	 printf(" > ");
	 c = getc( stdin );
	 gets( ibuf = instr );
	 printf("\n");
	 done = 0;
	 c = tolower( c );
	 if( c == 'l' || c == 'd' || c == 'x' || c == 'p' ) {
	    while( isspace( *ibuf ) ) ++ibuf;
	 }
	 fflush( stdin );
	 switch( c ) {
	 case 'f' :	// select and perform fit
	    if( rflag ) {
	       printf( "Toggle out of residual mode first\n" );
	       break;
	    }
	    //	    printf("\t0\t%s\n", gfit.name );
	    printf("\t0\t%s\n", micro.name );
	    //	    printf("\t1\t%s\n", fmicro.name );
	    //	    printf("\t1\t%s\n", bsmicro.name );
	    printf("\t1\t%s\n", fsmicro.name );
	    printf("\t2\t%s\n", parmicro.name );
	    printf("\t3\t%s\n", bsmicro2.name );
	    printf("\t4\t%s\n", bsmicro.name );
	    printf("\t5\t%s\n", yukmicro.name );
	    printf("\t> " );
	    gets( ibuf = instr );
	    switch( atoi( ibuf ) ) {
	    default :
	       printf("\tInvalid fit selection\n");
	       lcfit = &micro;
	       break;
	    case 0 : if (lcfit != &micro) { ++tamper; }
	       lcfit = &micro; break;
	    case 1 : if (lcfit != &fsmicro) { ++tamper; }
	       lcfit = &fsmicro; break;
	    case 2 : if (lcfit != &parmicro) { ++tamper; }
	       lcfit = &parmicro;
	       
	       double eb, el;
	       char ecliptic[256], fts[256], *scrp;

	       // Problem is, clc title is "MACHO object 1.3450.2557" and we
	       //   only want the f.t.s.  Also in lcfit.C.
	       strcpy(ecliptic, lcs[0]->title);
	       scrp = strrchr(ecliptic, ' ');      // ' ' makes \  an int
	       *scrp = 0;
	       ++scrp;
	       strcpy(fts, scrp);
	       
	       sprintf (ecliptic, "ParallaxFit_%s_ecliptic.l", fts);
	       el = clifpar.get(ecliptic, 270.);                  		 // Use rd2lb.gmn
	       printf("\nUsing %-20s = %.7g\n", ecliptic, el );
	       
	       sprintf (ecliptic, "ParallaxFit_%s_ecliptic.b", fts);
	       eb = clifpar.get(ecliptic, -5.);                  		 // Use rd2lb.gmn
 	       printf("Using %-20s = %.7g\n\n", ecliptic, eb );
	       
	       break;
	    case 3 : if (lcfit != &bsmicro2) { ++tamper; }
	       lcfit = &bsmicro2; break;
	    case 4 : if (lcfit != &bsmicro) { ++tamper; }
	       lcfit = &bsmicro; break;
	    case 5 : if (lcfit != &yukmicro) { ++tamper; }
	       lcfit = &yukmicro; break;
	    }
	    if( !lcfit )
	       break;
	    if( tamper ) {
	       pf.eval( lcs );
	       lcfit->init( lcs, pf );
	       tamper = 0; 
	    }
	    lcfit->iopar( 1 );
	    lcfit->ifit( 0, 0, 1 );
	    lcfit->iopar( );
	    dt = (thi - tlo) * 0.001;
	    for( m = p = 0; m < lcs.count( ); ++m ) 
	       if( (*lcs[m]).phtpkg == CommonLc::Theory ) {
		  (*lcs[m]).resize( 1001 );
		  for( n = 0; n < 1001; ++n ) 
		     (*lcs[m]).date[n] = tlo + dt*n;
		  lcfit->FitLc( *lcs[m], p );
		  ++p;
	       }
	    ++nfits;
	    vhi = -( vlo = Infinity );
	    plot.config( lcs );
	    plot.xnot = 0.2;
	    plot.xlimlo = tlo;
	    plot.xlimhi = thi;
	    for( dp = 0; dp < lcs.count( ); ++dp ) {
	       plot.ylim[1][dp] += 0.3 * (plot.ylim[1][dp] - plot.ylim[0][dp]);
	       plot.ylim[0][dp] = 0;
	    }
	    plot.draw( lcs );
	    ++done;
	    break;
	 case 'r' :	// residual
	    if( rflag ) {
	       lcs_residual( lcs, lcfit, 1 );		// correct data
	       rflag = 0;
	    } else { 
	       lcs_residual( lcs, lcfit, 0 );		// take residual
	       rflag = 1;
	       vhi = -( vlo = Infinity );
	    }
	    plot.config( lcs );				// reset limits
	    plot.xnot = 0.2;
	    plot.draw( lcs );
	    break;
	 case 'l' :	// change limits
	    if( narg = sscanf( ibuf, "%lf %lf %lf %lf",
			    &tlo, &thi, &vlo, &vhi ) );

	    plot.config( lcs );
	    plot.xnot = 0.2;
	    if( narg < 2 ) {
	       tlo = plot.xlimlo;
	       thi = plot.xlimhi;
	    } else {
	       plot.xlimlo = tlo;
	       plot.xlimhi = thi;
	    }
	    
	    if( lcfit ) {
	       dt = (thi - tlo) * 0.001;
	       for( m = p = 0; m < lcs.count( ); ++m ) 
		  if( (*lcs[m]).phtpkg == CommonLc::Theory ) {
		     (*lcs[m]).resize( 1001 );
		     for( n = 0; n < 1001; ++n ) 
			(*lcs[m]).date[n] = tlo + dt*n;
		     if( rflag ) 
			for( n = 0; n < 1001; ++n ) 
			   (*lcs[m]).value[n] = 0;
		     else
			lcfit->FitLc( *lcs[m], p );
		     ++p;
		  }
	    }
	    if( narg >= 4 ) {
	       plot.ylimlo = vlo;
	       plot.ylimhi = vhi;
	    }
	    plot.draw( lcs );
	    
	    break;
	 case ' ' :	// refresh plot
	    //vhi = -( vlo = Infinity );
	    plot.draw( lcs );
	    break;
	 case 'd' :	// write Lc data
	    narg = sscanf( ibuf, "%d %d %d", &l, &st, &ct );
	    
	    if( narg == 0 || l < 0 ) {	// write files
	       for( n = pd = pt = 0; n < lcs.count( ); ++n ) {
		  dlc = (*lcs[n]).phtpkg != CommonLc::Theory;
		  sprintf( cbuf, "%s%c%d.dat", starid,
			   dlc ? 'd' : 't', dlc ? pd : pt );
		  if( dlc ) ++pd;		// data curve
		  else ++pt;			// theory curve
		  printf("writing %s\n", cbuf );
		  if( (outfp = fopen( cbuf, "w") ) == 0 ) {
		     fprintf( stderr, "Cannot open %s to write\n", cbuf );
		     continue;
		  } else {
		     (*lcs[n]).print( stdout, 0, 0 );
		     (*lcs[n]).print( outfp );
		     fclose( outfp );
		  }
	       }
	       break;
	    }
	    switch( narg ) {	// subsets are printed to screen
	    default :
	       printf( "No arguments read\n");
	       break;
	    case 1 :
	       (*lcs[l]).print( stdout );
	       ++done;
	       break;
	    case 2 :
	       (*lcs[l]).print( stdout, st );
	       ++done;
	       break;
	    case 3 : case 4 :
	       (*lcs[l]).print( stdout, st, ct );
	       break;
	    }
	    break;

	 case 'p' :	// replot
	    narg = sscanf( ibuf, "%s %s", cbuf, outfi );
	    switch( narg ) {
	    default :
	       plot.print( lcs, "postscript" );
	       break;
	    case 1 :
	       plot.print( lcs, cbuf );
	       break;
	    case 2 :
	       plot.print( lcs, cbuf, outfi );
	       break;
	    }

	 case 'a' :
	    printf("Adding %d Lc's to fit set\n", hlcs.count( ) );
	    tamper = 1;
	    while( hlcs.count( ) ) {
	       lcs.append( hlcs[0] );
	       printf("%s rehooked\n",
		      Passband::info( (*hlcs[0]).filter )->name );
	       hlcs.unhook( hlcs[0] );
	    }
	    plot.config( lcs );
	    plot.xlimlo = tlo;
	    plot.xlimhi = thi;
	    plot.draw( lcs ); 
	    break;
	 case 's' :
	    for( n = 0; n < lcs.count( ); ++n ) {
	       if( lcs[n]->phtpkg == CommonLc::Theory )
		  continue;
	       if( M.eval( *lcs[n] ) )
		  continue;
	       printf("\n%s\n", (*Passband::info( lcs[n]->filter ) ).name );
	       M.print( stdout, 2 );
	    }
	    break;
	 case 'x' :
	    n = atoi( ibuf );
	    if( n < 0 || n >= lcs.count( ) ) {
	       printf("Invalid Lc #\n" );
	       break;
	    }
	    /*--- find data curve n ---*/
	    pd = -1;
	    for( m = p = 0; m < lcs.count( ); ++m ) {
	       if( lcs[m]->phtpkg == CommonLc::Theory ) 
		  continue;
	       if( p == n ) pd = m;
	       ++p;
	    }
	    if( pd >= 0 ) { 
	       hlcs.append( lcs[pd] );
	       printf("%s unhooked\n",
		      Passband::info( (*lcs[pd]).filter )->name );
	       lcs.unhook( lcs[pd] );
	       ++tamper;
	    }
	    /*--- find theory curve n ---*/
	    pt = -1;
	    for( m = p = 0; m < lcs.count( ); ++m ) {
	       if( lcs[m]->phtpkg == CommonLc::Theory ) {
		  if( p == n ) pt = m;
		  ++p;
	       }
	    }
	    if( pt >= 0 ) { 
	       hlcs.append( lcs[pt] );
	       printf("%s unhooked\n",
		      Passband::info( (*lcs[pt]).filter )->name );
	       lcs.unhook( lcs[pt] );
	       ++tamper;
	    }
	    plot.config( lcs );
	    plot.xlimlo = tlo;
	    plot.xlimhi = thi;
	    plot.draw( lcs );
	    break;

	 case 'q' :
	    exit( 0 );
	 default:
	    break;
	 }
	 fflush( stdin );
      } while( !done );			// waiting for directions
   } while( 1 );
   return 0;
}

/*--- internal functions ---*/
int lcs_residual( LcSet &lcs, GeneralLcFit *lcfit, int undo ) {
   
   int *ti, *di;
   di = new int[lcs.count( )];		// data indices
   ti = new int[lcs.count( )];		// theory indices
   int pd = 0;
   int pt = 0;
   int m, n;
   FILE *outfp;
   char outfi[256];

   if( !lcfit ) {
      printf("No valid fit\n");
      return SiteErr::errMissing;
   }
   
   for( m = pd = pt = 0; m < lcs.count( ); ++m ) 
      if( (*lcs[m]).phtpkg == CommonLc::Theory ) ti[pt++] = m;
      else di[pd++] = m;

   if( pd != pt ) {
      printf("Unmatched LcSet d%d t%d\n", pd, pt);
      return SiteErr::errProgram;
   }
   
   for( m = 0; m < pd; ++m ) {
      (*lcs[ti[m]]).resize( (*lcs[di[m]]).nobs );
      for( n = 0; n < (*lcs[di[m]]).nobs; ++n ) 
	 (*lcs[ti[m]]).date[n] = (*lcs[di[m]]).date[n];	
      lcfit->FitLc( *lcs[ti[m]], m );
      if( undo ) 
	 for( n = 0; n < (*lcs[di[m]]).nobs; ++n ) 
	    (*lcs[di[m]]).value[n] += (*lcs[ti[m]]).value[n];
      else {		// do
	 for( n = 0; n < (*lcs[di[m]]).nobs; ++n ) {
	    (*lcs[di[m]]).value[n] -= (*lcs[ti[m]]).value[n];
	    (*lcs[ti[m]]).value[n] = 0;
	 }

	 // print out residual Lc data files
	 //	 sprintf( outfi, "%d.%d.%d_r%d.dat", field, tile, seq, m );
	 //	 if( (outfp = fopen( outfi, "w") ) == 0 ) {
	 //	    fprintf( stderr, "Cannot open %s to write\n", outfi );
	 //	    continue;
	 //	 } else {
	 //	    (*lcs[di[m]]).print( outfp );
	 //	    fclose( outfp );
	 //	 }
      }	
   }			// loop over passbands
   delete [] di;
   delete [] ti;
   return SiteErr::errNone;
}

CommonLc* read_ascii(const char *lcdir, const char *object, const char *filter)
{
   char buf[1024];
   lcdir_star_dir( lcdir, object, 'r', buf );
   sprintf( buf+strlen(buf), "/ascii.%s", filter );
   FILE *fi = fopen( buf, "r" );
   if ( !fi ) return 0;

   // initialize a new common lightcurve

   CommonLc *lc = new CommonLc( buf, 0 );
   lc->refJD =juliantime( timemacho(0) );
   lc->mag0 = NaN;
   lc->norm1 = 1;
   lc->norm0 = 0;
   lc->units = CommonLc::Magnitude;
   lc->phtpkg = CommonLc::Undefined;
   if ( strcmp( filter, "MJUO.r" ) == 0 )
      lc->filter = Passband::R_MJUO;
   else if ( strcmp( filter, "MJUO.R" ) == 0 )
      lc->filter = Passband::R_MJUO;
   else if ( strcmp( filter, "MSO74.U" ) == 0 )
      lc->filter = Passband::U_MSO74;
   else if ( strcmp( filter, "MSO74.B" ) == 0 )
      lc->filter = Passband::B_MSO74;
   else if ( strcmp( filter, "MSO74.V" ) == 0 )
      lc->filter = Passband::V_MSO74;
   else if ( strcmp( filter, "MSO74.R" ) == 0 )
      lc->filter = Passband::R_MSO74;
   else if ( strcmp( filter, "MSO74.I" ) == 0 )
      lc->filter = Passband::I_MSO74;
   else if ( strcmp( filter, "CTIO.r" ) == 0 )
      lc->filter = Passband::R_CTIO;
   else
      lc->filter = Passband::Invalid;

   // count lines in data file

   int nobs = 0;
   while ( fgets(buf,1024,fi) )
      if (buf[0] != '#' && buf[0] != '!') ++nobs;
   lc->resize( nobs );
   rewind( fi );

   // second pass: fill in the lightcurve

   nobs = 0;
   int line = 1;
   while ( nobs < lc->nobs && fgets(buf,1024,fi) ) {
      double date, mag, err;
      if (  buf[0] != '#' && buf[0] != '!' ) {

	 // added by acb...
	 char *bp;	
	 bp = strtok(buf, " ");
	 date = atof( bp );
	 bp = strtok(NULL, " ");
	 mag = atof( bp );
	 bp = strtok(NULL, " ");
	 err = atof( bp );
	 
	 lc->date[ nobs ] = date;
	 lc->value[ nobs ] = mag;
	 lc->error[ nobs ] = err;
	 lc->obsid[ nobs ] = line;
	 ++nobs;
      }
      ++line;
   }
   lc->nobs = nobs;

   fclose( fi );
   return lc;
}

int read_ascii_set( LcSet &set,  const char *lcdir, const char *starid )
{
   char path[1024];
   if( !lcdir_star_dir( lcdir, starid, 'r', path ) )
       return SiteErr::errParams;

   DIR *dir;
   struct dirent *d;
   if ( !(dir = opendir(path)) ) {
      printf("Cannot open %s\n", path );
      return SiteErr::errMissing;
   }

   while ( d = readdir(dir) )
      if ( strncmp( d->d_name, "ascii.", 6 ) == 0 ) {
	 CommonLc *clc = read_ascii( lcdir, starid, d->d_name+6 );
	 if ( clc ) set.append( clc );
      }

   closedir( dir );
   return 0;
}

