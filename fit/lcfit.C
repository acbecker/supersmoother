#include <math.h>
#include <ctype.h>
#include <string.h>
#include "lcfit.H"

/*--- GeneralLcFit Methods --------------------------------------------------*/

GeneralLcFit::GeneralLcFit( ClifParameters &CP ) : cp( CP ) {
   sprintf( name, "Invalid Fit" );
   maxnobs = maxnpar = Tnobs = npar = mfit = niter = 0;	// int
   DOF = PDOF = 0;					// int
   ptype = cnstrt = dof = pdof = 0;			// int *
   pname = 0;						// char **
   Chi2 = PChi2 = 0;					// double
   par = norm = ptmp = fvalue = 0;		 	// double *
   chi2 = pchi2 = mxdchi2 = mxdcdate = 0;		// double *
   covar = alpha = dvdp = 0;				// double **
   lamda = NaN;						// double
   
   /*--- convergence criteria ---*/
   lamda = maxdchi2 = maxfchi2 = NaN;
   maxit = nchi = 0;

}

GeneralLcFit::~GeneralLcFit( void ) {
   free_matrix( covar );
   free_matrix( alpha );
   delete [] name;
   while( lcs.count( ) )
      lcs.unhook( lcs[0] );
   delete [] par;
   delete [] ptmp;
   delete [] ptype;
   free_matrix( pname );
   delete fvalue;
   delete [] chi2;
   delete [] pchi2;
   delete [] mxdchi2;
   delete [] mxdcdate;
   delete [] cnstrt;
   delete [] dof;
   delete [] pdof;
   free_matrix( dvdp );
}

int GeneralLcFit::load( LcSet &s, int inpar ) {
   
   int n, p;
   int smaxnobs = 0;
   int merr = 1, lerr = 0;
   double dJD;

   /*--- Note: GeneralLcFit only keeps copies of Lc pointers, it
    ---- doesn't maintain any Lc memory
    ---*/
   while( lcs.count( ) )			// clear out old data
      lcs.recycle( lcs[0] );

   npar = inpar;
   Tnobs = 0;
   if( s.count( ) == 0 )
      return SiteErr::errMissing;
   
   for( p = 0; p < s.count( ); ++p )
      if( (*s[p]).phtpkg != CommonLc::Theory ) {
	 Tnobs += (*s[p]).nobs;
	 smaxnobs = smaxnobs < (*s[p]).nobs ? (*s[p]).nobs : smaxnobs;
	 lcs.append( s[p] );
      }
   if( npar > maxnpar || smaxnobs > maxnobs ) {
      free_matrix( dvdp );
      matrix( dvdp, npar > maxnpar ? npar : maxnpar,
	      maxnobs > smaxnobs ? maxnobs : smaxnobs );
   }
   if( npar > maxnpar ) {
      free_matrix( covar );
      free_matrix( alpha );
      delete [] par;
      delete [] ptmp;
      delete [] ptype;
      free_matrix( pname );

      matrix( covar, npar, npar );
      matrix( alpha, npar, npar );
      par = new double[npar];
      ptmp = new double[npar];
      ptype = new int[npar];
      matrix( pname, npar, 64 );
   }
   if( smaxnobs > maxnobs ) {
      delete [] fvalue;
      fvalue = new double[smaxnobs];
   }
   if( lcs.count( ) > maxnpb ) {
      delete [] chi2;
      delete [] pchi2;
      delete [] mxdchi2;
      delete [] mxdcdate;
      delete [] cnstrt;
      delete [] dof;
      delete [] pdof;
      chi2 = new double[lcs.count( )];
      pchi2 = new double[lcs.count( )];
      mxdchi2 = new double[lcs.count( )];
      mxdcdate = new double[lcs.count( )];
      cnstrt = new int[lcs.count( )];
      dof = new int[lcs.count( )];
      pdof = new int[lcs.count( )];
   }
   
   maxnobs = maxnobs < smaxnobs ? smaxnobs : maxnobs;
   maxnpar = maxnpar < npar ? npar : maxnpar;
   maxnpb = maxnpb < lcs.count( ) ? lcs.count( ) : maxnpb;

   /*--- align all lc's with date of lcs[0] ---*/
   if( lcs.count( ) ) {
      refJD = (*lcs[0]).refJD;
      for( p = 0; p < lcs.count( ); ++p ) 
	 if( (*lcs[p]).refJD != refJD ) {
	    dJD = (*lcs[p]).refJD - refJD;
	    for( n = 0; n < (*lcs[p]).nobs; ++n )
	       (*lcs[p]).date[n] -= dJD;
	 }
   }
   for( n = 0; n < npar; ++n ) {			// Is this necessary?
      par[n] = NaN;
      ptype[n] = Guess;					// auto guess
   }
   return SiteErr::errNone;
}

int GeneralLcFit::init_par( int parno, double parval, int type) {
   if( parno < 0 || parno >= npar )
      return SiteErr::errParams;
   par[parno] = parval;
   ptype[parno] = type;
   return SiteErr::errNone;
}

int GeneralLcFit::iopar( int prompt, int stats ) {
   int n;
   char instr[256], *ibuf;
   char vbuf[256];
   fflush( stdin );
   printf("%s Fit\n", name );
   if( niter > 0 && stats ) {			// print fitpars
      for( n = 0; n < lcs.count( ); ++n ) {
	 printf("\t%-20s  %3d measurements, %2d constraints\n",
		Passband::info(lcs[n]->filter)->name, lcs[n]->nobs, cnstrt[n]);
	 printf("%-20s: %.5g / %d d.o.f. = %.5g\n", "chi2",
		chi2[n], dof[n], chi2[n] / dof[n] );
	 printf("%-20s: %.5g / %d d.o.f. = %.5g\n", "peak chi2",
		pchi2[n], pdof[n], pchi2[n] / pdof[n] );
	 printf("%-20s: %.5g on %.6g\n", "hi-sigma pt",
		mxdchi2[n], mxdcdate[n] );
      }
      printf("\tCombined statistics\n" );
      printf("%-20s: %.5g / %d d.o.f. = %.5g\n", "chi2",
	     Chi2, DOF, Chi2 / DOF );
      printf("%-20s: %.5g / %d d.o.f. = %.5g\n", "peak chi2",
	     PChi2, PDOF, PChi2 / PDOF );
      printf("\n");
   }
   printf("Fit Parameters\n" );
   if( prompt ) {
      printf("Specify fit parameter values and types\n" );
      printf("<value>\t\tfit parameter\n" );
      printf("<:[value]>\tfit parameter\n" );
      printf("<=[value]>\tfixed parameter\n" );
      printf("<*>\t\tauto-guess (when available)\n" );
      printf("<cr>\t\tretain current parameter status\n\n" );
   }

   for( n = 0; n < npar; ++n ) {
      printf("%-20s", pname[n] );
      switch ( ptype[n] ) {
      default :
	 printf(" -- type error --\n" );
	 return SiteErr::errParams;
      case Fixed :
	 sprintf(vbuf, "= %.5g", par[n] );
	 break;
      case Fit :
	 if( niter > 0 )
	    sprintf(vbuf, ": %.7g (%.7g)", par[n], sqrt(covar[n][n]) );
	 else
	    sprintf(vbuf, ": %.7g", par[n] );
	 break;
      case Guess :
	 sprintf(vbuf, "* %.7g", par[n] );
	 break;
      }

      if( norm && (n+lcs.count()) >= norm-par ) {		// print adu normalization too
	 
	 if( norm && n >= norm-par ) {
	    sprintf( vbuf+strlen(vbuf), " <%.6g>", lcs[n-(norm-par)]->norm1);
	 }
	 else {
	    sprintf( vbuf+strlen(vbuf), " <%.6g>", par[n] * lcs[(n+lcs.count())-(norm-par)]->norm1);
	 }
      }	
      
      if( prompt ) {				// accept new param values
	 printf("%-24s ? ", vbuf );
	 gets( ibuf = instr );
	 while( isspace( *ibuf ) && *ibuf) ++ibuf;

	 /*--- parse input string ---*/
	 switch (*ibuf) {
	 case '=' :					// Fixed
	    ++ibuf;					// skip '='
	    while( isspace( *ibuf ) && !*ibuf )
	       ++ibuf;
	    if( *ibuf == 0 )				// keep old par
	       ptype[n] = Fixed;
	    else if( isdigit( *ibuf ) || *ibuf == '.' || *ibuf == '-' )
	       init_par( n, atof( ibuf ), Fixed );	
	    else {
	       printf("Input format error\n" );
	       --n;
	    }
	    break;
	 case '*' :					// Guess
	    init_par( n, NaN, Guess );
	    break;
	 case 0 :					// no input
	    break;
	 case ':' :					// change to Fit
	    ptype[n] = Fit;
	    ++ibuf;					// must preceed Fit
	    while( isspace( *ibuf ) )
	       ++ibuf;
	    if( !*ibuf )				// retain old value
	       break;
	 default :					// Fit
	    if( isdigit( *ibuf ) || *ibuf == '.' || *ibuf == '-' ) {
	       init_par( n, atof( ibuf ), Fit );
	    } else {
	       printf("Input format error\n" );
	       --n;
	    }
	    break;
	 }
      } else {					// just show current value
	 printf("%s\n", vbuf );
      }
   }
   printf("\n");
   return SiteErr::errNone;
}

int GeneralLcFit::mrqfit( int verbose ) {

   int c, n, conv;
   double* ochi2 = new double[nchi];
   double deltachi2;
   int mrqerr = 0;
   char msg[256];
   char instr[256], *ibuf = instr;

   niter = conv = 0;
   lamda = -1;					// initialize
   for( c = 0; c < nchi; ++c )
      ochi2[c] = Infinity;

   if( mrqerr = GLFmrqmin( ) )
      return mrqerr;

   for( n = 0; n < maxit; ++n ) {			// iterate 
	 
      ochi2[n%nchi] = Chi2;
      if( mrqerr = GLFmrqmin( ) )
	 break;
      if( verbose ) {
	 sprintf(msg, "GeneralLcFit::mrqfit iteration %d  %.5g (%.5g) %.5g",
		 n, Chi2, ochi2[n%nchi] - Chi2,
		 Chi2/(Tnobs - npar) );
	 show_dvector( par, npar, msg );
	 if( verbose > 1 )
	    show_dvector( ochi2, nchi, "Chi2 History" );
      }

      /*--- check for convergence over entire chi2 history buffer ---*/
      for( conv = 1, c = 0; c < nchi; ++c ) {
	 conv &= ( deltachi2 = ochi2[c] - Chi2 ) >= 0; 
	 conv &= deltachi2  < maxdchi2;
	 conv &= deltachi2 / Chi2 < maxfchi2;
      }
      if( conv )
	 break;
   }
   delete [] ochi2;
   if( !conv && verbose ) 
      printf("GeneralLcFit::mrqfit: did not converge\n" );
   lamda = 0.0;				       // tidy up and
   if( mrqerr = GLFmrqmin( ) )		       // sort covariance matrix
      return mrqerr;
   /* niter > 0 ==> fit accomplished, niter == maxit ==> not converged yet */
   niter = n;
   return SiteErr::errNone;
}

int GeneralLcFit::renormalize( ) {
   int m, n, p;
   if( !norm )
      return SiteErr::errNone;
   for( p = 0; p < lcs.count( ); ++p ) {
      (*lcs[p]).renorm(norm[p], 0.0);	
      for( n = 0; n < npar; ++n ) {		// renorm covariance matrix?
	 m = npar - lcs.count( ) + p;		// for baseline flux->amp
	 covar[m][n] /= norm[p];               	// correct row
	 covar[n][m] /= norm[p]; 		// correct column
      }
      norm[p] = 1.0;	           		
   }
   return SiteErr::errNone;
}

int GeneralLcFit::ifit( int prompt, int verbose, int renorm ) {
   char c;
   int parfault = 0;
   if( verbose )
      iopar( );
   for( int p = 0; p < npar; ++p )
      parfault |= isnan( par[p] ) || ptype[p] == Guess;
   if( parfault && prompt )
      iopar( 1 );
   else if( parfault )
      return SiteErr::errParams;
   
   do {
      Chi2 = Infinity;	
      mrqfit( verbose );
      if( renorm )
	 renormalize( );
      if( verbose )
	 iopar( );
      if( prompt ) {
	 printf("\nReset parameters (y/n)? ");
	 c = getc( stdin );
	 fflush( stdin );
	 printf("\n");
	 prompt = ( tolower( c ) == 'y' );
	 if( prompt )
	    iopar( prompt );
      }
   } while( prompt );
   return SiteErr::errNone;
}

/*--- GeneralLcFit NR stuff ---*/

int GeneralLcFit::GLFmrqmin( void ) {
   int j,k,l,m;
   int mrqminerr = 0;
   static double ochi2, *beta, *dpar, **onedpar, *ftmp;

   if (lamda < 0.0) {				// initialize fit
      for( mfit = j = 0; j < npar; j++ ) 	// count unfixed params
	 if (ptype[j])
	    mfit++;
      beta = new double[npar];
      matrix( onedpar, mfit, 1 );
      dpar = *onedpar;				// slightly flakey
      lamda=0.001;
      if( mrqminerr = GLFmrqcof( alpha, beta ) )
	 return mrqminerr;
      ochi2=Chi2;
      for( j = 0; j < npar; j++ )
	 ptmp[j] = par[j];
   }
  for( j = l = 0; l < npar; l++  ) if( ptype[l] ) {		// j <= l
      for( k = m = 0; m < npar; m++ ) if( ptype[m] ) {		// k <= m
	 covar[j][k] = alpha[j][k];
	 k++;
      }
      covar[j][j] = alpha[j][j] * (1.0 + lamda);
      onedpar[j][0] = beta[j];
      ++j;
   }
   if( mrqminerr = cgaussj( covar, mfit, onedpar, 1) )	// est dpar
      return mrqminerr;
   if (lamda == 0.0) {					// last call
      //if( mrqminerr = GLFmrqcof( covar, dpar ) )	// optional??
      // return mrqminerr;
      GLFcovsrt( );
      free_matrix( onedpar );
      delete [] beta;
      return mrqminerr;
   }

   /*--- try new parameters ---*/
   for( j = l = 0; l < npar; l++ ) if (ptype[l])		// j <= l
      ptmp[l] = par[l] + dpar[j++];
   ftmp = par; par = ptmp; ptmp = ftmp;		// switch to new pars
   if( norm )
      norm = par + npar - lcs.count( );		// baseline parameters (npar-3)
   if( mrqminerr = GLFmrqcof( covar, dpar ) )
      return mrqminerr;
   if( Chi2 < ochi2 ) {					// success
      lamda *= 0.1;
      ochi2 = Chi2;
      for( j = l = 0; l < npar; l++ ) if( ptype[l] ) {		// j <= l
	 for( k = m = 0; m < npar; m++ )  if( ptype[m] ) {	// k <= m
	    alpha[j][k] = covar[j][k];				// transcribe
	    k++;
	 }
	 beta[j] = dpar[j];					// transcribe
	 ++j;
      }
   } else {							// failure
      ftmp = par; par = ptmp; ptmp = ftmp;			// switch back
      if( norm )
	 norm = par + npar - lcs.count( ); // norms are last lcs.count() pars
      lamda *= 10.0;						// change step
      Chi2=ochi2;						// retain best
   }
   return mrqminerr;
}

int GeneralLcFit::GLFmrqcof( double **Alpha, double *Beta ) {

   /*--------------------------------------------------------------------
     This is the standard NR GLFmrqcof modified to for the GeneralLcFit
     class.  The indexing has been changed to C standard 0 offset.
     Calculates Alpha and Beta defined on NR pp 682,3
     *-------------------------------------------------------------------*/

   int i, j, k, l, m, p;
   double wt, ierr2, dv, c, b, dc, pk;
   CommonLc *lc;

   for( j = 0; j < mfit; j++ ) {
      for( k = 0; k <= j; k++ )
	 Alpha[j][k]=0.0;
      Beta[j]=0.0;
   }
   Chi2 = PChi2 = 0;
   DOF = PDOF = 0;
   for( p = 0; p < lcs.count( ); ++p ) {		// passbands
      chi2[p] = pchi2[p] = 0;
      dof[p] = pdof[p] = 0;
      mxdchi2[p] = -Infinity;
      lc = lcs[p];
      eval( lc->date, fvalue, dvdp, lc->nobs, p );	// evaluate Lc
      pk = pkmin * (b =  norm ? norm[p] : 1.0);
      for( i = 0; i < lc->nobs; i++ ) {			// obs
	 ierr2 = 1.0 / (lc->error[i] * lc->error[i]);
	 dv = lc->value[i] - fvalue[i];
	 for( j = l = 0; l < npar; l++ )		// j is subset index
	    if ( ptype[l] ) {
	       wt = dvdp[l][i] * ierr2;	 
	       for( k = m = 0; m <= l; m++ )		// k is subset index
		  if ( ptype[m] )
		     Alpha[j][k++] += wt * dvdp[m][i];
	       Beta[j] += dv * wt;
	       ++j;
	    }
	 chi2[p] += c = dv * dv * ierr2;
	 ++dof[p];
	 /*--- Accumulate fit statistics on last pass ---*/
	 //if( lamda == 0 ) {
	 if( 1 ) {
	    if( fvalue[i] > pk ) { 		// peak chi^2
	       pchi2[p] += c;
	       ++pdof[p];
	    }
	    dc = (lc->value[i] - b);
	    dc = dc * dc * ierr2;
	    if( dc > mxdchi2[p] ) {
	       mxdchi2[p] = dc;
	       mxdcdate[p] = lc->date[i];
	    }
	 }
      }
      Chi2 += chi2[p];			// cumulative chi^2
      PChi2 += pchi2[p];		// cumulative peak chi^2
      DOF += dof[p];
      PDOF += pdof[p];
      dof[p] -= cnstrt[p];
      pdof[p] -= cnstrt[p];
   }
   DOF -= npar - nfxpar;
   PDOF -= npar - nfxpar;
   for( j = 1; j < mfit; j++ )
      for( k = 0; k < j; k++ )		// leave diagonal alone
	 Alpha[k][j] = Alpha[j][k]; 	// fill in symetric side
   return SiteErr::errNone;
}

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
int GeneralLcFit::GLFcovsrt( ) {
   int i,j,k;
   double temp;
   /*--- expand compacted covariance matrix ---*/
   for( i = mfit; i < npar; i++ ) // zero pad matrix
      for( j = 0; j < i; j++ )
	 covar[i][j] = covar[j][i] = 0.0;
   
   k = mfit-1;
   for( j = npar-1; j >= 0; j-- )
      if( ptype[j] ) {		// fit par index
	 for( i = 0; i < npar; i++ ) // swap cols
	    SWAP(covar[i][k],covar[i][j]);
	 for( i = 0; i < npar; i++ ) // swap rows
	    SWAP(covar[k][i],covar[j][i]);
	 k--;
      }
   return 0;
}

int cgaussj(double **a, int n, double **b, int m) {
   int i, icol, irow, j, k, l, ll;
   double big, dum, pivinv, temp;
   
   int *indxc = new int[n];
   int *indxr = new int[n];
   int *ipiv = new int[n];
   
   for( j = 0; j < n; j++ )
      ipiv[j] = 0;
   
   for( i = 0; i < n; i++ ) {
      big = 0.0;
      for( j = 0; j < n; j++ ) { // find largest elt
	 if( ipiv[j] != 1 ) {
	    for( k = 0; k < n; k++ ) {
	       if( ipiv[k] == 0 ) {
		  if( fabs( a[j][k] ) >= big ) {
		     big=fabs( a[j][k] );
		     irow=j;
		     icol=k;
		  }
	       } else if( ipiv[k] > 1 ) {
		  printf("cgaussj: piv > 1\n");
		  return SiteErr::errMath;
	       }
	    }
	 }
      }
      ++(ipiv[icol]);
      if( irow != icol ) { // swap rows in a & b
	 for( l = 0; l < n; l++ )  SWAP( a[irow][l], a[icol][l] );
	 for( l = 0; l < m; l++ )  SWAP( b[irow][l], b[icol][l] );
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if( a[icol][icol] == 0.0) {
	 printf("cgaussj: zero diagonal\n");
	 return SiteErr::errMath;
      }
      pivinv = 1.0 / a[icol][icol];
      a[icol][icol] = 1.0;
      for( l = 0; l < n; l++ ) a[icol][l] *= pivinv;
      for( l = 0; l < m; l++ ) b[icol][l] *= pivinv;
      for( ll = 0; ll < n; ll++ ) {
	 if(ll != icol) {
	    dum = a[ll][icol];
	    a[ll][icol] = 0.0;
	    for(l = 0; l < n; l++ ) a[ll][l] -= a[icol][l]*dum;
	    for(l = 0; l < m; l++ ) b[ll][l] -= b[icol][l]*dum;
	 }
      }
   }
   for( l = n - 1; l >= 0; l-- ) {
      if( indxr[l] != indxc[l] )
	 for( k = 0; k < n; k++ )
	    SWAP( a[k][indxr[l]], a[k][indxc[l]] );
   }
   delete [] ipiv;
   delete [] indxr;
   delete [] indxc;
   return SiteErr::errNone;
}

/*--- matrix utilities ---*/
void free_matrix( char*** h ) {
   if( !h ) return;
   if( !*h ) return;
   delete [] **h;
   delete [] *h; *h = 0; // explicitely set matrix invalid
   return;
}
void free_matrix( char** &m ) { free_matrix( (char***) &m ); };
void free_matrix( short** &m ) { free_matrix( (char***) &m ); };
void free_matrix( int** &m ) { free_matrix( (char***) &m ); };
void free_matrix( float** &m ) { free_matrix( (char***) &m ); };
void free_matrix( double** &m ) { free_matrix( (char***) &m ); };

int matrix( char*** m, int r, int c, int nbytes ) {
   if( !m )
      return SiteErr::errParams; // you fool!
   if( *m )
      free_matrix( m );
   *m = new char* [r];
   if( !*m )
      return SiteErr::errMemory;
   (*m)[0] = new char[ r * c * nbytes ];
   if( !(*m)[0] ) {
      delete [] *m;
      return SiteErr::errMemory;
   }
   for( int n = 1; n < r; ++n ) 
      (*m)[n] = (*m)[0] + n * c * nbytes;
   return SiteErr::errNone;
}
int matrix( char** &m, int r, int c ) {
   return matrix( (char***) &m, r, c, sizeof( char ) );
}
int matrix( short** &m, int r, int c ) {
   return matrix( (char***) &m, r, c, sizeof( short ) );
}
int matrix( int** &m, int r, int c ) {
   return matrix( (char***) &m, r, c, sizeof( int ) );
}
int matrix( float** &m, int r, int c ) {
   return matrix( (char***) &m, r, c, sizeof( float ) );
}
int matrix( double** &m, int r, int c ) {
   return matrix( (char***) &m, r, c, sizeof( double ) );
}

void show_dmatrix( double **m, int R, int C, char *name ) {
   
   if( name )
      printf("%s:\n", name );
   
   if( !m ) {
      printf("Null matrix handle\n");
      return;
   }
   for( int r = -1; r < R; ++r ) {
      if( r >= 0 ) printf("\n%2d: ", r );
      else printf("    ");
      if( r >= 0 && !m[r] ) {
	 printf("Null row pointer\n" );
	 continue;
      }
      for( int c = 0; c < C; ++c ) {
	 if( r == -1) printf("%10d: ", c );
	 else printf("%11.5g ", m[r][c] );
      }
   }
   printf("\n");
}

void show_dvector( double *v, int N, char *name ) {
   int n;
   if( name )
      printf("%s:\n", name );
   
   if( !v) {
      printf("Null pointer\n");
      return;
   }
   for( n = 0; n < N; ++n ) printf("%10d: ", n );
   printf("\n");
   for( n = 0; n < N; ++n ) printf("%11.5g ", v[n] );
   printf("\n");
}

/*--- MicroFit Methods ------------------------------------------------------*/

int MicroFit::init( LcSet &s ) {
   int err;
   int p;
   
   /*--- resize MicroFit and point lcs at data ---*/
   if( err = load( s, 3 + 2 * s.count() ) ) // 3 params + normalizations
      return err;
   npar = 3 + 2 * lcs.count( );	// some Lc's may have been cut
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( ); // baseline parameters
   
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
   
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Standard Microlensing" );
   sprintf( pname[0], "U_min" );
   sprintf( pname[1], "T_zero" );
   sprintf( pname[2], "T^hat" );
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 5;	// constraints on each passband
      sprintf( pname[p+3], "%-12s f",
	       Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+3+lcs.count()], "%-12s",
	       Passband::info( (*lcs[p]).filter )->name );
   }
   return 0;
}

int MicroFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;
   
   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      init_par( p+3, 1.0, Fixed );
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+3+lcs.count(), M.average( ), Fit );
   }

   double A = pf.rhidev.max;
   init_par( 0, sqrt(2) * sqrt(A/sqrt(A*A-1) - 1), Fit);
   init_par( 1, pf.rhidev.mdate, Fit);
   init_par( 2, pf.rhidev.duration, Fit);
   
   return err;
}

int MicroFit::eval( double *date, double *v, double **dVdP, int nobs, int pb ) {
   
   int t, n;
   int npb = lcs.count();

   double u, A, dAdu, u2, squ4, dt, f, B;
   
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = umin
   // par[1] = t0
   // par[2] = that
   // par[3..3+npb] = f 	lensed fraction of passband i-3
   // par[3+npb..2*npb+3] = B	baseline norm of passband i-3-npb
   
   /* u = sqrt(umin^2 + (2(t-t0)/that)^2 */
   /* where that = v/(2 Re) */
   
   if( pb >= lcs.count( ) )
      for( t = 0; t < nobs; ++t )
	 v[t] = 1.0;
    
   f = pb >= 0 ? par[3+pb] : 1;        // lensed fraction
   B = pb >= 0 ? par[3+npb+pb] : 1;    // baseline normalization

   for( t = 0; t < nobs; ++t ) { // observations
      u = 2 * (dt = (date[t] - par[1]) ) / par[2];
      u = sqrt(u2 = (par[0] * par[0] + u*u) );
      v[t] = B * (f * (A = (u2 + 2) / (u * (squ4 = sqrt(u2+4) ) ) ) + 1 - f);
      
      if( dVdP ) {
	 dAdu = 2 / squ4 - (u2 + 2) * (squ4 + u2 / squ4 ) / (u2 * ( u2 + 4 ) );
	 dVdP[0][t] = B * f * dAdu * par[0] / u; 			// dFdu
	 dVdP[1][t] = -B * f * dAdu * 4 * dt / (u * par[2] * par[2]);	// dFdt0
	 dVdP[2][t] = dt * dVdP[1][t] /  par[2];			// dFdthat
	 
	 for(n = 0; n < npb; ++n) {
	    dVdP[n+3][t] = pb == n ? B*(A-1) : 0;          	  // Lensed fraction f
	    dVdP[n+npb+3][t] = pb == n ? f*(A-1) + 1 : 0;         // Baseline B
	 }
      }
   }
   return SiteErr::errNone;
}


/*--- MicroFinSrcFit Methods ------------------------------------------------*/
#include "FiniteSourceAmp.H"
 
int MicroFinSrcFit::init( LcSet &s ) {
   int p;
   int err = 0;
 
   /*--- resize MicroFinSrcFit and point lcs at data ---*/
   if( err = load( s, 4 + 2 * s.count() ) )         // 4 params + fracs & normalizations
      return err;
   npar = 4 + 2 * lcs.count( );                     // some Lc's may have been cut
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( );            // points to unlensed baseline
 
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
 
 
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Finite Source Microlensing" );
   sprintf( pname[0], "U_min" );
   sprintf( pname[1], "T_zero" );
   sprintf( pname[2], "T^hat" );
   sprintf( pname[3], "U_star" );
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 6;                    // constraints on each passband
      sprintf( pname[p+4], "%-12s f",
	       Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+4+lcs.count()], "%-12s",
	       Passband::info( (*lcs[p]).filter )->name );
   }
   return 0;
}
 
int MicroFinSrcFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;

   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      init_par( p+4, 1.0, Fixed );
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+4+lcs.count(), M.average( ), Fit );
   }

   double A = pf.rhidev.max;
   init_par( 0, sqrt(2) * sqrt(A/sqrt(A*A-1) - 1), Fit);
   init_par( 1, pf.rhidev.date, Fit);
   init_par( 2, pf.rhidev.duration, Fit);
   init_par( 3, 0.1, Fit);
   
   return err;
}
int MicroFinSrcFit::eval(double *date, double *v, double **dVdP, int nobs, int pb ){
 
   int t, n;
   int npb = lcs.count();
   
   double u, A, dAdu, dAdustar, dAdlimb, u2, squ4, dt, f, B, w, limba, limbb;
   double err0, err1, err2;
 
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = umin
   // par[1] = t0
   // par[2] = that
   // par[3] = ustar
   // par[4..4+npb] = f 	lensed fraction of passband i-4
   // par[4+npb..2*npb+4] = B	baseline norm of passband i-4-npb
   
   /* u = sqrt(umin^2 + (2(t-t0)/that)^2 */
   /* where that = v/(2 Re) */
  
   if( pb >= lcs.count( ) ) { // unrecognized pb ==> 1's ??
      for( t = 0; t < nobs; ++t )
         v[t] = 1.0;
      return 0;
   }
   
   if ( (par[0] < 0.0) || (par[3] < 0.0) ) {
      v[0] = Infinity;
      return 0;
   }
 
   if ( pb < 0 || pb >= lcs.count( ) ) {
      limba = limbb = 0;
   } else {
      w = Passband::info(lcs[pb]->filter)->midpoint;
      if ( w <= 3600.0 ) {
         limba = cp.get( "FiniteSourceFit.limb.darkening.a.U", 0.0);
         limbb = cp.get( "FiniteSourceFit.limb.darkening.b.U", 0.0);
      }
      else if ( w > 3600.0 && w <= 4800.0 ) {
         limba = cp.get( "FiniteSourceFit.limb.darkening.a.B", 0.0);
         limbb = cp.get( "FiniteSourceFit.limb.darkening.b.B", 0.0);
      }
      else if ( w > 4800.0 && w <= 6000.0 ) {
         limba = cp.get( "FiniteSourceFit.limb.darkening.a.V", 0.0);
         limbb = cp.get( "FiniteSourceFit.limb.darkening.b.V", 0.0);
      }
      else if ( w > 6000.0 && w <= 8000.0 ) {
         limba = cp.get( "FiniteSourceFit.limb.darkening.a.R", 0.0);
         limbb = cp.get( "FiniteSourceFit.limb.darkening.b.R", 0.0);
      }
      else if ( w > 8000.0 ) {
         limba = cp.get( "FiniteSourceFit.limb.darkening.a.I", 0.0);
         limbb = cp.get( "FiniteSourceFit.limb.darkening.b.I", 0.0);
      }
   }
   
   f = pb >= 0 ? par[4+pb] : 1;        // lensed fraction
   B = pb >= 0 ? par[4+npb+pb] : 1;    // baseline normalization
      

   for( t = 0; t < nobs; ++t ) {                // observations
      u = 2. * (dt = (date[t] - par[1]) ) / par[2];
      u = sqrt(u2 = (par[0] * par[0] + u*u) );
      squ4 = sqrt(u2+4);
      
      v[t] = B * (f * (A = FS_Amp(u, par[3], limba, limbb) ) + 1 - f);
      
      if( dVdP ) {
         dAdu = FS_dAmp(FS_Amp, 0, double(u), double(par[3]), double(limba),
                        double(limbb), 0.001, &err0);
         dAdustar = FS_dAmp(FS_Amp, 1, double(u), double(par[3]),
         double(limba), double(limbb), 0.001, &err1); 
 
/*       dAdu = FS_dAmp(FS_Amp, 0, double(u), double(par[3]), double(limba),
                        double(limbb), double(u), &err0);
         dAdustar = FS_dAmp(FS_Amp, 1, double(u), double(par[3]),
                        double(limba), double(limbb), double(par[3]), &err1); *//*       printf("%.6f %.6f %.6f %.6f %.6f %.6f\n", u, dAdu, err0, par[3], dAdustar, err1); */
 
         dVdP[0][t] = B * f * dAdu * par[0] / u;                      // dFdu
         dVdP[1][t] = -B * f * dAdu * 4 * dt / (u * par[2] * par[2]); // dFdt0
         dVdP[2][t] = dt * dVdP[1][t] /  par[2];                      // dFdthat
         dVdP[3][t] = B * f * dAdustar;                               // dFdustar
	 
         for(n = 0; n < lcs.count( ); ++n) {
	    dVdP[n+4][t] = pb == n ? B*(A-1) : 0;             // Lensed fraction f
	    dVdP[n+npb+4][t] = pb == n ? f*(A-1) + 1 : 0;     // Baseline B
	 }
      }
   }
   return SiteErr::errNone;
}


/*--- MicroParFit Methods ------------------------------------------------*/
 
int MicroParFit::init( LcSet &s ) {
   int p;
   int err = 0;
 
   /*--- resize MicroParFit and point lcs at data ---*/
   if( err = load( s, 5 + 2 * s.count() ) )         // 5 params + fracs & normalizations
      return err;
   npar = 5 + 2 * lcs.count( );                     // some Lc's may have been cut
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( );        	    // points to unlensed baseline

   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
 
 
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Parallax Microlensing" );
   sprintf( pname[0], "U_min" );
   sprintf( pname[1], "T_zero" );
   sprintf( pname[2], "T^hat" );
   sprintf( pname[3], "1/V_t" );
   sprintf( pname[4], "Theta" );
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 7;                    // constraints on each passband
      sprintf( pname[p+5], "%-12s f",
	       Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+5+lcs.count()], "%-12s",
	       Passband::info( (*lcs[p]).filter )->name );
   }
   return 0;
}
 
int MicroParFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;

   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      init_par( p+5, 1.0, Fixed );
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+5+lcs.count(), M.average( ), Fit );
   }

   double A = pf.rhidev.max;
   init_par( 0, sqrt(2) * sqrt(A/sqrt(A*A-1) - 1), Fit);
   init_par( 1, pf.rhidev.mdate, Fit);
   init_par( 2, pf.rhidev.duration, Fit);
   init_par( 3, 0.000001, Fit);
   init_par( 4, -1, Fit);

   return err;
}

int MicroParFit::eval(double *date, double *v, double **dVdP, int nobs, int pb ){
 
   int t, n;
   int npb = lcs.count();
   
   double u, A, dAdu, u2, squ4, dt, f, B;
   double ecc, tperi, el, eb, tc, wo, au, ebrad;
   double w, wt, paralpha, du2dparalpha;
   double du2dumin, dudumin, du2dt0, dudt0, dwdthat, dparalphadw, du2dw, du2dthat, dudthat;
   double dparalphad1vt, du2d1vt, dud1vt, du2dtheta, dudtheta;
   double err0, err1, err2;
 
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = umin
   // par[1] = t0
   // par[2] = that
   // par[3] = 1/v_t
   // par[4] = theta
   // par[5..5+npb] = f 	lensed fraction of passband i-5
   // par[5+npb..2*npb+5] = B	baseline norm of passband i-5-npb
  
   if( pb >= lcs.count( ) ) { // unrecognized pb ==> 1's ??
      for( t = 0; t < nobs; ++t )
         v[t] = 1.0;
      return 0;
   }

   char ecliptic[256], fts[256], *scrp;

   // Problem is, clc title is "MACHO object 1.3450.2557" and we
   //   only want the f.t.s.  Also in fitml.C.
   strcpy(ecliptic, lcs[0]->title);
   scrp = strrchr(ecliptic, ' ');      // ' ' makes \  an int
   *scrp = 0;
   ++scrp;
   strcpy(fts, scrp);
   
   sprintf (ecliptic, "ParallaxFit_%s_ecliptic.l", fts);
   el = cp.get(ecliptic, 270.);                  		 // Use rd2lb.gmn

   sprintf (ecliptic, "ParallaxFit_%s_ecliptic.b", fts);
   eb = cp.get(ecliptic, -5.);                  		 // Use rd2lb.gmn


   ecc = cp.get("ParallaxFit.earth.eccentricity", 0.0167);       // E3 of Astronomical Almanac
   tperi = cp.get("ParallaxFit.machodate.perihelion", 366.313);  // Fri Jan 01 23:30:42 US/Pacific 1993
   tc = cp.get("ParallaxFit.machodate.vernal.equinox", 444.612); // Sun Mar 21 06:41:17 US/Pacific 1993

   wo = 2. * acos(-1.) / 365.256363;  // Sidereal year, C1 of Astronomical Almanac
   au = 149597870.691 / 86400.0;
   ebrad = eb * acos(-1.) / 180.;

   tc += (el - 180.) * (365.256363 / 360.); 

   f = pb >= 0 ? par[5+pb] : 1;        // lensed fraction
   B = pb >= 0 ? par[5+npb+pb] : 1;    // baseline normalization
      
   for( t = 0; t < nobs; ++t ) {                // observations

      dt = (date[t] - par[1]);
      
      w = 2. / par[2];
      wt = wo * (date[t] - tc) + 2. * ecc * sin(wo * (date[t] - tperi));
      paralpha = w * au * par[3] * (1. - ecc * cos(wo * (date[t] - tperi)));

      u2 = par[0]*par[0] + w*w * (date[t] - par[1])*(date[t] - par[1]) + paralpha*paralpha * sin(wt)*sin(wt) +
	 2. * paralpha * sin(wt) * (w * (date[t] - par[1]) * sin(par[4]) + par[0] * cos(par[4])) +
	 paralpha*paralpha * sin(ebrad)*sin(ebrad) * cos(wt)*cos(wt) +
	 2. * paralpha * sin(ebrad) * cos(wt) * (w * (date[t] - par[1]) * cos(par[4]) - par[0] * sin(par[4]));

      u = sqrt(u2);
      v[t] = B * (f * (A = (u2 + 2.) / (u * (squ4 = sqrt(u2+4.) ) ) ) + 1 - f);
      
      if( dVdP ) {
	 dAdu = 2. / squ4 - (u2 + 2.) * (squ4 + u2 / squ4 ) / (u2 * ( u2 + 4. ) );

	 du2dparalpha = 2. * paralpha * sin(wt)*sin(wt) +
	    2. * sin(wt) * (w * (date[t] - par[1]) * sin(par[4]) + par[0] * cos(par[4])) +
	    2. * paralpha * sin(ebrad)*sin(ebrad) * cos(wt)*cos(wt) +
	    2. * sin(ebrad) * cos(wt) * (w * (date[t] - par[1]) * cos(par[4]) - par[0] * sin(par[4]));

	 du2dumin = 2. * par[0] +
	    2. * paralpha * sin(wt) * cos(par[4]) -
	    2. * paralpha * sin(ebrad) * cos(wt) * sin(par[4]);
	 dudumin = du2dumin / 2. / u;

	 du2dt0 = -1. * (2. * w*w * (date[t] - par[1]) +
			 2. * paralpha * sin(wt) * (w * sin(par[4])) +
			 2. * paralpha * sin(ebrad) * cos(wt) * (w * cos(par[4])));
	 dudt0 = du2dt0 / 2. / u;

	 dwdthat = -1. * w / par[2];
	 dparalphadw = paralpha / w;
	 du2dw = 2. * w * (date[t] - par[1])*(date[t] - par[1]) +
	    2. * paralpha * sin(wt) * ((date[t] - par[1]) * sin(par[4])) +
	    2. * paralpha * sin(ebrad) * cos(wt) * ((date[t] - par[1]) * cos(par[4])) +
	    du2dparalpha * dparalphadw;
	 du2dthat = du2dw * dwdthat;
	 dudthat = du2dthat / 2. / u;

	 dparalphad1vt = paralpha / par[3];
	 du2d1vt = du2dparalpha * dparalphad1vt;
	 dud1vt = du2d1vt / 2. / u;

	 du2dtheta = 2. * paralpha * sin(wt) * (w * (date[t] - par[1]) * cos(par[4]) - par[0] * sin(par[4])) +
	    2. * paralpha * sin(ebrad) * cos(wt) * (-1. * w * (date[t] - par[1]) * sin(par[4]) - par[0] * cos(par[4]));
	 dudtheta = du2dtheta / 2. / u;


         dVdP[0][t] = B * f * dAdu * dudumin;                        // dFdumin
         dVdP[1][t] = B * f * dAdu * dudt0;                          // dFdt0
         dVdP[2][t] = B * f * dAdu * dudthat;                        // dFdthat
         dVdP[3][t] = B * f * dAdu * dud1vt;                         // dFd1vt
         dVdP[4][t] = B * f * dAdu * dudtheta;                       // dFdtheta
	 
         for(n = 0; n < npb; ++n) {
	    dVdP[n+5][t] = pb == n ? B*(A-1) : 0;          // Lensed fraction f
	    dVdP[n+npb+5][t] = pb == n ? f*(A-1) + 1 : 0;  // Baseline B
	 }
      }
   }
   return SiteErr::errNone;
}


/*--- MicroBinSrcFit2 Methods ------------------------------------------------*/

/* Form of the equations taken from Dominik, A & A 329, p361.  I set
   eccentricity parameter epsilon and inclination parameter gamma equal
   to zero */
   
int MicroBinSrcFit2::init( LcSet &s ) {
   int p;
   int err = 0;
 
   /*--- resize MicroBinSrcFit2 and point lcs at data ---*/
   if( err = load( s, 10 + 2 * s.count() ) )         // 10 params + fracs & normalizations
      return err;
   npar = 10 + 2 * lcs.count( );                     // some Lc's may have been cut
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( );        	    // points to unlensed baseline

   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
 
 
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Binary Source Microlensing" );
   sprintf( pname[0], "U_min" );
   sprintf( pname[1], "T_zero" );
   sprintf( pname[2], "T^hat" );
   sprintf( pname[3], "Angle alpha" );	// Between source axis and lens velocity
   sprintf( pname[4], "Lum frac" );
   sprintf( pname[5], "Mass frac" );
   sprintf( pname[6], "SemiMajor Axis" );
   sprintf( pname[7], "Orbital Period" );
   sprintf( pname[8], "Orbital Phase at T_0" );
   sprintf( pname[9], "Orbital Inclination" );	// Beta in Dominik
   
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 12;                    // constraints on each passband
      sprintf( pname[p+10], "%-12s f",
	       Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+10+lcs.count()], "%-12s",
	       Passband::info( (*lcs[p]).filter )->name );
   }
   return 0;
}
 
int MicroBinSrcFit2::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;

   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      init_par( p+10, 1.0, Fixed );
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+10+lcs.count(), M.average( ), Fit );
   }

   double A = pf.rhidev.max;
   init_par( 0, sqrt(2) * sqrt(A/sqrt(A*A-1) - 1), Fit);
   init_par( 1, pf.rhidev.mdate, Fit);
   init_par( 2, pf.rhidev.duration, Fit);
   init_par( 3, 0, Fit);
   init_par( 4, 0.5, Fit);
   init_par( 5, 0.5, Fit);
   init_par( 6, 1, Fit);
   init_par( 7, pf.rhidev.duration/5, Fit);
   init_par( 8, 0, Fit);
   init_par( 9, 0, Fit);
   
   return err;
}

int MicroBinSrcFit2::eval(double *date, double *v, double **dVdP, int nobs, int pb ){
 
   int t, n;
   int npb = lcs.count();
   
   double A, A1, A2, u12, u22, u1, u2, dt, f, B;
   double squ14, squ24, zeta, ysrc1, ysrc2, ylens1, ylens2;
   double dAdu1, dAdu2, du12dylens1, du12dylens2, du22dylens1, du22dylens2;
   double du12dysrc1, du12dysrc2, du22dysrc1, du22dysrc2;
   double dylens1dumin, dylens2dumin, du12dumin, du22dumin, du1dumin, du2dumin;
   double dylens1dt0, dylens2dt0, dzetadt0, dysrc1dzeta, dysrc2dzeta;
   double du12dt0, du22dt0, du1dt0, du2dt0;
   double dylens1dthat, dylens2dthat, du12dthat, du22dthat, du1dthat, du2dthat;
   double dylens1dbinalpha, dylens2dbinalpha, du12dbinalpha, du22dbinalpha, du1dbinalpha, du2dbinalpha;
   double du12dmf, du22dmf, du1dmf, du2dmf;
   double dysrc1dsma, dysrc2dsma, du12dsma, du22dsma, du1dsma, du2dsma;
   double dzetadper, du12dper, du22dper, du1dper, du2dper;
   double dzetadphase, du12dphase, du22dphase, du1dphase, du2dphase;
   double dysrc2dinc, du12dinc, du22dinc, du1dinc, du2dinc; 
   
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = U_min
   // par[1] = T_zero
   // par[2] = T^hat
   // par[3] = Angle binalpha
   // par[4] = Lum frac
   // par[5] = Mass frac
   // par[6] = SemiMajor Axis
   // par[7] = Orbital Period
   // par[8] = Orbital Phase at T_0
   // par[9] = Orbital Inclination
   // par[10..10+npb] = f 	   lensed fraction of passband i-10
   // par[10+npb..2*npb+10] = B	   baseline norm of passband i-10-npb
  
   if( pb >= lcs.count( ) ) { // unrecognized pb ==> 1's ??
      for( t = 0; t < nobs; ++t )
         v[t] = 1.0;
      return 0;
   }
   
   f = pb >= 0 ? par[10+pb] : 1;        // lensed fraction
   B = pb >= 0 ? par[10+npb+pb] : 1;	// baseline normalization
      
   for( t = 0; t < nobs; ++t ) {                // observations

      dt = (date[t] - par[1]);

      ylens1 = cos(par[3]) * 2.*(date[t] - par[1])/par[2] - sin(par[3]) * par[0];
      ylens2 = sin(par[3]) * 2.*(date[t] - par[1])/par[2] + cos(par[3]) * par[0];

      zeta = 2. * M_PI * (date[t] - par[1]) / par[7] + par[8];
      
      ysrc1 = par[6] * cos(zeta);
      ysrc2 = par[6] * cos(par[9]) * sin(zeta);

      u12 = ylens1*ylens1 + ysrc1*ysrc1 * (1-par[5])*(1-par[5]) +
	 ylens2*ylens2 + ysrc2*ysrc2 * (1-par[5])*(1-par[5]) +
	 -2. * (1-par[5]) * (ylens1*ysrc1 + ylens2*ysrc2);
      u1 = sqrt(u12);

      u22 = ylens1*ylens1 + ysrc1*ysrc1 * (-par[5])*(-par[5]) +
	 ylens2*ylens2 + ysrc2*ysrc2 * (-par[5])*(-par[5]) +
	 -2. * (-par[5]) * (ylens1*ysrc1 + ylens2*ysrc2);
      u2 = sqrt(u22);


      v[t] = B * f * par[4] * (A1 = (u12 + 2.) / (u1 * (squ14 = sqrt(u12 + 4.) ) ) ) +
	     B * f * (1 - par[4]) * (A2 = (u22 + 2.) / (u2 * (squ24 = sqrt(u22 + 4.) ) ) ) +
	     B * (1-f);
      
      
      if( dVdP ) {
	 
	 dAdu1 = 2. / squ14 - (u12 + 2.) * (squ14 + u12 / squ14 ) / (u12 * ( u12 + 4. ) );
	 dAdu2 = 2. / squ24 - (u22 + 2.) * (squ24 + u22 / squ24 ) / (u22 * ( u22 + 4. ) );

	 du12dylens1 = 2. * ylens1 - 2. * (1-par[5]) * ysrc1;
	 du12dylens2 = 2. * ylens2 - 2. * (1-par[5]) * ysrc2;
	 du22dylens1 = 2. * ylens1 - 2. * (-par[5]) * ysrc1;
	 du22dylens2 = 2. * ylens2 - 2. * (-par[5]) * ysrc2;

	 du12dysrc1 = 2. * ysrc1 * (1-par[5])*(1-par[5]) - 2. * (1-par[5]) * ylens1;
	 du12dysrc2 = 2. * ysrc2 * (-par[5])*(-par[5]) - 2. * (1-par[5]) * ylens2;
	 du22dysrc1 = 2. * ysrc1 * (1-par[5])*(1-par[5]) - 2. * (-par[5]) * ylens1;
	 du22dysrc2 = 2. * ysrc2 * (-par[5])*(-par[5]) - 2. * (-par[5]) * ylens2;

	 dylens1dumin = -sin(par[3]);
	 dylens2dumin = cos(par[3]);
	 du12dumin = du12dylens1 * dylens1dumin + du12dylens2 * dylens2dumin;
	 du22dumin = du22dylens1 * dylens1dumin + du22dylens2 * dylens2dumin;
	 du1dumin = du12dumin / 2. / u1;
	 du2dumin = du22dumin / 2. / u2;

	 dylens1dt0 = -2. * cos(par[3]) / par[2];
	 dylens2dt0 = -2. * sin(par[3]) / par[2];
	 dzetadt0 = -2. * M_PI / par[7];
	 dysrc1dzeta = -par[6] * sin(zeta);
	 dysrc2dzeta = par[6] * cos(par[9]) * cos(zeta);
	 du12dt0 = du12dylens1 * dylens1dt0 + du12dylens2 * dylens2dt0 +
	    du12dysrc1 * dysrc1dzeta * dzetadt0 +
	    du12dysrc2 * dysrc2dzeta * dzetadt0;
	 du22dt0 = du22dylens1 * dylens1dt0 + du22dylens2 * dylens2dt0 +
	    du22dysrc1 * dysrc1dzeta * dzetadt0 +
	    du22dysrc2 * dysrc2dzeta * dzetadt0;
	 du1dt0 = du12dt0 / 2. / u1;
	 du2dt0 = du22dt0 / 2. / u2;

	 dylens1dthat = -cos(par[3]) * 2.*(date[t] - par[1]) / (par[2]*par[2]);
	 dylens2dthat = -sin(par[3]) * 2.*(date[t] - par[1]) / (par[2]*par[2]);
	 du12dthat = du12dylens1 * dylens1dthat + du12dylens2 * dylens2dthat;
	 du22dthat = du22dylens1 * dylens1dthat + du22dylens2 * dylens2dthat;
	 du1dthat = du12dthat / 2. / u1;
	 du2dthat = du22dthat / 2. / u2;

	 dylens1dbinalpha = -sin(par[3]) * 2.*(date[t] - par[1])/par[2] - cos(par[3]) * par[0];
	 dylens2dbinalpha =  cos(par[3]) * 2.*(date[t] - par[1])/par[2] - sin(par[3]) * par[0];
	 du12dbinalpha = du12dylens1 * dylens1dbinalpha + du12dylens2 * dylens2dbinalpha;
	 du22dbinalpha = du22dylens1 * dylens1dbinalpha + du22dylens2 * dylens2dbinalpha;
	 du1dbinalpha = du12dbinalpha / 2. / u1;
	 du2dbinalpha = du22dbinalpha / 2. / u2;

	 du12dmf = (ysrc1*ysrc1 + ysrc2*ysrc2) * -2. * (1-par[5]) + 2. * (ylens1*ysrc1 + ylens2*ysrc2);
	 du22dmf = (ysrc1*ysrc1 + ysrc2*ysrc2) * 2. * par[5] + 2. * (ylens1*ysrc1 + ylens2*ysrc2);
	 du1dmf = du12dmf / 2. / u1;
	 du2dmf = du22dmf / 2. / u2;

	 dysrc1dsma = cos(zeta);
	 dysrc2dsma = cos(par[9]) * sin(zeta);
	 du12dsma = du12dysrc1 * dysrc1dsma + du12dysrc2 * dysrc2dsma;
	 du22dsma = du22dysrc1 * dysrc1dsma + du22dysrc2 * dysrc2dsma;
	 du1dsma = du12dsma / 2. / u1;
	 du2dsma = du22dsma / 2. / u2;

	 dzetadper = -2. * M_PI * (date[t] - par[1]) / (par[7]*par[7]);
	 du12dper = du12dysrc1 * dysrc1dzeta * dzetadper +
	    du12dysrc2 * dysrc2dzeta * dzetadper;
	 du22dper = du22dysrc1 * dysrc1dzeta * dzetadper +
	    du22dysrc2 * dysrc2dzeta * dzetadper;
	 du1dper = du12dper / 2. / u1;
	 du2dper = du22dper / 2. / u2;

	 dzetadphase = 1.;
	 du12dphase = du12dysrc1 * dysrc1dzeta * dzetadphase +
	    du12dysrc2 * dysrc2dzeta * dzetadphase;
	 du22dphase = du22dysrc1 * dysrc1dzeta * dzetadphase +
	    du22dysrc2 * dysrc2dzeta * dzetadphase;
	 du1dphase = du12dphase / 2. / u1;
	 du2dphase = du22dphase / 2. / u2;

	 dysrc2dinc = -par[6] * sin(par[9]) * sin(zeta);
	 du12dinc = du12dysrc2 * dysrc2dinc;
	 du22dinc = du22dysrc2 * dysrc2dinc;
	 du1dinc = du12dinc / 2. / u1;
	 du2dinc = du22dinc / 2. / u2;
	 
         dVdP[0][t] = B * f * ( dAdu1 * du1dumin  + dAdu2 * du2dumin);    // dFdumin
         dVdP[1][t] = B * f * ( dAdu1 * du1dt0    + dAdu2 * du2dt0);      // dFdt0
         dVdP[2][t] = B * f * ( dAdu1 * du1dthat  + dAdu2 * du2dthat);    // dFdthat
         dVdP[3][t] = B * f * ( dAdu1 * du1dbinalpha + dAdu2 * du2dbinalpha);   // dFdbinalpha
         dVdP[4][t] = B * f * ( A1 - A2 );		  		  // dFdlfrac
         dVdP[5][t] = B * f * ( dAdu1 * du1dmf    + dAdu2 * du2dmf);      // dFdmfrac
         dVdP[6][t] = B * f * ( dAdu1 * du1dsma   + dAdu2 * du2dsma);     // dFdsmaxis
	 dVdP[7][t] = B * f * ( dAdu1 * du1dper   + dAdu2 * du2dper);     // dFdperiod
	 dVdP[8][t] = B * f * ( dAdu1 * du1dphase + dAdu2 * du2dphase);   // dFdphase
	 dVdP[9][t] = B * f * ( dAdu1 * du1dinc   + dAdu2 * du2dinc);     // dFdinclination
	 
         for(n = 0; n < npb; ++n) {
	    dVdP[n+10][t] = pb == n ?     B * ( par[4] * A1 + (1 - par[4]) * A2 - 1) : 0;      // f
	    dVdP[n+npb+10][t] = pb == n ? f * ( par[4] * A1 + (1 - par[4]) * A2 - 1) + 1 : 0;  // B
	 }
      }
   }
   return SiteErr::errNone;
}


/*--- MicroYukFit Methods ------------------------------------------------*/
 
int MicroYukFit::init( LcSet &s ) {
   int p;
   int err = 0;
 
   /*--- resize MicroYukFit and point lcs at data ---*/
   if( err = load( s, 5 + 2 * s.count() ) )         // 5 params + normalizations
      return err;
   npar = 5 + 2 * lcs.count( );                     // some Lc's may have been cut
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( );        	    // points to unlensed baseline

   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
 
 
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "PPN (Yukawa) Microlensing" );
   sprintf( pname[0], "U_min" );
   sprintf( pname[1], "T_zero" );
   sprintf( pname[2], "T^hat" );
   sprintf( pname[3], "Alpha" );
   sprintf( pname[4], "Lambda" );
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 7;                    // constraints on each passband
      sprintf( pname[p+5], "%-12s f",
	       Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+5+lcs.count()], "%-12s",
	       Passband::info( (*lcs[p]).filter )->name );
   }
   return 0;
}
 
int MicroYukFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;

   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      init_par( p+5, 1.0, Fixed );
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+5+lcs.count(), M.average( ), Fit );
   }

   double A = pf.rhidev.max;
   init_par( 0, sqrt(2) * sqrt(A/sqrt(A*A-1) - 1), Fit);
   init_par( 1, pf.rhidev.mdate, Fit);
   init_par( 2, pf.rhidev.duration, Fit);
   init_par( 3, 0.01, Fit);
   init_par( 4, 1, Fit);

   return err;
}

int MicroYukFit::eval(double *date, double *v, double **dVdP, int nobs, int pb ){
 
   int t, n;
   int npb = lcs.count();

   double u, A, dAdu, u2, squ4, dt, f, B;
   double dA, ddAdu, ddAda, ddAdl;
   
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = umin
   // par[1] = t0
   // par[2] = that
   // par[3] = alpha
   // par[4] = lambda
   // par[5..5+npb] = f 	lensed fraction of passband i-5
   // par[5+npb..2*npb+5] = B	baseline norm of passband i-5-npb
  
   if( pb >= lcs.count( ) )
      for( t = 0; t < nobs; ++t )
	 v[t] = 1.0;
    
   f = pb >= 0 ? par[5+pb] : 1;        // lensed fraction
   B = pb >= 0 ? par[5+npb+pb] : 1;    // baseline normalization

   for( t = 0; t < nobs; ++t ) { // observations
      u = 2. * (dt = (date[t] - par[1]) ) / par[2];
      u = sqrt(u2 = (par[0] * par[0] + u*u) );

      A = (u2 + 2.) / (u * (squ4 = sqrt(u2+4.) ) );

      dA = ( 1. / ( sqrt(4. + u2) * par[4] ) ) *
	 exp( -1. * sqrt(u2 + 4.) / (2. * par[4]) ) * par[3] *
	 ( sinh(u / (2. * par[4])) - cosh(u / (2. * par[4])) *
	   (2. + u2) / (u * (4. + u2)) );

				  
      
      v[t] = B * (f * (A + dA) + 1 - f);

      if( dVdP ) {

	 dAdu = 2. / squ4 - (u2 + 2.) * (squ4 + u2 / squ4 ) / (u2 * ( u2 + 4. ) );

	 ddAdu = ( 1. / ( u2 * (4. + u2)*(4. + u2) * par[4]*par[4] ) ) *
	    exp( -1 * sqrt(4. + u2) / (2. * par[4]) ) * par[3] *
	    ( sinh(u/(2. * par[4])) * -1 * u * ( 4. + u2*u2 + u2 * (5. + sqrt(4. + u2) * par[4])) +
	      cosh(u/(2. * par[4])) *
	      (8 * par[4] + u2*u2 * (sqrt(4 + u2) + par[4]) + u2 * (3. * sqrt(4. + u2) + 2*par[4])) );
	    
	 ddAdl = ( 1. / ( u * (4. + u2) * par[4]*par[4]*par[4] ) ) *
	    exp( -1 * sqrt(4. + u2) / (2. * par[4]) ) * par[3] *
	    ( sinh( u /(2. * par[4]) ) * u * ( 3. + u2 - sqrt(4. + u2) * par[4] ) -
	      cosh( u /(2. * par[4]) ) *
	      ( sqrt(4. + u2) + u2 * ( sqrt(4. + u2) - par[4] ) - 2. * par[4] ) );

	 ddAda = ( 1. / ( u * (4. + u2) * par[4] ) ) *
	    exp( -1. * sqrt(4. + u2) / (2. * par[4]) ) *
	    ( sinh( u / (2. * par[4]) ) * u * sqrt(4. + u2) -
	      cosh( u / (2. * par[4]) ) * (2. + u2) );



	 dVdP[0][t] = B * f * (dAdu + ddAdu) * par[0] / u; 				// dFdumin
	 dVdP[1][t] = -B * f * (dAdu + ddAdu)  * 4 * dt / (u * par[2] * par[2]);	// dFdt0
	 dVdP[2][t] = dt * dVdP[1][t] /  par[2];					// dFdthat
	 dVdP[3][t] = B * f * ddAda;							// dFdalpha
	 dVdP[4][t] = B * f * ddAdl; 							// dFdlambda
	 
	 for(n = 0; n < npb; ++n) {
	    dVdP[n+5][t] = pb == n ? B*(A+dA-1) : 0;          	  // Lensed fraction f
	    dVdP[n+npb+5][t] = pb == n ? f*(A+dA-1) + 1 : 0;         // Baseline B
	 }
      }
   }
   return SiteErr::errNone;

}



/*--- OUTDATED FITS ---------------------------------------------------------*/
/*--- GaussianFit Methods ---------------------------------------------------*/

int GaussianFit::init( LcSet &s ) {
   int p;
   int err = 0;
   
   /*--- resize GaussianFit and point lcs at data ---*/
   if( err = load( s, 3 + s.count() ) ) // 3 params + normalizations
      return err;
   npar = 3 + lcs.count( );	// some Lc's may have been cut
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( ); // baseline parameters
   
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 200 );
   nchi = cp.get( "GeneralLcFit.nchi", 10) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
   
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Gaussian + 1" );
   sprintf( pname[0], "Amplitude" );
   sprintf( pname[1], "Center" );
   sprintf( pname[2], "Width" );
   for( p = 0; p < lcs.count( ); ++p ) {
      sprintf( pname[p+3], "%s", Passband::info( (*lcs[p]).filter )->name );
      cnstrt[p] = 4;	// constraints on passband
   }
   for( p = 0; p < npar; ++p ) 
      init_par( p, NaN, Guess );

   return err;
}

int GaussianFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;
   
   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+3, M.average( ), Fit );
   }
   init_par( 0, pf.rhidev.max - 1, Fit );
   init_par( 1, pf.rhidev.mdate, Fit );
   init_par( 2, pf.rhidev.duration, Fit );
   return err;
}

int GaussianFit::eval( double *date, double *v, double **dVdP, int nobs, int pb){
   
   double B, arg, ex, fac;
   int n, t;
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = amplitude
   // par[1] = t0
   // par[2] = sigma
   // par[3..npar] = B_i 	(baseline normalization of passband n-3)
   
   // v = B_i ( Amp * exp( -((t-t0)/sigma)^2 ) + 1 )
   
   B = pb >= 0 ? par[3+pb] : 1; // Baseline normalization
   for( t = 0; t < nobs; ++t ) { // observations
      arg = (date[t]-par[1])/par[2];
      ex = exp(-arg*arg);
      v[t] = B * ( par[0] * ex + 1 );
      fac = B * par[0] * ex * 2 * arg;
      if( dVdP ) {
	 dVdP[0][t] = B * ex;
	 dVdP[1][t] = fac / par[2];
	 dVdP[2][t] = fac * arg / par[2];
	 for(n = 0; n < lcs.count( ); ++n)	
	    dVdP[n+3][t] = pb == n ? par[0] * ex + 1: 0;
      }
   }
   return SiteErr::errNone;
}


/*--- MicroBlendFit Methods -------------------------------------------------*/

int MicroBlendFit::init( LcSet &s ) {
   int p;
   int err = 0;
   
   /*--- resize MicroBlendFit and point lcs at data ---*/
   if( err = load( s, 3 + 2 * s.count() ) )	// 3 params + fracs & norms
      return err;
   npar = 3 + 2 * lcs.count( ); // some Lc's may have been cut
   
   /*--- set norm if you which peak stats to be normalized ---*/
   norm = par + npar - lcs.count( ); // baseline parameters 
   
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
   
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Blended Microlensing" );
   sprintf( pname[0], "U_min" );
   sprintf( pname[1], "T_zero" );
   sprintf( pname[2], "T^hat" );
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 5;	// constraints on each passband
      sprintf( pname[p+3], "%s f",
	      Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+3+lcs.count()], "%s",
	      Passband::info( (*lcs[p]).filter )->name );
   }
   
   return err;
}

int MicroBlendFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;
   
   /*--- guess parameter values ---*/
   for( p = 0; p < lcs.count( ); ++p ) {
      init_par( p+3, 1.0, Fit );
      M.eval( *lcs[p], 0.1, 0.7 );
      init_par( p+3+lcs.count(), M.average( ), Fit );
   }
   double A = pf.rhidev.max;
   init_par( 0, sqrt(2) * sqrt(A/sqrt(A*A-1) - 1), Fit);
   init_par( 1, pf.rhidev.mdate, Fit);
   init_par( 2, pf.rhidev.duration, Fit);

   return err;
}

int MicroBlendFit::eval( double *date, double *v, double **dVdP,
			int nobs, int pb ) {
   
   int t, n;
   double u, A, dAdu, u2, squ4, dt, B, f;
   int npb = lcs.count();
   
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = umin
   // par[1] = t0
   // par[2] = that
   // par[3..3+npb) = f_i 	(lensed fraction from passband i-3)
   // par[3+npb..2npb+3) = B_i	(baseline norm of passband i-3-npb)
   
   /* u = sqrt(umin^2 + (2(t-t0)/that)^2 */
   /* where that = v/(2 Re) */
   
   if( pb >= lcs.count( ) )
      for( t = 0; t < nobs; ++t )
	 v[t] = 1.0;
   
   
   f = pb >= 0 ? par[3+pb] : 1;         // lensed fraction
   B = pb >= 0 ? par[3+npb+pb] : 1;	// Baseline normalization
   for( t = 0; t < nobs; ++t ) {        // observations
      u = 2 * (dt = (date[t] - par[1]) ) / par[2];
      u = sqrt(u2 = (par[0] * par[0] + u*u) );
      v[t] = B*(f*(A=(u2+2)/(u*(squ4=sqrt(u2+4)))) + 1-f);
      if( dVdP ) {
	 dAdu = 2/squ4-(u2+2)*(squ4+u2/squ4)/(u2*(u2+4));
	 dVdP[0][t] = B*f*dAdu*par[0]/u;                // dFdu
	 dVdP[1][t] = -B*f*dAdu*4*dt/(u*par[2]*par[2]); // dFdt0
	 dVdP[2][t] = dt*dVdP[1][t]/par[2];	        // dFdthat
	 for(n = 0; n < npb; ++n) {
	    dVdP[n+3][t] = pb == n ? B*(A-1) : 0;       // Lensed fraction f
	    dVdP[n+npb+3][t] = pb == n ? f*(A-1)+1 : 0; // Baseline B
	 }
      }
   }
   return SiteErr::errNone;
}


/*--- MicroBinSrcFit Methods ------------------------------------------------*/

int MicroBinSrcFit::init( LcSet &s ) {
   int p;
   int err = 0;
   
   /*--- resize MicroBinSrcFit and point lcs at data ---*/
   if( err = load( s, 6 + 2 * s.count() ) )     // 3 params + fracs & norms
      return err;
   npar = 6 + 2 * lcs.count( ); // some Lc's may have been cut
   
   /*--- set norm if you want stats to be normalized ---*/
   norm = par + npar - lcs.count( ); // baseline parameters 
   
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
   
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Double Source Microlensing (old)" );
   sprintf( pname[0], "U_min(1)" );
   sprintf( pname[1], "T_zero(1)" );
   sprintf( pname[2], "T^hat(1)" );
   sprintf( pname[3], "U_min(2)" );
   sprintf( pname[4], "T_zero(2)" );
   sprintf( pname[5], "T^hat(2)-T^hat(1)" );
   for( p = 0; p < lcs.count( ); ++p ) {
      sprintf( pname[p+6], "%s f(1)",
              Passband::info( (*lcs[p]).filter )->name );
      sprintf( pname[p+6+lcs.count()], "%s",
              Passband::info( (*lcs[p]).filter )->name );
   }
   return err;
}

int MicroBinSrcFit::init( LcSet &s, PeakFinder &pf ) {
   int p;
   int err;
   if( err = init( s ) )
      return err;

   /*--- guess parameter values ---*/
   if( ptype[5] == Guess ) {
      init_par( 5, 0.0, Fixed ); // assume equal t^hats
      nfxpar = 1;        // this is no longer an indep constraint
   }
   for( p = 0; p < lcs.count( ); ++p ) {
      cnstrt[p] = 8 - nfxpar; // constraints on each passband
      if( ptype[p+6] == Guess )
         init_par( p+6, 1.0, Fit );
      if( ptype[p+6+lcs.count()] == Guess ) {
	 M.eval( *lcs[p], 0.1, 0.7 );
	 init_par( p+6+lcs.count(), M.average( ), Fit );
      }
   }
   
   return err;
}
   
int MicroBinSrcFit::eval( double *date, double *v, double **dVdP,
			 int nobs, int pb ) {
   
   int t, n;
   double u_1, A_1, dA_1du, u2_1, squ4_1, dt_1, B, f, th2, tmp;
   double u_2, A_2, dA_2du, u2_2, squ4_2, dt_2;
   int npb = lcs.count();
   
   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = umin_1			// source 1 param
   // par[1] = t0_1			// source 1 param
   // par[2] = that_1			// source 1 param
   // par[3] = umin_2			// source 2 param
   // par[4] = t0_2			// source 2 param
   // par[5] = that_2 - that_1		// source 2 param
   // par[6..6+npb) = f_i 	(source_1 fraction in passband i)
   // par[6+npb..6+2npb) = B_i	(baseline norm of passband i)
   
   /* u = sqrt(umin^2 + (2(t-t0)/that)^2 */
   /* where that = v/(2 Re) */
   
   if( pb >= lcs.count( ) ) { // unrecognized pb ==> 1's ??
      for( t = 0; t < nobs; ++t )
	 v[t] = 1.0;
      return 0;
   }
   
   f = pb >= 0 ? par[6+pb] : 1; // source_1 fraction
   B = pb >= 0 ? par[6+npb+pb] : 1;	// Baseline normalization
   for( t = 0; t < nobs; ++t ) { // observations
      u_1 = 2 * (dt_1 = (date[t] - par[1]) ) / par[2];
      u_1 = sqrt(u2_1 = (par[0] * par[0] + u_1*u_1) );
      
      u_2 = 2 * (dt_2 = (date[t] - par[4]) ) / (th2 = par[5]+par[2]);
      u_2 = sqrt(u2_2 = (par[3] * par[3] + u_2*u_2) );
      
      v[t] = B*f*(A_1=(u2_1+2)/(u_1*(squ4_1=sqrt(u2_1+4))))
	 +B*(1-f)*(A_2=(u2_2+2)/(u_2*(squ4_2=sqrt(u2_2+4))));
      if( dVdP ) {
	 dA_1du = 2/squ4_1-(u2_1+2)*(squ4_1+u2_1/squ4_1)/(u2_1*(u2_1+4));
	 dA_2du = 2/squ4_2-(u2_2+2)*(squ4_2+u2_2/squ4_2)/(u2_2*(u2_2+4));
	 dVdP[0][t] = B*f*dA_1du*par[0]/u_1; // dFdu_1
	 dVdP[1][t] = -B*f*dA_1du*4*dt_1/(u_1*par[2]*par[2]); // dFdt0_1
	 dVdP[2][t] = dt_1*dVdP[1][t]/par[2];	      // dFdthat_1
	 dVdP[3][t] = B*f*dA_2du*par[3]/u_2;    // dFdu_2
	 dVdP[4][t] = -B*f*dA_2du*4*dt_2/(u_1*th2*th2); // dFdt0_2
	 tmp = dt_2*dVdP[4][t]/th2;		// dFdthat_2
	 dVdP[5][t] = 1.0/(1.0/tmp - 1.0/dVdP[2][t]); // dF/ddthat
	 for(n = 0; n < npb; ++n) {
	    dVdP[n+6][t] = pb == n ? B*(A_1-A_2) : 0; // dA/df
	    dVdP[n+npb+6][t] = pb == n ? f*A_1+(1-f)*A_2 : 0; // Amp
	 }
      }
   }
   return SiteErr::errNone;
}

/*--- SineFit Methods ---------------------------------------------------*/

int SineFit::init( LcSet &s ) {
   int p;
   int err = 0;
   
   /*--- resize SineFit and point lcs at data ---*/
   if( err = load( s, 3 ) ) // 3 params
      return err;
   npar = 3;

   /*--- set norm if you wish peak stats to be normalized ---*/
   norm = 0;
   
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
   
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Sine" );
   sprintf( pname[0], "Amplitude" );
   sprintf( pname[1], "Period" );
   sprintf( pname[2], "Phase" );
   for( p = 0; p < lcs.count( ); ++p ) 
      cnstrt[p] = 3;	// constraints on each passband

   for( p = 0; p < npar; ++p ) {
      par[p] = NaN;
      ptype[p] = Guess;					// auto guess
   }
   return err;
}

int SineFit::eval( double *date, double *v, double **dVdP, int nobs, int pb) {
   
   double arg, a, s, c, p, f;
   int t;

   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = amplitude
   // par[1] = period
   // par[2] = phase
   
   // v = amplitude *  sin( 2_PI * (t - phase)  / period )
   a = par[0];
   f = 1.0/par[1];
   p = par[2];
   for( t = 0; t < nobs; ++t ) { // observations
      v[t] = a * (s = sin( arg = 2.0 * M_PI * (date[t]-p) * f ) );
      if( dVdP ) {
	 dVdP[0][t] = s;
	 dVdP[1][t] = - arg * f * a * ( c = cos( arg ) );
	 dVdP[2][t] = - 2 * M_PI * f * a * c;
      }
   }
   return SiteErr::errNone;
}

/*--- CosineFit Methods ---------------------------------------------------*/

int CosineFit::init( LcSet &s ) {
   int p;
   int err = 0;
   
   /*--- resize CosineFit and point lcs at data ---*/
   if( err = load( s, 3 ) ) // 3 params
      return err;
   npar = 3;

   /*--- set norm if you wish peak stats to be normalized ---*/
   norm = 0;
   
   /*--- set clif parameters ---*/
   maxdchi2 = cp.get( "GeneralLcFit.maxdchi2", double(0.001) );
   maxfchi2 = cp.get( "GeneralLcFit.maxfchi2", double(1e-6) );
   maxit = cp.get( "GeneralLcFit.maxit", 100 );
   nchi = cp.get( "GeneralLcFit.nchi", 6) ;
   pkmin = cp.get( "GeneralLcFit.peak.min", double(1.05) );
   
   /*--- overwrite generic parameter names ---*/
   sprintf( name, "Cosine" );
   sprintf( pname[0], "Amplitude" );
   sprintf( pname[1], "Period" );
   sprintf( pname[2], "Phase" );
   for( p = 0; p < lcs.count( ); ++p ) 
      cnstrt[p] = 3;	// constraints on each passband

   for( p = 0; p < npar; ++p ) {
      par[p] = NaN;
      ptype[p] = Guess;					// auto guess
   }
   return err;
}

int CosineFit::eval( double *date, double *v, double **dVdP, int nobs, int pb) {
   
   double arg, a, s, c, p, f;
   int t;

   /*--- Fill in v[] and dvdp[][] if internal Lc ---*/
   // par[0] = amplitude
   // par[1] = period
   // par[2] = phase
   
   // v = amplitude *  cos( 2_PI * (t - phase)  / period )
   a = par[0];
   f = 1.0 / par[1];
   p = par[2];
   for( t = 0; t < nobs; ++t ) { // observations
      arg = 2.0 * M_PI * f * (date[t]-p);
      v[t] = a * (c = cos( arg ) );
      if( dVdP ) {
	 dVdP[0][t] = c;
	 dVdP[1][t] = arg * f * a * ( s = sin( arg ) );
	 dVdP[2][t] = 2 * M_PI * f * a * s;
      }
   }
   return SiteErr::errNone;
}

/*-------------- PeakFinder Methods ---------------------------------------*/
PeakFinder::PeakFinder( ClifParameters &cp ) {
   fsingle = f0 = f1 = 0;
   maxn = 0;
   penalty = cp.get( "PeakFinder.point.penalty", double(4.0) );
   nfilters = ( int ) cp.get( "PeakFinder.num.filters", double(10.0) );
}

PeakFinder::~PeakFinder( void ) {
   delete [] fsingle;
   delete [] f0;
   delete [] f1;
}

int PeakFinder::init( int npts ) {
   if( npts > maxn ) {
      delete [] fsingle;
      delete [] f0;
      delete [] f1;
      if( !(fsingle = new element[npts]) )
	 return SiteErr::errMemory;
      if( !(f0 = new element[npts]) )
	 return SiteErr::errMemory;
      if( !(f1 = new element[npts]) )
	 return SiteErr::errMemory;
   }
   maxn = npts;
   element newpf;
   newpf.date = newpf.mdate = newpf.duration = 0.0;
   newpf.max = - ( newpf.min = Infinity );
   newpf.dev =  newpf.pdev = newpf.ndev = -Infinity;
   newpf.hipoint = newpf.lopoint = -Infinity;
   newpf.Fn = newpf.npts = newpf.np = 0;

   hidev = newpf;
   hi2dev = newpf;
   rhidev = newpf;
   rhi2dev = newpf;
   lodev = newpf;
   lo2dev = newpf;
   rlodev = newpf;
   rlo2dev = newpf;

   return SiteErr::errNone;
}   

int PeakFinder::eval( LcSet &lcs, double *N ) {

   /*--- currently maximum number of passbands is 64 ---*/
   double *date[64];			// current position in each date vector
   double *dlim[64];			// end of date vectors
   double norm[64];			// normalization for each value vector
   int next[64];			// index of next pb's
   int nnext, P, i, n;
   double d;
   double lastdate = NaN, lodate;
   double rlod, rhid;
   int hipeak, rhipeak, lopeak, rlopeak;

   /*--- resize the filter and find norms ---*/
   int Nobs = 0;
   int npb = lcs.count( );		// number of passbands finished
   int p;
   for( p = 0; p < lcs.count( ); ++p ) {
      Nobs += lcs[p]->nobs;
      if( N ) {
	 norm[p] = N[p];
      } else {
	 M.eval( *lcs[p], 0.1, 0.9 );
	 norm[p] =  M.average( );
      }
      date[p] = lcs[p]->nobs ? lcs[p]->date : 0;
      if( !date[p] )
	 --npb;
      dlim[p] = lcs[p]->date + lcs[p]->nobs;
   }
   if( Nobs <= 2*lcs.count( ) )
      return SiteErr::errMissing;
   init( Nobs );
   element *f, *fp, *fpp, *f1d;
   f = fsingle;			// load single point filter in fsingle
   fpp = f0;			// reset f0;
   fp = f1;			// reset f1;

   /*--- make the single point deviance vector ---*/
   for( n = 0; n < Nobs; ++n, ++f, ++fp, ++fpp ) {
      if( npb == 0 )				// all passbands finished
	 break;
      lodate = Infinity;
      nnext = 0;

      /*--- look through all passbands for next point ---*/
      for( p = 0; p < lcs.count( ); ++p ) {
	 if( !date[p] )
	    continue;
	 else if( *date[p] < lodate ) {
	    lodate = *date[p];
	    nnext = 1;
	    next[0] = p;
	 } else if( *date[p] == lodate ) {
	    next[nnext++] = p;		// keep track of simultaneous points
	 }
      }
      
      if( lodate == Infinity || nnext == 0 ) 	// out of points 
	 break;					// this shouldn't happen

      /*--- single time slice deviance calculation ---*/
      f->dev = f->pdev = f->ndev = - nnext * penalty;
      f->date = lodate;
      f->duration = n == 0 ? 0.0 : lodate - lastdate; 	// important definition
      f->npts = nnext;			// number of indep points in bin
      f->np = 0;			// number of positive points in bin
      f->max = 0.0;
      for( p = 0; p < nnext; ++p ) {	// loop through simultaneous points
	 P = next[p];			// next passband
	 i = date[P] - lcs[P]->date;	// current index in passband P
	 if( ++date[P] >= dlim[P] ) {	// increment date vector P
	    date[P] = 0;
	    --npb;
	 }
	 f->max += lcs[P]->units == CommonLc::Linear ?
	    lcs[P]->value[i] / norm[P] :	// Amplification
	    lcs[P]->value[i] - norm[P];		// dMag
	 d = ( lcs[P]->value[i] - norm[P] ) / lcs[P]->error[i];		// chi
	 if( d > 0 ) {			// positive excursion
	    d = d*d;			// chi2
	    f->pdev += d;
	    f->np++;
	 } else {			// negative excursion
	    d = d*d;			// chi2
	    f->ndev += d;
	 }
	 f->dev += d;
      }
      /*--- extrema are treated as averages in filter #1 ---*/
      f->min = f->max /= nnext;		// avg of all simultaneous points
      f->hipoint = f->pdev;
      f->lopoint = f->ndev;
      //printf("%9.3f %9.4f %9.3f %9.3f %3d %3d  p: %8.3f  n: %8.3f\n",
      //     f->date, f->duration, f->max, f->min,
      //     f->npts, f->np, f->pdev, f->ndev );
      lastdate = lodate;
   }
   npoints = n;			// use for loop index to set

   /*--- begin Fibonacci recursion for extrema ---*/
   rlod = -Infinity;		// intermediate result, robust lo dev
   rhid = -Infinity;		// intermediate result, robust hi dev 
   fsize = pfsize = 1;		// set first and second filter sizes
   int F, Fn;
   for( F = 2; F < nfilters; ++F ) {

      /*--- set filter pointers to correct buffers ---*/
      f = ( F % 2 ? f1 : f0 ) + npoints - 1;
      fp = ( F % 2 ? f0 : f1 ) + npoints - 1;
      fpp = f - fsize;	
      f1d = fsingle + npoints - fsize;		// one pt duration fix
      Fn = pfsize + fsize;			// Fibonacci_n
      int minn = Fn - 1;			// first non-zero element

      /*--- handle special cases --*/
      if( F == 2 ) {
	 fp = fsingle + npoints - 1;
	 fpp = fp - 1;
      } else if ( F == 3 ) {
	 fpp = fsingle + npoints - 3;
      }

      //printf("Fibonacci[%d]: %d points\n", F, minn + 1 );

      /*--- loop backward so as not to clobber fpp early --*/
      for( n = npoints - 1; n >= minn; --n ) {

	 /*--- calculate composite filter element ---*/
	 f->Fn = Fn;
	 f->date = fsingle[n].date;
	 f->mdate = fsingle[n - Fn/2].date;  // date of midpoint
	 /*--- note the f1d single point duration correction ---*/
	 f->duration = F == 2 ? fp->duration : F == 3 ?
	    fp->duration + f1d->duration :
	    fp->duration + fpp->duration + f1d->duration;
	 f->max = fp->max > fpp->max ? fp->max : fpp->max;
	 f->min = fp->min < fpp->min ? fp->min : fpp->min;
	 f->dev = fp->dev + fpp->dev;
	 f->pdev = fp->pdev + fpp->pdev;
	 f->ndev = fp->ndev + fpp->ndev;
	 f->hipoint = fp->hipoint > fpp->hipoint ? fp->hipoint : fpp->hipoint;
	 f->lopoint = fp->lopoint > fpp->lopoint ? fp->lopoint : fpp->lopoint;
	 f->npts = fp->npts + fpp->npts;
	 f->np = fp->np + fpp->np;

	 if( f->pdev > hidev.pdev ) {
	    hidev = *f;
	    hipeak = n;
	 }
	 if( F > 2 && (f->pdev - f->hipoint) > rhid ) {
	    rhidev = *f;
	    rhid = rhidev.pdev - rhidev.hipoint;
	    rhipeak = n;
	 }
	 if( f->ndev > lodev.ndev ) {
	    lodev = *f;
	    lopeak = n;
	 }
	 if( F > 2 && (f->ndev - f->lopoint) > rlod ) {
	    rlodev = *f;
	    rlod = rlodev.ndev - rlodev.lopoint;
	    rlopeak = n;
	 }
	 /*--- step back and continue ---*/
	 --f;
	 --fp;
	 --fpp;
	 --f1d;
      }
      pfsize = fsize;
      fsize = minn + 1;
   }

   /*--- begin Fibonacci recursion for 2nd most extreme bins ---*/
   rlod = -Infinity;		// intermediate result
   rhid = -Infinity;		// intermediate result

   /*--- exclude prev peak dates in search ---*/
   int hthi = hipeak + 1;
   int htlo = hipeak - hidev.Fn - 1;
   int rhthi = rhipeak + 1;
   int rhtlo = rhipeak - rhidev.Fn - 1;
   int lthi = lopeak + 1;
   int ltlo = lopeak - lodev.Fn - 1;
   int rlthi = rlopeak + 1;
   int rltlo = rlopeak - rlodev.Fn - 1;
   int nlo, nhi;

   fsize = pfsize = 1;			// set first and second filter sizes

   for( F = 2; F < nfilters; ++F ) {

      /*--- set filter pointers to correct buffers ---*/
      f = ( F % 2 ? f1 : f0 ) + npoints - 1;
      fp = ( F % 2 ? f0 : f1 ) + npoints - 1;
      fpp = f - fsize;	
      f1d = fsingle + npoints - fsize;		// one pt duration fix
      Fn = pfsize + fsize;			// Fibonacci_n
      int minn = Fn - 1;			// lowest non-zero element

      /*--- handle special cases --*/
      if( F == 2 ) {
	 fp = fsingle + npoints - 1;
	 fpp = fp - 1;
      } else if ( F == 3 ) {
	 fpp = fsingle + npoints - 3;
      }

      /*--- loop backward so as not to overwrite fpp early --*/
      for( n = npoints - 1; n >= minn; --n ) {

	 /*--- calculate composite filter element ---*/
	 f->Fn = Fn;
	 f->date = fsingle[n].date;
	 f->mdate = fsingle[n - Fn/2].date;  // date of midpoint
	 /*--- note the f1d single point duration correction ---*/
	 f->duration = F == 2 ? fp->duration : F == 3 ?
	    fp->duration + f1d->duration :
	    fp->duration + fpp->duration + f1d->duration;
	 f->max = fp->max > fpp->max ? fp->max : fpp->max;
	 f->min = fp->min < fpp->min ? fp->min : fpp->min;
	 f->dev = fp->dev + fpp->dev;
	 f->pdev = fp->pdev + fpp->pdev;
	 f->ndev = fp->ndev + fpp->ndev;
	 f->hipoint = fp->hipoint > fpp->hipoint ? fp->hipoint : fpp->hipoint;
	 f->lopoint = fp->lopoint > fpp->lopoint ? fp->lopoint : fpp->lopoint;
	 f->npts = fp->npts + fpp->npts;
	 f->np = fp->np + fpp->np;
	 nlo = n - f->Fn - 1;
	 nhi = n + 1;

	 /*--- keep track of 2nd most extreme elements ---*/
	 if( f->pdev > hi2dev.pdev && (nhi < htlo || nlo > hthi ) )
	    hi2dev = *f;
	 if( F > 2 && ( f->pdev - f->hipoint) > rhid &&
	     ( nhi < rhtlo || nlo > rhthi ) ) {
	    rhi2dev = *f;
	    rhid = rhi2dev.pdev - rhi2dev.hipoint;
	 }
	 if( f->ndev > lo2dev.ndev && ( nhi < ltlo || nlo > lthi ) )
	    lo2dev = *f;
	 if( F > 2 && (f->ndev - f->lopoint) > rlod &&
	     ( nhi < rltlo || nlo > rlthi ) ) {
	    rlo2dev = *f;
	    rlod = rlo2dev.ndev - rlo2dev.lopoint;
	 }
	 /*--- step back and continue ---*/
	 --f;
	 --fp;
	 --fpp;
	 --f1d;
      }
      pfsize = fsize;
      fsize = minn + 1;
   }
}

int PeakFinder::print( FILE *fp ) {
   fprintf( fp, "PeakFinder: %d filters on %d points\n", nfilters, npoints );
   fprintf(fp,"                 value     peak       date     length   obs"
	   "   pts contrib\n");
   fprintf( fp,"  hidev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	   hidev.pdev, hidev.max, hidev.mdate, hidev.duration,
	   hidev.Fn, hidev.npts, hidev.np );
   fprintf( fp," hi2dev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	    hi2dev.pdev, hi2dev.max, hi2dev.date, hi2dev.duration,
	    hi2dev.Fn, hi2dev.npts, hi2dev.np );
   fprintf( fp," rhidev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	  (rhidev.pdev - rhidev.hipoint - 1),
	    rhidev.max, rhidev.date, rhidev.duration,
	    rhidev.Fn, rhidev.npts, rhidev.np - 1);
   fprintf( fp,"rhi2dev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	  (rhi2dev.pdev - rhi2dev.hipoint),
	    rhi2dev.max, rhi2dev.date, rhi2dev.duration,
	    rhi2dev.Fn, rhi2dev.npts, rhi2dev.np - 1);
   fprintf( fp,"  lodev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	    lodev.ndev, lodev.min, lodev.date, lodev.duration,
	    lodev.Fn, lodev.npts, lodev.npts - lodev.np );
   fprintf( fp," lo2dev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	    lo2dev.ndev, lo2dev.min, lo2dev.date, lo2dev.duration,
	    lo2dev.Fn, lo2dev.npts, lo2dev.npts - lo2dev.np );
   fprintf( fp," rlodev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	  (rlodev.ndev - rlodev.lopoint),
	    rlodev.min, rlodev.date, rlodev.duration,
	    rlodev.Fn, rlodev.npts, rlodev.npts - rlodev.np - 1 );
   fprintf( fp,"rlo2dev = %12.4f %8.3f @ %8.3f x %8.3f %5d %5d %7d\n",
	  (rlo2dev.ndev - rlo2dev.lopoint),
	    rlo2dev.min, rlo2dev.date, rlo2dev.duration,
	    rlo2dev.Fn, rlo2dev.npts, rlo2dev.npts - rlo2dev.np - 1 );
   return 0;
}

