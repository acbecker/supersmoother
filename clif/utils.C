/* written by John Doug Reynolds, March 1998 */

#include "libclif.H"
#include <libsite.H>
#include <string.h>
#include <math.h>

char*	duplicate( const char *str )
{
   return str ? strcpy(new char [strlen(str)+1],str) : 0;
}

/* --- isort ---------------------------------------------------------------

   Adapted from the _Numerical Recipes in C_ implementation of Quicksort.

   Writes the list of indices representing the input array as sorted in
   ascending order.  The order of the input array is not changed.  In the
   unlikely event that the internal stack is too small, an error message
   is printed on stderr and the return value is 1; otherwise the return
   value is 0.

   This implementation does not appear to be optimized for arrays of
   more than 10000 elements.  It may be possible to improve performance
   for much larger arrays by sorting both the indices and a duplicate
   array of values.

   A generalization along the lines of qsort() would be a nice.

   John Doug Reynolds, June 1996

*/

inline void swap( int &i, int &j )
{
   int k;
   k = i; i = j; j = k;
}

int isort( int idx[], const double val[], int len )
{
   const int stacksize = 50;
   int stack[stacksize], top = 0;
   int p0 = 0, p1 = len - 1;
   int h, i, j, k;
   double a;

   for ( i = 0; i < len; ++i ) idx[i] = i;

   while ( 1 ) {
      if ( p1 - p0 < 7 ) {
	 for ( j = p0 + 1; j <= p1; ++j ) {
	    a = val[h=idx[j]];
	    for ( i = j-1; i >= 0 && val[idx[i]] > a; --i ) idx[i+1] = idx[i];
	    idx[i+1] = h;
	 }
	 if ( top == 0 ) break;
	 p1 = stack[--top];
	 p0 = stack[--top];
      }
      else {
	 k = (p0 + p1) >> 1;
	 swap( idx[k], idx[p0+1] );
	 if ( val[idx[p0+1]] > val[idx[p1]] ) swap( idx[p0+1], idx[p1] );
	 if ( val[idx[p0]]   > val[idx[p1]] ) swap( idx[p0],   idx[p1] );
	 if ( val[idx[p0+1]] > val[idx[p0]] ) swap( idx[p0+1], idx[p0] );
	 i = p0 + 1;
	 j = p1;
	 a = val[h=idx[p0]];
	 while ( 1 ) {
	    while ( val[idx[++i]] < a );
	    while ( val[idx[--j]] > a );
	    if ( j < i ) break;
	    swap( idx[i], idx[j] );
	 }
	 idx[p0] = idx[j];
	 idx[j] = h;
	 if ( top + 2 > stacksize ) {
	    fprintf( stderr, "ERROR: isort internal stack overflow\n" );
	    return 1;
	 }
	 if ( p1 - i + 1 >= j - p0 ) {
	    stack[top++] = i;
	    stack[top++] = p1;
	    p1 = j - 1;
	 }
	 else {
	    stack[top++] = p0;
	    stack[top++] = j - 1;
	    p0 = i;
	 }
      }
   }

   return 0;
}

int isort( int idx[], const int val[], int len )
{
   const int stacksize = 50;
   int stack[stacksize], top = 0;
   int p0 = 0, p1 = len - 1;
   int h, i, j, k;
   int a;

   for ( i = 0; i < len; ++i ) idx[i] = i;

   while ( 1 ) {
      if ( p1 - p0 < 7 ) {
	 for ( j = p0 + 1; j <= p1; ++j ) {
	    a = val[h=idx[j]];
	    for ( i = j-1; i >= 0 && val[idx[i]] > a; --i ) idx[i+1] = idx[i];
	    idx[i+1] = h;
	 }
	 if ( top == 0 ) break;
	 p1 = stack[--top];
	 p0 = stack[--top];
      }
      else {
	 k = (p0 + p1) >> 1;
	 swap( idx[k], idx[p0+1] );
	 if ( val[idx[p0+1]] > val[idx[p1]] ) swap( idx[p0+1], idx[p1] );
	 if ( val[idx[p0]]   > val[idx[p1]] ) swap( idx[p0],   idx[p1] );
	 if ( val[idx[p0+1]] > val[idx[p0]] ) swap( idx[p0+1], idx[p0] );
	 i = p0 + 1;
	 j = p1;
	 a = val[h=idx[p0]];
	 while ( 1 ) {
	    while ( val[idx[++i]] < a );
	    while ( val[idx[--j]] > a );
	    if ( j < i ) break;
	    swap( idx[i], idx[j] );
	 }
	 idx[p0] = idx[j];
	 idx[j] = h;
	 if ( top + 2 > stacksize ) {
	    fprintf( stderr, "ERROR: isort internal stack overflow\n" );
	    return 1;
	 }
	 if ( p1 - i + 1 >= j - p0 ) {
	    stack[top++] = i;
	    stack[top++] = p1;
	    p1 = j - 1;
	 }
	 else {
	    stack[top++] = p0;
	    stack[top++] = j - 1;
	    p0 = i;
	 }
      }
   }

   return 0;
}

int	fold( class CommonLc &lc1
	      , const CommonLc &lc2, double period, double date0 )
{
   if ( isnan(period) || isnan(date0) || period <= 0 )
      return SiteErr::errParams;

   lc1.copy(lc2,1); // copy header only
   lc1.refJD += date0; // shift origin to date0

   if ( lc2.nobs <= 0 ) return 0;
   if ( lc1.resize(lc2.nobs) ) return SiteErr::errMemory;

   int i, nobs = lc1.nobs;
   int *index = new int [nobs];
   double *phase = new double [nobs];
   double invper = 1.0 / period;
   for ( i = 0; i < nobs; ++i ) {
      double cycle = invper * (lc2.date[i] - date0);
      phase[i] = fmod(cycle,1.0) + (cycle < 0 ? 1 : 0);
   }

   isort( index, phase, nobs );

   for ( i = 0; i < nobs; ++i ) {
      int j = index[i];
      lc1.date[i] = phase[j];
      lc1.value[i] = lc2.value[j];
      lc1.error[i] = lc2.error[j];
      lc1.obsid[i] = lc2.obsid[j];
   }

   delete [] phase;
   delete [] index;
   return 0;
}

/*-------------- LcMerge ----------------------------------------------------*/
/*-------------- written by Mark Pratt  Nov 1995 ----------------------------*/

int LcMerge( const CommonLc &lc1, const CommonLc &lc2, CommonLc &lc,
	     const ClifParameters &par )
{
   int verbose = par.get( "LcMerge.verbose", 0 );

   const CommonLc *lcp;
   char msg[256];
   lc.erase( );
   
   if( lc1.norm0 != lc2.norm0 ||
       lc1.norm1 != lc2.norm1 ||
       lc1.filter != lc2.filter ||
       lc1.units != lc2.units ||
       lc1.refJD != lc2.refJD ) 
      return SiteErr::errFormat;

   int N, sodpair;

   /* It is ok to merge different Sodophot types of the same filter */
   sodpair = lc1.phtpkg == CommonLc::Sodophot ||
      lc1.phtpkg == CommonLc::pSodophot || lc1.phtpkg == CommonLc::oSodophot;
   sodpair &= lc2.phtpkg == CommonLc::Sodophot ||
      lc2.phtpkg == CommonLc::pSodophot || lc2.phtpkg == CommonLc::oSodophot;

   if( (lc1.phtpkg != lc2.phtpkg) && !sodpair )
      return SiteErr::errFormat;	// missmatched photometry packages

   /*--- Check for zero length Lc's ---*/
   if( !lc1.nobs && !lc2.nobs )	{			// don't bother
      return SiteErr::errNone;
   }

   lc.resize( N = lc1.nobs + lc2.nobs );   // potential overallocation, bfd
   /*--- fill in merged header ---*/
   lc.ra = lc1.ra;
   lc.dec = lc1.dec;
   lc.refJD = lc1.refJD;
   lc.norm0 = lc1.norm0;
   lc.norm1 = lc1.norm1;
   lc.units = lc1.units;
   lc.filter = lc1.filter;
   lc.phtpkg = sodpair ? CommonLc::Sodophot : lc1.phtpkg;
   if ( lc1.title ) {
      lc.title = new char [ strlen(lc1.title) + 1 ];
      if ( lc.title ) strcpy( lc.title, lc1.title );
   }
   
   int t, t1, t2, T, s;
   double Sv = 0.0, Svv = 0.0, v;
   double ierr2, S = 0.0;
   int nmat = 0;

   for( T = t1 = t2 = t = 0; t < N; ++t ) {

      if( t1 == lc1.nobs && t2 == lc2.nobs )		// Ooops!
	 break;
      else if( t1 == lc1.nobs )				// no more lc1 obs
	 s = 1;				
      else if( t2 == lc2.nobs )				// no more lc2 obs
	 s = -1;
      else
         s = lc1.obsid[t1] > lc2.obsid[t2] ? 1 :
         lc1.obsid[t1] < lc2.obsid[t2] ? -1 : 0;


      switch (s) {
      case 1 :				// add lc2(t2) to lc
	 lcp = &lc2; T = t2++; break;
      case -1 : 			// add lc1(t1) to lc
	 lcp = &lc1; T = t1++; break;
      case 0 : 				// duplicate observation
	 lcp = &lc1; T = t1++;		// take lc1(t1) values as default
	 ++t2;	
	 // Keep track of stddev in offset as diagnostic
	 if ( verbose ) {
	    ierr2 = 1 / (lc1.error[t1] * lc1.error[t1]);
	    v = lc.units == CommonLc::Magnitude ?
	       (lc1.value[t1] - lc2.value[t2]) : lc1.value[t1] / lc2.value[t2];
	    S += ierr2;
	    Sv += v * ierr2;
	    Svv += v * v * ierr2;
	    ++nmat;
	 }
	 break;
      default :
	 return SiteErr::errProgram;
      }
      lc.date[t] = lcp->date[T];
      lc.value[t] = lcp->value[T];
      lc.error[t] = lcp->error[T];
      lc.obsid[t] = lcp->obsid[T];
   }

   if( nmat ) {
      printf("LcMerge: Average offset: %.4f (%0.4f) over %d measurments\n",
		Sv/S, sqrt(S*Svv-Sv*Sv)/S, nmat );
   }
   lc.nobs = t;
   return SiteErr::errNone;
}

void LcMerge( LcSet &s, const ClifParameters &par ) {
   CommonLc *lc1, *lc2;
   CommonLc *lcp = s.newlc( );
   int merged;
   int n1, n2;
   for( n1 = 0; n1 < s.count( );  ) {
      lc1 = s[n1];
      merged = 0;
      for( n2 = n1 + 1; n2 < s.count( ); ++n2 ) {
	 lc2 = s[n2];
	 if( !LcMerge( *lc1, *lc2, *lcp, par ) ) {
	    merged = 1;
	    s.recycle( lc1 );
	    s.recycle( lc2 );
	    s.append( lcp );
	    lcp = s.newlc( );
	    break;
	 } 
      }
      if( !merged ) ++n1;
   }
   s.append( lcp );
   s.recycle( lcp );
}

/*-------------- LcMoment Methods ------------------------------------------*/
/*-------------- written by Mark Pratt  Sep 1996 ---------------------------*/
LcMoment::LcMoment( void ) {
   maxn = 0;
   sindex = 0;
   reset( 0 );
}
LcMoment::~LcMoment( void ) {
   delete [] sindex;
}
int LcMoment::reset( int npts ) {
   ok = sorted = N = nlo = nhi = 0;
   S = Sv = Sdv2 = Sdv3 = Sdv4 = Se = 0;
   if( npts > maxn ) {
      delete [] sindex;
      if( !( sindex = new int[(int)(1.2 * npts)]) ) {
	 return SiteErr::errMemory;
      }
      maxn = (int)(1.2 * npts);
   }
   return SiteErr::errNone;
}
int LcMoment::eval( CommonLc &lc, double plimlo, double plimhi ) {
   if( lc.nobs <= 0 ) {
      return SiteErr::errMissing;
   }
   
   if( reset( lc.nobs ) )
      return SiteErr::errMemory;

   nlo = int(plimlo * lc.nobs);
   nhi = int(plimhi * lc.nobs);
   nlo = nlo < 0 ? 0 : nlo >= lc.nobs ? lc.nobs-1 : nlo;
   nhi = nhi < 0 ? 0 : nhi >= lc.nobs ? lc.nobs-1 : nhi;
   
   if( nlo > nhi ) {
      return SiteErr::errParams;
   }

   if(nlo > 0 || nhi < lc.nobs-1 ) {
      if( isort( sindex, lc.value, lc.nobs  ) ) 
	 return SiteErr::errMemory;
      sorted = 1;
   }

   int m, n;
   double *v = lc.value;
   double *e = lc.error;
   double ie2, v0, dv, ve2;

   if( sorted ) { 
      for( n = nlo; n <= nhi; ++n ) {
	 m = sindex[n];
	 S += ie2 = 1.0 / (e[m] * e[m]);
	 Sv += v[m] * ie2;
	 Se += e[m];
	 ++N;
      }
      v0 = Sv/S;
      for( n = nlo; n <= nhi; ++n ) {
	 m = sindex[n];
	 ie2 = 1.0 / (e[m] * e[m]);
	 dv = v[m] - v0;
	 Sdv2 += ve2 = dv * dv * ie2;
	 Sdv3 += ve2 * dv;
	 Sdv4 += ve2 * dv * dv;
      }
   } else {
      for( n = 0; n < lc.nobs; ++n, ++v, ++e ) {
	 S += ie2 = 1.0 / (*e * *e);
	 Sv += *v * ie2;
	 Se += *e;
	 ++N;
      }
      v0 = Sv/S;
      v = lc.value;
      e = lc.error;
      for( n = 0; n < lc.nobs; ++n, ++v, ++e ) {
	 ie2 = 1.0 / (*e * *e);
	 dv = *v - v0;
	 Sdv2 += ve2 = dv * dv * ie2;
	 Sdv3 += ve2 * dv;
	 Sdv4 += ve2 * dv * dv;
      }
   }
   ok = 1;
   return SiteErr::errNone;
}
int LcMoment::pindex( double percentile ){
   if( !ok || !sorted || percentile < 0.0 || percentile > 1.0)
      return -1;
   int i = int( percentile * _size);
   return sindex[i >= _size ? _size-1 : i];
}
double LcMoment::moment( int m ){
   if( !ok ) return NaN;
   /*--- error weighted moments ---*/
   switch( m ) {
   case 0 :
      return S / N;
   case 1 :
      return Sv / S;
   case 2 :
      return Sdv2 / S;
   case 3:
      return Sdv3 / S;
   case 4:
      return Sdv4 / S;
   default:
      return NaN;
   }
}
double LcMoment::average( void ){
   if( !ok ) return NaN;
   return (Sv/S);
}
double LcMoment::chisq( void ) {
   if( !ok ) return NaN;
   return Sdv2;
}
double LcMoment::skewness( void ){
   if( !ok ) return NaN;
   return moment(3) / pow(moment(2), 1.5);
}
double LcMoment::kurtosis( void ){
   if( !ok ) return NaN;
   return moment(4) / pow(moment(2),2.0);
}
int LcMoment::print( FILE *fp, int verbose ) {
   if( maxn == 0 ) {		// print header
      fprintf( fp,  "#    N          S         Sv       Sdv2"
	       "       Sdv2       Sdv3       Sdv4         Se\n" );
   } else if( !ok ) {
      fprintf( fp, "#LcMoment::print invalid\n");
	 return SiteErr::errMissing;
   }
   
   if( verbose ) {
      fprintf( fp,"%-10s: %-12d %-10s: %-12g %-10s: %-12g\n",
	       "N", N, "S", S, "Sv", Sv ); 
      fprintf( fp,"%-10s: %-12g %-10s: %-12g %-10s: %-12g\n",
	       "Sdv2", Sdv2, "Sdv3", Sdv3, "Sdv4", Sdv4 );
      if( verbose > 1 ) {
	 fprintf( fp,"%-10s: %-12g %-10s: %-12g %-10s: %-12g\n",
		  "average", average( ), "error", sqrt(1/S),
		  "avg error", Se/N );
	 fprintf( fp,"%-10s: %-12g %-10s: %-12g %-10s: %-12g\n",
		  "chi2/df",  chisq( )/N, "skewness", skewness( ),
		  "kurtosis", kurtosis( ) );
      }
   } else {
      fprintf( fp, "%6d %10g %10g %10g %10g %10g %10g\n",
	       N, S, Sv, Sdv2, Sdv3, Sdv4, Se );
   }
   return 0;
}

