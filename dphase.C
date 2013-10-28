#include <string.h>
#include <math.h>
#include <libclif.H>
#include "phase.H"

#define NSSTAB 		10000

/*--- FORTRAN prototypes ---*/
extern "C" {
    int	supdbl_(int*,double*
                ,double*,double*,int*,double*,double*,double*,double*);
    int	dqrdc_(double*,int*,int*,int*,double*,double*,double*,int*);
    int	dqrsl_(double*,int*,int*,int*,double*
               ,double*,double*,double*,double*,double*,double*,int*,int*);
    int	percur_(int*,int*,double*,double*,double*,int*,double*
                ,int*,int*,double*,double*,double*,double*,int*,int*,int*);
    int	splev_(double*,int*,double*,int*,double*,double*,int*,int*);
}

/*--- SSTab methods ------------------------------------------------------*/

SSTab::SSTab( void ) {
    ok = maxnobs = ndata = 0;
    pindex = dindex = nexti = 0;
    F = 0;
    dv = xv = yv = wv = fv = dfv = work = 0;
    BASS = clifpar.get( "ReimannPhaser.BASS", 0.0 );
}

SSTab::~SSTab( void ) {
    delete [] dv;
    delete [] xv;
    delete [] yv;
    delete [] wv;
    delete [] fv;
    delete [] dfv;
    delete [] work;
    delete [] pindex;
    delete [] dindex;
    delete [] nexti;
}

int SSTab::reset( int npts ) {
    ok = 0;
    if( !nexti ) 
        if( !(nexti = new int[NSSTAB]) )
            return SiteErr::errMemory;
    
    if( npts <= maxnobs ) {
        ndata = npts;
        return 0;
    }
    delete [] dv;
    delete [] xv;
    delete [] yv;
    delete [] wv;
    delete [] fv;
    delete [] dfv;
    delete [] work;
    delete [] pindex;
    delete [] dindex;
    if( !(dv = new double[npts]) )
        return SiteErr::errMemory;
    if( !(xv = new double[npts]) )
        return SiteErr::errMemory;
    if( !(yv = new double[npts]) )
        return SiteErr::errMemory;
    if( !(wv = new double[npts]) )
        return SiteErr::errMemory;
    if( !(fv = new double[npts]) )
        return SiteErr::errMemory;
    if( !(dfv = new double[npts]) )
        return SiteErr::errMemory;
    if( !(work = new double[npts*7]) )
        return SiteErr::errMemory;
    if( !(pindex = new int[npts]) )
        return SiteErr::errMemory;
    if( !(dindex = new int[npts]) )
        return SiteErr::errMemory;
    maxnobs = ndata = npts;
    return 0;
}

/*------------------------------------------------------------------------*
  SSTab::init folds lc by period and calculates a supersmoothed curve
  and its derivative w.r.t. phase at the sampled phases.  To facilitate
  evaluation at arbitrary phases, an index lookup table is created to
  quickly locate the bracketing sampled phases.
  *------------------------------------------------------------------------*/
int SSTab::init( const CommonLc &lc, double period ) {
    
    if( lc.nobs <= 0 ) 
        return SiteErr::errMissing;
    
    if( reset( lc.nobs ) )
        return SiteErr::errMemory;
    
    int i, j, k, two=2;
    double span = 0.0;
    
    F = 1.0/period;
    double *pp = (double *) work;	// 
    for( i = 0; i < ndata; ++i ) {
        pp[i] = fmod( lc.date[i] * F, 1.0);		// phase [0,1]
        if( pp[i] < 0.0 || pp[i] >= 1.0 ) {
            printf("fmod out of range %10.4f %10.4f %10.5f\n", F, lc.date[i], pp[i] ); 
        }
    }
    isort( pindex, pp, ndata );
    
    /*--- load double precision buffers ---*/
    for( i = 0; i < ndata; ++i ) {
        j = pindex[i];
        dindex[j] = i;
        xv[i] = pp[j];
        dv[i] = lc.date[j];
        yv[i] = lc.value[j];
        wv[i] = 1.0 / lc.error[j];
    }
    
    /*--- fortran supersmoother call ---*/
    supdbl_(&ndata, xv, yv, wv, &two, &span, &BASS, fv, work);
    
    /*--- calculate derivative of supersmoothed curve ---*/
    for( i = 1; i < ndata-1; ++i ) 
        dfv[i] = (fv[i+1]-fv[i-1])/(xv[i+1]-xv[i-1]);
    dfv[0] = (fv[1]-fv[ndata-1])/(xv[1]-xv[ndata-1]+1.0);	//first
    dfv[ndata-1] = (fv[0]-fv[ndata-2])/(xv[0]-xv[ndata-2]+1.0);	//last
    
    /*--- create index -> phase map ---*/
    for( i = 0, j = 0; i < ndata; ++i ) {
        k = int( NSSTAB * xv[i] );
        for( ; j < k; ++j )
            nexti[j] = i;
    }
    for( ; j < NSSTAB; ++j )
        nexti[j] = 0;
    
    if( 0 ) for( i = 0; i < ndata; ++i )
        printf("ss: %3d %10.4f %10.5f %5d  %8.3f %8.3f %8.3f\n",
               i, dv[i], xv[i], int( NSSTAB * xv[i] ),
               yv[i], fv[i], dfv[i] ); 
    ok = 1;
    return 0;
}

int SSTab::init( const CommonLc &lc, double period, double p0) {
    
    if( lc.nobs <= 0 ) 
        return SiteErr::errMissing;
    
    if( reset( lc.nobs ) )
        return SiteErr::errMemory;
    
    int i, j, k, two=2;
    double span = 0.0;
    
    F = 1.0/period;
    double *pp = (double *) work;	// 
    for( i = 0; i < ndata; ++i ) {
        pp[i] = fmod( (lc.date[i]-p0) * F, 1.0);		// phase [0,1]
        if( pp[i] < 0.0 || pp[i] >= 1.0 ) {
            printf("fmod out of range %10.4f %10.4f %10.5f\n", F, lc.date[i], pp[i] ); 
        }
    }
    isort( pindex, pp, ndata );
    
    /*--- load double precision buffers ---*/
    for( i = 0; i < ndata; ++i ) {
        j = pindex[i];
        dindex[j] = i;
        xv[i] = pp[j];
        dv[i] = lc.date[j];
        yv[i] = lc.value[j];
        wv[i] = 1.0 / lc.error[j];
    }
    
    /*--- fortran supersmoother call ---*/
    supdbl_(&ndata, xv, yv, wv, &two, &span, &BASS, fv, work);
    
    /*--- calculate derivative of supersmoothed curve ---*/
    for( i = 1; i < ndata-1; ++i ) 
        dfv[i] = (fv[i+1]-fv[i-1])/(xv[i+1]-xv[i-1]);
    dfv[0] = (fv[1]-fv[ndata-1])/(xv[1]-xv[ndata-1]+1.0);	//first
    dfv[ndata-1] = (fv[0]-fv[ndata-2])/(xv[0]-xv[ndata-2]+1.0);	//last
    
    /*--- create index -> phase map ---*/
    for( i = 0, j = 0; i < ndata; ++i ) {
        k = int( NSSTAB * xv[i] );
        for( ; j < k; ++j )
            nexti[j] = i;
    }
    for( ; j < NSSTAB; ++j )
        nexti[j] = 0;
    
    if( 0 ) for( i = 0; i < ndata; ++i )
        printf("ss: %3d %10.4f %10.5f %5d  %8.3f %8.3f %8.3f\n",
               i, dv[i], xv[i], int( NSSTAB * xv[i] ),
               yv[i], fv[i], dfv[i] ); 
    ok = 1;
    return 0;
}

/*------------------------------------------------------------------------*
  SSTab::ss evaluates the supersmoother curve at date X.  This is an
  interpolation of the  neighboring sampled phases.
  *------------------------------------------------------------------------*/
double SSTab::ss( double X ) {		// interpolate ss curve
    if( !ok  )
        return NaN;
    static int I, N;
    static double p;
    
    p = fmod( X * F, 1.0 );
    N = int( NSSTAB * p );     		// map phase into index space
    if( N >= NSSTAB ) {			// should not happen
        printf("Warning index out of bounds %12g %d\n",p,N);
        N = 0;
        p -= 1.0;
    }
    I = nexti[N];
    if( I ) 					
        return (fv[I] * (p-xv[I-1]) + fv[I-1] * (xv[I]-p)) / (xv[I]-xv[I-1]);
    else		// wraparound
        return (fv[0] * (p-xv[ndata-1]+1.0) + fv[ndata-1] * (xv[0]-p))/
            (xv[0]-xv[ndata-1]+1.0);
}

/*------------------------------------------------------------------------*
  SSTab::eval evaluates the supersmoother curve and its derivative at
  date X.  This is an interpolation of the neighboring sampled phases.
  *------------------------------------------------------------------------*/
int SSTab::eval( double X, double *f, double *df ) {
    if( !ok  )
        return 1;
    static int I, N;
    static double p;
    
    p = fmod( X * F, 1.0 );
    N = int( NSSTAB * p );     		// map phase into index space
    if( N >= NSSTAB ) {			// should not happen
        printf("Warning index out of bounds %12g %d\n",p,N);
        N = 0;
        p -= 1.0/F;				// -epsilon
    }
    I = nexti[N];
    if( I ) {
        *f = (fv[I]*(p-xv[I-1]) + fv[I-1]*(xv[I]-p)) / (xv[I]-xv[I-1]);
        *df = (dfv[I] * (p-xv[I-1]) + dfv[I-1] * (xv[I]-p)) / (xv[I]-xv[I-1]);
    } else {		// wraparound
        *f = (fv[0] * (p-xv[ndata-1]+1.0) + fv[ndata-1] * (xv[0]-p))/
            (xv[0]-xv[ndata-1]+1.0);
        *df = (dfv[0]*(p-xv[ndata-1]+1.0) + dfv[ndata-1]*(xv[0]-p))/
            (xv[0]-xv[ndata-1]+1.0);
    }
    return 0;
}

/*------------------------------------------------------------------------*
  SSTab::fom (figure of merit) is average absolute value in sigmas of
  deviations from the supersmoothed curve.
  *------------------------------------------------------------------------*/
double SSTab::fom( void ) {
    if( !ok )
        return NaN;
    int i;
    double efn = 0;
    for( i = 0; i < ndata; ++i ) 
        efn += fabs((yv[i]-fv[i])*wv[i]);
    return efn/ndata;
}

/*------------------------------------------------------------------------*
  SSTab::chi2 per d.o.f of the supersmoother fit.
  *------------------------------------------------------------------------*/
double SSTab::chi2( void ) {
    if( !ok )
        return NaN;
    int i;
    double efn = 0;
    for( i = 0; i < ndata; ++i ) 
        efn += pow((yv[i]-fv[i])*wv[i], 2);
    return efn/ndata;
}


/*--- DPhase methods ------------------------------------------------------*/

DPhase::DPhase( ClifParameters &CP ) {
    obsperfit = (int) CP.get( "DPhase.Min.Obs.Bin", double(20.0) );
    maxfrac =  CP.get( "DPhase.Min.Obs.Bin", double(0.1) );
    if( (int) CP.get( "DPhase.Units.AU", double(1.0) ) )
        tunits = 86400 * 2.997924e10 / 1.495979e13;
    else
        tunits = 1.0;
}

DPhase::~DPhase( void ) {
}

int DPhase::derive( const CommonLc &dLc, double period ) {
    
    /*--- SSTab::init loads xv,dv,yf,wv,fv,dfv,nexti,dindex,pindex  ---*/
    init( dLc, period );
    if( !ok ) {
        return SiteErr::errProgram;
    }
    
    pLc.resize( dLc.nobs );
    ssLc.resize( dLc.nobs );
    dpLc.resize( dLc.nobs );
    
    lc_hdr_cpy( dLc, pLc );
    lc_hdr_cpy( dLc, ssLc );
    lc_hdr_cpy( dLc, dpLc );
    
    ssLc.phtpkg = CommonLc::Theory;
    dpLc.units = CommonLc::Linear;
    
    int i, j;
    
    /*--- fill in easy data ---*/
    for( i = 0; i < dLc.nobs; ++i ) {
        pLc.date[i] = ssLc.date[i] = xv[i];		// phase {0,1} units
        pLc.value[i] = yv[i];				// data values
        ssLc.value[i] = fv[i];				// smoothed values
        ssLc.error[i] = 0;
        j = pindex[i];
        pLc.obsid[i] = ssLc.obsid[i] = dLc.obsid[j];
        pLc.error[i] = dLc.error[j];
    }
    
    
    if( tunits != 1.0 )
        tunits = 86400 * 2.997924e10 / 1.495979e13;		// days -> AU
    else
        tunits = 1.0;
    
    /*----------------------------------------------------------------------*
      Calculate dphase values on groups of dLc points.  Groups should
      be less than MAXDPBIN of the timespan but contain more than
      MINDPOBS points.  Care must be taken to keep these value reasonable.
      *----------------------------------------------------------------------*/
    
    int nlast = -1;
    int ilo = 0, ihi = 0;
    int ndp = 0;
    double maxtbin =  maxfrac * ( dLc.date[dLc.nobs-1] - dLc.date[0] );
    while( ihi < dLc.nobs - obsperfit ) {
        ilo = nlast + 1;
        ihi = nlast + obsperfit;
        if( 0 && dLc.date[ihi] - dLc.date[ilo] > maxtbin ) {
            printf("dp bin too long, skipping obs %10.4f\n", dLc.date[ilo] );
            ++nlast;
            continue;
        }
        /*--- date[ilo] -- date[ihi] is an adequate bin ---*/
        
        if( !fitphase( ilo, ihi, ndp ) )
            ++ndp;
        nlast = ihi;
    }
    dpLc.nobs = ndp;
    return 0;
}

/*------------------------------------------------------------------------*
  fitphase fits the supersmoother curve to the region [ilo,ihi] of the
  data lightcurve (dLc).  The fit is in one parameter, delta date.
  This has an "expectation" value of ~0 for a periodic lightcurve.
  The result are written into the dphase lightcurve (dpLc) at location
  ndp.  The fit is performed with an adaptation of the
  Levenberg-Marquardt chi^2 minimization algorithm.
  *------------------------------------------------------------------------*/
int DPhase::fitphase( int ilo, int ihi, int ndp ) {
    
    if( !ok || ilo < 0 || ihi >= ndata )
        return 1;
    double alpha, beta;			// current guess
    double xalpha, xbeta, xchi2;		// last good guess
    double chi2, dchi2;
    double lambda;
    double dok, dtest;
    int niter;
    xchi2 = dchi2 = Infinity;
    lambda = 0.001;
    dtest = 0.0;
    niter = 0;
    while( dchi2 > 0.0001 && niter < 100) {
        alphabet( &alpha, &beta, &chi2, dtest, ilo, ihi );
        if( xchi2 < chi2 ) {		// reguess with smaller step 
            lambda *= 10.0;
            dtest = dok + xbeta / (xalpha * (1.0 + lambda) );
        } else {				// take bigger step next time 
            dok = dtest;			// update best guess
            dchi2 = xchi2 - chi2;		// dchi2 > 0 for improvement
            lambda /= 10.0;
            xalpha = alpha;		// remember current chi2 surface
            xbeta = beta;			// remember current chi2 surface
            xchi2 = chi2;
            dtest = dok + beta / (alpha * (1.0 + lambda) );
            ++niter;			// only successful steps count 
        }
    }
    dpLc.value[ndp] = fmod(dok, 1.0) * tunits;
    dpLc.error[ndp] = sqrt(1.0/alpha) * tunits;
    dpLc.date[ndp] = 0;
    int i;
    for( i = ilo; i <= ihi; ++i )
        dpLc.date[ndp] += dv[dindex[i]];
    dpLc.date[ndp] /= (ihi - ilo + 1);
    return 0;
}

/*------------------------------------------------------------------------*
  alphabet calculates alpha, beta and chi^2 of the region dLc[ilo,ihi]
  to the supersmoother fit offset by delta date = dd.
  alpha = Sum( (df/dd)^2 / sigma^2  )
  beta = Sum( (data - f) * (df/dd) / sigma^2  )
  f and df/dd are supersmoother interpolations at date + dd
  *------------------------------------------------------------------------*/
int DPhase::alphabet( double *alpha, double *beta, double *chi2, double dd,
                      int ilo, int ihi ) {
    static int i, j;
    static double ie2, f, df, dy;
    *alpha = *beta = *chi2 = 0;
    for( i = ilo; i <= ihi; ++i ) {
        j = dindex[i];
        ie2 = wv[j] * wv[j];
        eval( dv[j] + dd, &f, &df );
        *alpha += ie2 * df * df;
        *beta += ie2 * df * (dy = (yv[j] - f) );
        *chi2 += ie2 * dy * dy;
    }
    return 0;
}


/*--- FitPhase methods ------------------------------------------------------*/

FitPhase::FitPhase( ) {
    period = 0;
    dperiod = 0;
}

FitPhase::~FitPhase( void ) {
}

int FitPhase::fitphase(CommonLc &jLc0, CommonLc &hLc0, CommonLc &kLc0) {
    double alpha, beta;			// current guess
    double xalpha, xbeta, xchi2;		// last good guess
    double chi2, dchi2;
    double lambda;
    double dok, dtest;
    //float dlambda = 10;
    float dlambda = 2;
    int niter;
    int nitermax=1000;
    
    jLc = jLc0;
    hLc = hLc0;
    kLc = kLc0;
    ndata = jLc.nobs + hLc.nobs + kLc.nobs;
    reset(ndata);
    
    xchi2 = dchi2 = Infinity;
    lambda = 0.001;
    dtest = period;
    niter = 0;
    while( (dchi2 > 1e-7) && (niter < nitermax) ) {
        alphabet( &alpha, &beta, &chi2, dtest );
        
        if( xchi2 < chi2 ) {		// reguess with smaller step 
            lambda *= dlambda;
            dtest = dok + xbeta / (xalpha * (1.0 + lambda) );
        } else {				// take bigger step next time 
            dok = dtest;			// update best guess
            dchi2 = xchi2 - chi2;		// dchi2 > 0 for improvement
            lambda /= dlambda;
            xalpha = alpha;		// remember current chi2 surface
            xbeta = beta;			// remember current chi2 surface
            xchi2 = chi2;
            dtest = dok + beta / (alpha * (1.0 + lambda) );
            ++niter;			// only successful steps count 
        }
        //fprintf(stderr, "Iteration %d : p = %f, chi^2 = %f, dchi^2 = %f\n", niter, dtest, chi2, dchi2);
        
        
    }
    period  = dok;
    dperiod = sqrt(1.0/alpha);
    
    //fprintf(stderr, "Final Fit Period : %.7f +/- %.7f; chi^2 = %.3f; npts = %d\n", period, dperiod, chi2, ndata);
    
    return 0;
}

int FitPhase::fitphase(CommonLc &Lc0) {
    double alpha, beta;			// current guess
    double xalpha, xbeta, xchi2;		// last good guess
    double chi2, dchi2;
    double lambda;
    double dok, dtest;
    //float dlambda = 10;
    float dlambda = 2;
    int niter;
    int nitermax=1000;
    
    /* ACB - THIS IS A HORRIBLE HACK */
    jLc = Lc0;
    hLc = Lc0;
    kLc = Lc0;
    ndata = jLc.nobs + hLc.nobs + kLc.nobs;
    reset(ndata);
    
    xchi2 = dchi2 = Infinity;
    lambda = 0.001;
    dtest = period;
    niter = 0;
    while( (dchi2 > 1e-7) && (niter < nitermax) ) {
        alphabet( &alpha, &beta, &chi2, dtest );
        
        if( xchi2 < chi2 ) {		// reguess with smaller step 
            lambda *= dlambda;
            dtest = dok + xbeta / (xalpha * (1.0 + lambda) );
        } else {				// take bigger step next time 
            dok = dtest;			// update best guess
            dchi2 = xchi2 - chi2;		// dchi2 > 0 for improvement
            lambda /= dlambda;
            xalpha = alpha;		// remember current chi2 surface
            xbeta = beta;			// remember current chi2 surface
            xchi2 = chi2;
            dtest = dok + beta / (alpha * (1.0 + lambda) );
            ++niter;			// only successful steps count 
        }
        //fprintf(stderr, "Iteration %d : p = %f, chi^2 = %f, dchi^2 = %f\n", niter, dtest, chi2, dchi2);
        
        
    }
    period  = dok;
    dperiod = sqrt(1.0/alpha);
    
    //fprintf(stderr, "Final Fit Period : %.7f +/- %.7f; chi^2 = %.3f\n", period, dperiod, chi2);
    
    return 0;
}

int FitPhase::alphabet( double *alpha, double *beta, double *chi2, double dp) {
    static int i, j;
    static double ie2, f, df, dy;
    *alpha = *beta = *chi2 = 0;
    
    jSs.init(jLc, dp);
    hSs.init(hLc, dp);
    kSs.init(kLc, dp);
    
    //fprintf(stderr, "Caw %d %d %d %f\n", jSs.ndata, hSs.ndata, kSs.ndata, dp);
    
    for( i = 0; i < jSs.ndata; ++i ) {
        j = jSs.dindex[i];
        ie2 = jSs.wv[j] * jSs.wv[j];
        jSs.eval( jSs.dv[j], &f, &df );
        *alpha += ie2 * df * df;
        *beta += ie2 * df * (dy = (jSs.yv[j] - f) );
        *chi2 += ie2 * dy * dy;
    }
    
    for( i = 0; i < hSs.ndata; ++i ) {
        j = hSs.dindex[i];
        ie2 = hSs.wv[j] * hSs.wv[j];
        hSs.eval( hSs.dv[j], &f, &df );
        *alpha += ie2 * df * df;
        *beta += ie2 * df * (dy = (hSs.yv[j] - f) );
        *chi2 += ie2 * dy * dy;
    }
    
    for( i = 0; i < kSs.ndata; ++i ) {
        j = kSs.dindex[i];
        ie2 = kSs.wv[j] * kSs.wv[j];
        kSs.eval( kSs.dv[j], &f, &df );
        *alpha += ie2 * df * df;
        *beta += ie2 * df * (dy = (kSs.yv[j] - f) );
        *chi2 += ie2 * dy * dy;
    }
    
    return 0;
}


