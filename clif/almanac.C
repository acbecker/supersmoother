/* written by Welch and Reynolds, June 1996 */

#include "libclif.H"
#include <math.h>

static double jul_cen( double mday )
{
   // Returns the fraction of a Julian century since J2000.0
   return ( ( mday - 2451545.0 ) / 36525.0 ) ;
}

static double precess_m ( double t )
{
   // Returns the value of the precession constant M 
   // Explanatory Supplement to the Astronomical Almanac
   // page 106 (3.213-2)
   double dr = M_PI / 180.0 ;
   return ( ( 1.2812323 + ( 0.0003879 + 0.0000101 * t ) * t ) * t * dr );
}

static double precess_n ( double t )
{
   // Returns the value of the precession constant N
   // Explanatory Supplement to the Astronomical Almanac
   // page 106 (3.213-2)
   double dr = M_PI / 180.0 ;
   return ( ( 0.5567530 - ( 0.0001185 + 0.0000116 * t ) * t ) * t * dr );
}

static double ra_mid ( double t, double ra0, double dec0 )
{
   // Returns the value of the ra at a time halfway between the
   // equinox of date and equinox J2000.0
   return ( ra0 + 0.5*( precess_m(t) + precess_n(t)*sin(ra0)*tan(dec0)) );
}

static double dec_mid ( double t, double ra0, double dec0 )
{
   // Returns the value of the dec at a time halfway between the
   // equinox of date and equinox J2000.0
   return ( dec0 + 0.5*precess_n(t)*ra_mid(t, ra0, dec0) );
}

static double ra_date ( double t, double ra0, double dec0 )
{
   // Returns the value of the ra at the equinox of date
   return ( ra0 +
            precess_m(t) +
            precess_n(t)*sin(ra_mid(t,ra0,dec0))*tan(dec_mid(t,ra0,dec0)) );
}

static double dec_date ( double t, double ra0, double dec0 )
{
   // Returns the value of the dec at the equinox of date
   return ( dec0 + precess_n(t)*cos(ra_mid(t,ra0,dec0)) );      
}

static double sun_mean_anomaly ( double t )
{
   // Returns mean anomaly of the Sun in radians given
   // the fraction of a Julian century from J2000.0
   double dr = M_PI / 180.0 ;
   double g;

   g = 357.528 + 35999.050 * t;
   g = g - int( g / 360.0 ) * 360.0;
   if ( g < 0.0 ) { 
     g = g + 360.0;
   }
   return( g*dr );
}

static double sun_eccentric_anomaly ( double g )
{
   // Returns the eccentric anomaly of the Sun in radians
   // given the mean anomaly of the Sun in radians
   // Iterate twice on Kepler's equation
   double ecc = 0.0167086171540; 
   double ea;

   ea = g + ecc * sin ( g ); 
   ea = g + ecc * sin ( ea ); 
   ea = g + ecc * sin ( ea ); 
   return ( ea );
}

static double sun_distance ( double t )
{
   // Return the distance of the Sun in AU given the fraction
   // of a Julian century from J2000.0
   double ecc = 0.0167086171540; 
   double ea;

   ea = sun_eccentric_anomaly( sun_mean_anomaly( t ) );
   return( 1.0000011 * ( 1.0 - ecc * cos( ea ) ) );
}

static double sun_mean_longitude ( double t )
{
   // Returns mean longitude of the Sun in radians given
   // the fraction of a Julian century from J2000.0
   double dr = M_PI / 180.0 ;
   double l;

   l = 280.460 + 36000.770 * t;
   l = l - int( l / 360.0 ) * 360.0;
   if ( l < 0.0 ) {
     l = l + 360.0;
   }
   return( l*dr );
}

static double sun_ecliptic_long ( double t )
{
   // Returns ecliptic longitude of the Sun in radians
   // given the fraction of a Julian century from J2000.0
   double dr = M_PI / 180.0 ;

   return( sun_mean_longitude( t ) +
           (1.915 * sin ( sun_mean_anomaly( t ) ) +
            0.020 * sin ( 2.0 * sun_mean_anomaly( t ) ) ) * dr );
}

static double obliquity_of_ecliptic ( double t )
{
   // Returns the obliquity of the ecliptic in radians
   // for the date given the fraction of a Julian century
   // from J2000.0
   double dr = M_PI / 180.0 ;

   return( ( 23.4393 - 0.01300 * t ) * dr );
}

static double ecl_lat ( double t, double ra, double dec )
{
   // Returns the ecliptic latitude of date in radians
   // given the ra and dec for the equinox of date as
   // well as the fraction of a Julian century from J2000.0
   double oe = obliquity_of_ecliptic ( t );

   return ( asin( sin(dec)*cos(oe) - cos(dec)*sin(oe)*sin(ra) ) );
}

static double ecl_lng ( double t, double ra, double dec, double beta )
{
   // Returns the ecliptic longitude of date in radians
   // given the ra and dec for the equinox of date as
   // well as the fraction of a Julian century from J2000.0
   // and the ecliptic latitude of date
   double oe = obliquity_of_ecliptic ( t );
   double sl, cl, el;

   sl = (sin(dec)*sin(oe) + cos(dec)*cos(oe)*sin(ra))/cos(beta);
   cl = cos(dec)*cos(ra)/cos(beta);
   el = atan2( sl, cl );
   if ( el < 0.0 ) {
      el = el + 2.0 * M_PI;
   }
   return ( el );
}

static double light_time ( double t, double ra, double dec )
{
   // Returns heliocentric correction in days. 
   double dist = sun_distance ( t );
   double beta = ecl_lat ( t, ra, dec );
   double lambda = ecl_lng ( t, ra, dec, beta);
   double lambda_sun = sun_ecliptic_long ( t );

   return ( -0.005775*dist*cos(beta)*cos(lambda - lambda_sun) );
}

double heliocentric( double t, double ra, double dec )
{
   // Returns heliocentric Julian date given Julian date and
   // ra and dec in radians

   double t_century = jul_cen( t );
   double ra_equinox = ra_date( t_century, ra, dec );
   double dec_equinox = dec_date( t_century, ra, dec );

   return t + light_time( t_century, ra_equinox, dec_equinox );
}
