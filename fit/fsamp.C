#include <stdlib.h>
#include "FiniteSourceAmp.H"

int main( int argc, char *argv[] )
{
   if( argc < 3 ) {
      printf("\nUsage:\n\t%s [u] [u_star]\n", argv[0] );;
      exit(1);
   }

   double u = 0;
   double ustar = 0;
   double amp = 0;
   
   u = atof( argv[1] );
   ustar = atof( argv[2] );

   amp = FS_Amp(u, ustar, 0, 0);

   //printf("%f %f %f\n", u, ustar, amp);
   printf("%f\n", amp);
}
