#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <Macho.h>
#include <libsite.H>
#include <libclif.H>
#include "phase.H"
#include "fit/lcfit.H"

/* fits a single Lc, derived from fazAndy4 */
/* super basic : just phase a Lc and print out periods */

CommonLc *Lc;

static void usage(void);

int main( int argc, char *argv[] ) {
    
    clifpar.load( );
    ReimannPhaser rp( clifpar );
    
    int minpts = 20;
    int nobs, np, i, influx=0;
    double period, fixp0, per, p0=0;
    float sigma, tmin=-1e10, tmax=1e10;
    char whichper[8];
    
    char cmd[1024];
    char buf[1024];
    
    fixp0 = 0;
    period = 0;
    sigma  = 0;
    sprintf(whichper, "avg");
    sprintf(cmd, "%s", argv[0]);
    for (i = 1; i < argc; i++) 
        sprintf(cmd, "%s %s", cmd, argv[i]);
    sprintf(cmd, "%s\n", cmd);
    
    while(argc > 1 && (argv[1][0] == '-' || argv[1][0] == '+')) {
        switch (argv[1][1]) {
        case '?':
        case 'h':
            usage();
            exit(0);
            break;
        case 'f':
            influx = 1;
            break;
        case 'p':
            if(argc == 2) {
                fprintf(stderr, "-p requires a number");
                exit(1);
            } else {
                argc--; argv++;
                period = atof(argv[1]);
            }
            break;
        case 't':
            if(argc == 2) {
                fprintf(stderr, "-t requires a number");
                exit(1);
            } else {
                argc--; argv++;
                fixp0 = atof(argv[1]);
            }
            break;
        case 'm':
            if(argc == 2) {
                fprintf(stderr, "-m requires a number");
                exit(1);
            } else {
                argc--; argv++;
                minpts = atoi(argv[1]);
            }
            break;
        case 'k':
            if(argc == 2) {
                fprintf(stderr, "-k requires a number");
                exit(1);
            } else {
                argc--; argv++;
                rp.KNOT = atoi(argv[1]);
            }
            break;
        case 's':
            if(argc == 2) {
                fprintf(stderr, "-s requires a number");
                exit(1);
            } else {
                argc--; argv++;
                sigma = atof(argv[1]);
            }
            break;
            
        case 'n':
            if(argc == 2) {
                fprintf(stderr, "-n requires a number");
                exit(1);
            } else {
                argc--; argv++;
                rp.NMIN = atoi(argv[1]);
            }
            break;
        case 'i':
            if(argc == 2) {
                fprintf(stderr, "-i requires a number");
                exit(1);
            } else {
                argc--; argv++;
                rp.MINP = atof(argv[1]);
            }
            break;
        case 'a':
            if(argc == 2) {
                fprintf(stderr, "-a requires a number");
                exit(1);
            } else {
                argc--; argv++;
                rp.MAXP = atof(argv[1]);
            }
            break;
        case 'r':
            if(argc == 2) {
                fprintf(stderr, "-r requires a number");
                exit(1);
            } else {
                argc--; argv++;
                rp.RATE = atof(argv[1]);
            }
            break;
        case 'b':
            if(argc == 2) {
                fprintf(stderr, "-b requires a number");
                exit(1);
            } else {
                argc--; argv++;
                rp.BASS  = atof(argv[1]);
            }
            break;
        case '0':
            if(argc == 2) {
                fprintf(stderr, "-0 requires a number");
                exit(1);
            } else {
                argc--; argv++;
                tmin = atof(argv[1]);
            }
            break;
        case '1':
            if(argc == 2) {
                fprintf(stderr, "-1 requires a number");
                exit(1);
            } else {
                argc--; argv++;
                tmax = atof(argv[1]);
            }
            break;
            
        default:
            fprintf(stderr, "Unknown option %s\n", argv[1]);
            exit(1);
            break;
        }
        argc--;
        argv++;
    }
    char *infile = argv[1];
    
    /* read in ascii lightcurve */
    FILE *fi = fopen( infile, "r" );
    if ( !fi ) return 0;
    
    nobs = 0;
    int line = 1;
    float mag, err;
    double date;
    while ( fgets(buf,1024,fi) ) {
        if ( isalpha(buf[0]) || ispunct(buf[0]) ) continue;
        else ++nobs;
    }

    rewind( fi );
    
    Lc = new CommonLc( infile, 0 );
    if (influx) {
        Lc->norm1  = 0;
        Lc->units  = CommonLc::Linear;
    }
    else {
        Lc->norm1  = 25;
        Lc->units  = CommonLc::Magnitude;
    }
    Lc->resize( nobs );
    
    nobs = 0;
    line = 1;
    while ( fgets(buf,1024,fi) ) {
        if ( isalpha(buf[0]) || ispunct(buf[0]) ) continue;
        if ( sscanf( buf, "%lf %f %f", &date, &mag, &err) == 3 ) {
            if (date > tmin && date < tmax) {
                Lc->date[ nobs ]  = date;
                Lc->value[ nobs ] = mag;
                Lc->error[ nobs ] = err;
                Lc->obsid[ nobs ] = line;
                if ((nobs > 1) && (Lc->date[ nobs ] < Lc->date[ nobs-1 ])) {
                    fprintf(stderr, "Must be monotonically increasing!  %f %f\n", Lc->date[ nobs ], Lc->date[ nobs-1 ]);
                    fprintf(stderr, "  %s\n", buf);
                    fprintf(stderr, "  %f %f %f\n", date, mag, err);
                    exit(1);
                }
                ++nobs;
            }
        }
        ++line;
    }
    Lc->nobs = nobs;
    fclose( fi );
    
    if (nobs < minpts) {
        fprintf(stderr, "Too Few Points to Fit (%d < %d)\n", nobs, minpts);
        exit(1);
    }
    
    if (!(influx)) {
        fprintf(stdout, "Converting from Magnitudes to flux (if input in flux, use -f)\n");
        Lc->convert( CommonLc::Linear );
    }
    else {
        fprintf(stdout, "Leaving data in flux units\n");
    }
    
    
    sprintf(buf, "%s_Phase", Lc->title);
    
    if (period == 0.) {
        /* find the periods */
        np = 0;
        per = 0;
        if ( Lc ) {
            fprintf(stdout, "Phasing \n");
            //fprintf(stderr, "Data 0 : %lf %f %f\n", Lc->date[0], Lc->value[0], Lc->error[0]);
            if( rp.phase( *Lc ) ) {
                fprintf(stdout, "Could not determine period\n");
            } else {
                ++np;
                per = 1.0/rp.FF[0];
            }
            fprintf(stdout, "Period : %.8lf\n", per);
            fprintf(stdout, "  %13.9f %8.4f %8.4f\n  %13.9f %8.4f %8.4f\n  %13.9f %8.4f %8.4f\n  %13.9f %8.4f %8.4f\n  %13.9f %8.4f %8.4f\n",
                    1.0 / rp.FF[0], rp.mm[0], rp.TT[0],
                    1.0 / rp.FF[1], rp.mm[1], rp.TT[1], 
                    1.0 / rp.FF[2], rp.mm[2], rp.TT[2],
                    1.0 / rp.FF[3], rp.mm[3], rp.TT[3],
                    1.0 / rp.FF[4], rp.mm[4], rp.TT[4]);
            
        }
    }
}

static void usage(void)
{
    char **line;
    
    static char *msg[] = {
        "Usage: fazAndy [options] 2Mass_data_file",
        "",
        "Your options are:",
        "       -h      This message",
        "       -m      Minimum number of data points : 20",
        "       -p      Fix period (fast)",
        "       -t      Fix t0 (use with -p)",
        "       -s      Sigma clip output data",
        "       -n      Phaser NMIN : 20",
        "       -a      Phaser MAXP : -1 (full timeseries)",
        "       -i      Phaser MINP : 0.1 (halving this doubles runtime)",
        "       -r      Phaser RATE : 6.0 (doubling doubles runtime)",
        "       -b      Phaser BASS : 0.0",
        "       -0      Minimum Date",
        "       -1      Maximum Date",
        NULL,
    };
    
    for(line = msg;*line != NULL;line++) {
        fprintf(stderr,"%s\n",*line);
    }
}
