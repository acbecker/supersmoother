/* Module:	Macho.h
 *		$Id: Macho.h,v 1.57 1996/02/24 00:25:22 jdr Exp $
 * Purpose:	Provide global definitions for Macho software
 * Copyright:	1991, Regents of the University of California
 *		This software may not be distributed to others without
 *		permission of the author.
 * Author:	Tim S. Axelrod , Lawrence Livermore National Laboratory
 *		tsa@llnl.gov 
 * 		Roberta A. Allsman, Lawrence Livermore National Laboratory
 *		robynallsman@llnl.gov
 * Modified:	$Log: Macho.h,v $
 * Modified:	Revision 1.57  1996/02/24 00:25:22  jdr
 * Modified:	added MAX_STARTAG_RECS
 * Modified:
 * Revision 1.56  1994/03/06  05:41:12  jdr
 * added definitions for SensorLog structure and bulletin board item
 *
 * Revision 1.55  1993/12/27  03:26:53  robyn
 * Updating repository for new Sun compilers for C and Fortran.
 *
 * Revision 1.54  1993/07/30  19:42:28  robyn
 * Installed another AcquiredFrom Star source location: HST originate but
 * finally positioned using SOD - for Bulge and SMC fiducial generation.
 *
 * Revision 1.53  1993/05/17  19:24:03  robyn
 * update for TileHit mod * SodGrouper
 *
 * Revision 1.52  1993/04/28  18:41:55  robyn
 * Increasing the size of the Field overlap boundary
 *
 * Revision 1.51  1993/04/26  22:28:17  robyn
 * sodo integration
 *
 * Revision 1.50  1993/03/08  18:14:11  robyn
 * Final ObsState DB upgrade software installation
 *
 * Revision 1.49  1993/02/26  19:48:30  robyn
 * More Sod changes
 *
 * Revision 1.48  1992/08/24  10:14:11  robyn
 * Restructure ObsState DB
 *
 * Revision 1.47  1992/08/23  04:33:20  robyn
 * *** empty log message ***
 *
 * Revision 1.46  1992/08/12  05:12:19  robyn
 * Specify the Std enum values.
 *
 * Revision 1.45  1992/07/20  01:07:47  robyn
 * Added enum TemplateAvailability for FieldDB
 *
 * Revision 1.44  1992/07/19  01:24:29  robyn
 * Alter MAX_STD_RECS to correspond to expected # records per tile.
 *
 * Revision 1.43  1992/07/16  06:21:59  robyn
 * Increase DB record size for DapObs.
 *
 * Revision 1.42  1992/07/14  05:34:49  robyn
 * Changed SideOfPier enum
 * Changed # records in Field DB and Std DB
 *
 * Revision 1.41  1992/07/07  05:08:59  tsa
 * 1. Moved enum typedefs from MachoBB.h to Macho.h
 *
 * Revision 1.40  1992/07/07  04:55:13  robyn
 * Changed Pier to SideOfPier
 *
 * Revision 1.39  1992/06/05  11:44:55  robyn
 * Change StdType ename names again
 *
 * Revision 1.38  1992/06/05  11:40:12  robyn
 * Alter names of StdType enum values.
 *
 * Revision 1.37  1992/06/05  11:30:07  robyn
 * Change typedef declaration for "c" compilation.
 *
 * Revision 1.36  1992/06/05  06:07:18  robyn
 * Modified space allocated fro overlapping tile list.
 *
 * Revision 1.35  1992/06/05  01:14:27  robyn
 * He changed his mind on the FIELD_DIMENSION_EDGE; now 30"
 *
 * Revision 1.34  1992/06/05  01:13:07  robyn
 * Change FIELD_DIMENSION to FIELD_EDGE to indicate size boundary around
 * 	region encompassing all a Field's chunks.
 *
 * Revision 1.33  1992/06/04  07:56:34  robyn
 * Changed StdType definition.
 *
 * Revision 1.32  1992/05/27  04:51:33  robyn
 * Change RA_DEC_PRECISION to 1000. for force float arithmetic
 *
 * Revision 1.31  1992/05/22  00:12:56  robyn
 * Moved enum defns from RTGeneric.H to Macho.h
 *
 * Revision 1.30  1992/05/13  11:41:11  robyn
 * Repaired RA min -> RA arc min use
 *
 * Revision 1.29  1992/05/12  03:38:56  tsa
 * 1. redefined RA_DEC_PRECISION, MAX_VALID_RA, MAX_VALID_DEC
 *
 * Revision 1.28  1992/04/10  23:52:17  tsa
 * 1. Increased MAX_ALARM_MSG (length of alarm message)
 *
 * Revision 1.27  1992/03/04  04:03:54  robyn
 * Reduced # records for ObsState since couldn't memory map > 58 records.
 * Added libmglobal.a containing malloc.caltech.o to get around mem align error
 * 	in system malloc
 *
 * Revision 1.26  1992/03/03  05:25:08  robyn
 * Reducing the file extension size.
 *
 * Revision 1.25  1992/03/01  01:56:22  robyn
 * Remove extraneous =
 *
 * Revision 1.24  1992/03/01  01:55:07  robyn
 * Add DIR_MODE and FILE_MODE to use for modes when creating files and directories
 *
 * Revision 1.23  1992/02/24  21:32:49  robyn
 * Added MAX_SODIN_BYTES
 *
 * Revision 1.22  1992/02/21  09:31:22  tsa
 * 1. Added MAX_SEEING_CHUNKS
 *
 * Revision 1.21  1992/02/15  03:08:18  robyn
 * Added new DB: Domain
 *
 * Revision 1.20  1992/02/14  05:09:24  robyn
 * Added INVALID_USHORT, tweeked value of INVALID_SHORT
 *
 * Revision 1.19  1992/02/12  22:29:19  robyn
 * Changed INVALID_MAG value to +29999
 *
 * Revision 1.18  1992/02/12  22:10:39  robyn
 * Added INVALID_MAG to designate illegal magnitude thruout DB.
 *
 * Revision 1.17  1992/02/11  21:48:49  robyn
 * Added INVALID_LONG as indicator that an int*4 was not loaded with "real" data.
 *
 * Revision 1.16  1992/02/09  00:53:23  robyn
 * Modified MAX_VALID_DEC, again. Fixed comments about it.
 *
 * Revision 1.15  1992/02/09  00:40:29  robyn
 * MOdified MAX_VALID_DEC.
 *
 * Revision 1.14  1992/02/07  07:39:21  robyn
 * Added RA/DEC max values and precision.
 *
 * Revision 1.13  1992/02/06  04:58:16  robyn
 * Added realistic MAX_RECS for each DB type.
 *
 * Revision 1.12  1992/02/05  00:28:19  robyn
 * Max PSF and Max Norm Star count once again tweeked--to 50 max.
 *
 * Revision 1.11  1992/02/03  00:42:33  robyn
 * Lowered default file extension size to 250Kbytes.
 *
 * Revision 1.10  1992/02/01  22:33:19  robyn
 * Added value defining precision of RA to tenth of arc second.
 *
 * Revision 1.9  1992/01/29  22:11:05  robyn
 * Change MAX_PSF_STARS into MAX_PSF_NORM_STARS to cover cumulative allocation
 * of PSF and Normalization star lists.
 *
 * Revision 1.8  1992/01/27  07:00:54  robyn
 * Change MAX_CHUNKS to fix Bennett's 512 pixel/chunk (172 min needed.
 * Change MAX_PSF_STARS to smaller
 *
 * Revision 1.7  1992/01/24  00:19:02  robyn
 * Change MAX_PIXELS to MAX_X_PIXELS MAX_Y_PIXELS
 *
 * Revision 1.6  1992/01/22  04:22:08  tsa
 * 1. Added new Error code definitions.
 *
 * Revision 1.5  1992/01/21  23:59:17  robyn
 * Testing CVS
 *
 * Revision 1.4  1992/01/21  23:27:57  robyn
 * Added MAX_PIXELS
 *
 * Revision 1.3  1992/01/19  23:51:29  robyn
 * Change name of INVALID.
 *
 * Revision 1.2  1992/01/19  05:39:53  tsa
 * Added #define for alarm message length
 *
 * Revision 1.1  1992/01/16  23:04:33  robyn
 * Initial revision
 *
 */

/*---------------------------------------------------------------------------*/
/*			    Macho.h	 				     */

#ifndef _Macho_h
#define	_Macho_h

/* ---------------- #defines of parameters ------------------------- */

#define MAX_ALARM_MSG	256
#define MAX_PATH_LEN    256
#define NAME_SIZE       256
#define MAX_FILE_NAME	256
	/* max # of characters allowed in filename */

#define MAX_COLORS	2
	/* number of color bands under review--e.g. red and blue */
typedef enum    {Red,Blue} Color;

typedef enum    {East=-1,West=1} SideOfPier;
	/* side of pier where telescope pointed  */

typedef enum {HasNoTemplate=0, HasDapOnly=1, HasSodOnly=2, HasDapSod=3}
	TemplateAvailability;
	/* indicates template availability status for a Field */

typedef enum {NoSlot=0, SlotDapOnly=1, SlotSodOnly=2, SlotDapSod=3}
	SlotField;
	/* indicates ObsTiles should be slotted for the Field */

typedef enum { NewObs, ReObs } ObsType;

typedef enum { DAP, SOD, NoPhot } PhotType;

#define MAX_CHUNKS	128
	/* max # allowed of (red + blue) chunks defined for any field */

#define MAX_SEEING_CHUNKS	1
	/* max # allowed of seeing chunks for any field */

#define MAX_TEMPLATE_TILE_OVERLAP 200
	/* max # of tiles whose stars are included in a field's templates */

#define MAX_PSF_STARS 50
	/* max # of stars needed for point spread function use by SOD */
#define MAX_NORM_STARS 50
	/* max # of stars needed for magnitude normalization use by SOD */

#define MAX_SODIN_BYTES 5000
	/* max # of bytes in Sod Input Parameter File stored in ObsState */
	/* stripped of comments */

#define MAX_DESCRIPTION 80
	/* max # characters allowed in descriptive string for sensors */
	/* this parameter is no longer used by sensors struct or code */

#define	SENSOR_BB_NAME		"SensorLog"
	/* bulletin board item name for dome sensors data; see MachoBB.h */
#define	MAX_SENSOR_NAME		10
#define	MAX_SENSOR_ALARMS	920
	/* maximum array sizes for strings in SensorLog */

#define INVALID_MAG   +29999
	/* value recognized a invalid photometry in db analysis processes */
#define INVALID_USHORT 59999
	/* value recognized a invalid thruout DB for ids, coords, etc */
#define INVALID_SHORT -29999
	/* value recognized a invalid thruout DB for ids, coords, etc */
#define INVALID_LONG -99999999
	/* value recognized a invalid thruout DB for ids, coords, etc */

#define MAX_VALID_RA	1296000000
#define MAX_VALID_DEC	324000000
	/* max valid integer for ra: +23:59:59.99999   dec: +-90:0:0.0 */
#define RA_DEC_PRECISION 1000.
	/* multiplicative factor to convert arc seconds to int values  - if you
	   change this value, MAX_VALID_RA and MAX_VALID_DEC must be changed
	   also!
	   */
	
#define FIELD_EDGE	60
	/* Field is defined as minimum region encompassing its chunks + 60 arc
	   second edge all around */

#define MAX_X_PIXELS	1024
#define MAX_Y_PIXELS	2048
	/* valid pixel value range:: x:(0,1023) y:(0,2047) */

#define NCAMSHORTS 28
	/* #shorts allocated to cam_config in ObsDesc structure*/
#define NCCDVARS 13
	/* there are 13 cam_config variables; to be numbered 0 through 12 */


#define FILE_MODE (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)
	/* always want plain files created with rw-rw--- */
#define DIR_MODE (S_IRWXU | S_IRWXG | S_ISGID)
	/* always want directories created with rwxrws-- */

#define FILE_EXTENSION_BLOCK_COUNT 120000
	/* add this many bytes every time file needs extension */

#define FILE_SUMMARY_SIZE 4000
	/* # of bytes in DB file summary block */

#define SMALL_FILE_SUMMARY_SIZE 512
	/* # of bytes in DB file summary block */

#define MAX_FIELD_RECS	500
#define MAX_STD_RECS	100
#define MAX_BP_RECS	1000
#define MAX_DOMAIN_RECS	10
#define MAX_TILE_RECS   35000
#define MAX_CCD_RECS	20
#define MAX_STAR_RECS	5000
#define MAX_TILEHIT_RECS 500
#define MAX_OBS_STATE_RECS 200
#define MAX_PHOTOBS_RECS 150
#define MAX_STARTAG_RECS 3000
	/* max # records allocated in DB TOC block; if too small requires 
	   read->write to new file with larger TOC block */

/* Error codes */

#define	ERRWARN_BBOD		4 
#define	ERRFATAL_BBOD		5


/* --------------- Epoch designators for Star Coordinates ------------------ */
#ifndef _Epoch_db_
#define  _Epoch_db_
typedef enum { Epoch_J2000,Epoch_B1950 } StarEpoch;
#endif /* _Epoch_db_ */

/* --------------- Binary value indicating sign ---------------------------- */
#ifndef _Sign_db_
#define  _Sign_db_
typedef enum { Positive, Negative } Sign;
#endif /* _Sign_db_ */

/* --------------- Identifies fiducial selection type for a star ------------*/
#ifndef _Fiducial_db_
#define  _Fiducial_db_
#define STD_ASTROMETRY 2
#define STD_PHOTOMETRY 1
typedef enum { StdTypeNone=0,StdTypePhotometry=1,StdTypeAstrometry=2,
	StdTypeAstroPhoto=3 } StdType;
#endif /* _Fiducial_db_ */

/* -------------- Identifies DB use of Star ---------------------------------*/
#ifndef _StarClassification_db_
#define  _StarClassification_db_
typedef enum { ClassificationNone, ClassificationFiducial } StarClassification;
#endif /* _StarClassification_db_ */

/*-------------- Identifies where star was originally acquired --------------*/
#ifndef _StarAcquired_db_
#define  _StarAcquired_db_
typedef enum { AcquiredOther, AcquiredCosmos, AcquiredDap, AcquiredSod,
		AcquiredHST_Sod} StarAcquired;
#endif /* _StarAcquired_db_ */

#endif /* _Macho_h */
