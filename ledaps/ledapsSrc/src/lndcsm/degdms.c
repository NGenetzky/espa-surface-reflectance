
/*****************************************************************************
	NAME:  DEGDMS

	PURPOSE:  To convert total degrees, total minutes, or total seconds to 
		  packed degress,minutes,seconds.

	PROGRAM HISTORY:
	PROGRAMMER	DATE	REASON
	----------	----	------
	G. Schmidt  3/2013  Borrowed code from geolib for LEDAPS

    RETURNS:
        false - error occured processing
        true -  successful completion

	ALGORITHM:
    Receive an angle in seconds, minutes, or degrees
	Convert it to DMS.  
	The angle is then checked to be sure it is within the limits 
	of its use(LAT, LON, or DEGREES).

******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "space.h"
#include "error.h"
#include "bool.h"

static int find_deg( double angle	/* Angle in total degrees */
     );

static int find_min( double angle	/* Angle in total degrees */
     );

static double find_sec( double angle	/* Angle in total degrees */
     );

bool degdms( double *deg, double *dms, char *code, char *check )
{			/* degdms */

    double tsec;
    double MAXDMS;
    double MAXMIN = 60060;
    double MAXSEC = 60;
    double MINDMS;
    double MINMIN = -60060;
    double MINSEC = -60;
    double tempmin;
    double tempsec;

    long tdeg;
    long tmin;
    long sign;

    if ( strcmp( check, "LAT" ) == 0 )
    {
	MAXDMS = 90000000;
	MINDMS = -90000000;
    }
    else if ( strcmp( check, "LON" ) == 0 )
    {
	MAXDMS = 180000000;
	MINDMS = -180000000;
    }
    else
    {
	MAXDMS = 360000000;
	MINDMS = 0;
    }

    if ( strcmp( code, "DMS" ) != 0 )
    {
	if ( strcmp( code, "MIN" ) == 0 )
	{
	    *deg = *deg / 60.0;
	}
	if ( strcmp( code, "SEC" ) == 0 )
	{
	    *deg = *deg / 3600.0;
	}
	tdeg = ( long ) find_deg( *deg );
	tmin = ( long ) find_min( *deg );
	tsec = find_sec( *deg );
	sign = 1;
	if ( *deg < 0 )
	    sign = -1;
	tdeg = abs( tdeg );
	*dms = ( ( tdeg * 1000000 ) + ( tmin * 1000 ) + tsec ) * sign;

/*  Check to be sure coordinate is valid
----------------------------------------*/
	if ( ( *dms > MAXDMS ) || ( *dms < MINDMS ) )
	{
        RETURN_ERROR("Invalid coordinate value", "c_degdms", false);
	}
    }
    else
    {
	*dms = *deg;

/*  Check to be sure coordinate is valid
----------------------------------------*/
	if ( ( *dms > MAXDMS ) || ( *dms < MINDMS ) )
	{
        RETURN_ERROR("Invalid coordinate value", "c_degdms", false);
	}

	if ( ( strcmp( check, "LAT" ) != 0 ) && ( strcmp( check, "LON" ) != 0 ) )
	{
	    if ( *dms <= 0 )
	    {
        RETURN_ERROR("Invalid coordinate value", "c_degdms", false);
	    }
	}

/*  parse out the minutes value from DMS and check against MAXMIN
-----------------------------------------------------------------*/
	tempmin = *dms - ( ( ( int ) ( *dms / 1000000 ) ) * 1000000 );
	if ( ( tempmin > MAXMIN ) || ( tempmin < MINMIN ) )
	{
        RETURN_ERROR("Invalid coordinate value", "c_degdms", false);
	}

/*  parse out the seconds value from DMS and check against MAXSEC
-----------------------------------------------------------------*/
	tempsec = *dms - ( ( ( int ) ( *dms / 1000 ) ) * 1000 );
	if ( ( tempsec > MAXSEC ) || ( tempsec < MINSEC ) )
	{
        RETURN_ERROR("Invalid coordinate value", "c_degdms", false);
	}
    }

    return ( true );
}				/* degdms */

/****************************************************************************
NAME:			find_deg, find_min, find_sec	

PURPOSE:  Extracts deg, min, or sec portions of an angle

PROGRAM HISTORY:
VERSION	 DATE	AUTHOR	   CODE/CONT   REASON
-------	 ----	------	   ---------   -----------------------------
  5.0	 5/89	D. Steinwand  CSB      LAS 5.0 (original) Development

COMPUTER HARDWARE AND/OR SOFTWARE LIMITATIONS:   must be run under TAE

PROJECT:  LAS

ALGORITHM DESCRIPTION:
	Extract portion of angle
	Return
*****************************************************************************/
static int find_deg( double angle	/* Angle in total degrees */
     )
{
    long sign;			/* Sign of angle */
    long deg;			/* Degrees portion of angle */
    long minute;		/* Minutes portion of angle */
    double sec;			/* Seconds portion of angle */

    deg = ( long ) angle;
    sign = 1;
    if ( deg < 0 )
	sign = -1;
    deg = ( long ) fabs( angle );
    minute = ( long ) ( ( fabs( angle ) - deg ) * 60.0 );
    sec = ( ( ( fabs( angle ) - deg ) * 60.0 ) - minute ) * 60.0;
    if ( sec >= 59.999 )
	minute++;
    if ( minute >= 60 )
	deg++;
    deg *= sign;
    return ( deg );
}

static int find_min( double angle	/* Angle in total degrees */
     )
{
    double sec;			/* Seconds portion of angle */
    long minute;		/* Minutes portion of angle */

    angle = fabs( angle );
    angle -= ( long ) angle;
    minute = ( long ) ( angle * 60.0 );
    sec = ( ( angle * 60.0 ) - minute ) * 60.0;
    if ( sec > 59.999 )
	minute++;
    if ( minute >= 60 )
	minute -= 60;
    return ( minute );
}

static double find_sec( double angle	/* Angle in total degrees */
     )
{
    long temp_angle;

    angle = fabs( angle );
    angle -= ( long ) angle;
    angle *= 60.0;
    angle -= ( long ) angle;
    angle *= 60.0;
    if ( angle > 59.999 )
	angle -= 60.0;
    temp_angle = ( long ) ( angle * 1000.0 );	/* Truncate to 0.001 sec */
    angle = temp_angle / 1000.0;
    return ( angle );
}
