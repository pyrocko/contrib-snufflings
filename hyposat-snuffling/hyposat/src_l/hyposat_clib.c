/*
@FORTRAN callable function to convert human time format to epochal time
and fill in any missing human format time fields, needs enough info to
compute julian date if julian date is not supplied.
Author: T. McElfresh
Modified for BUG use: Dierk Wand, Ruhr-University Bochum, 18. Jan. 1989
*/


/* Author: T. McElfresh */

/* FORTRAN callable function to convert "human" to epochal time */

/* Input date hour minute second */
/*    or year month/mname day hour minute second */
/*    or year doy hour minute second */

/* Assigns epochal time & other unfilled fields */

struct date_time{
	double epoch;
	long date;
	int year;
	int month;
	char mname[4];
	int day;
	int doy;
	int hour;
	int minute;
	float second;
};

int  fhtoe_(epoch,date,year,month,mname,day,doy,hour,minute,second)

	double *epoch;
	long *date;
	int *year;
	int *month;
	char mname[4];
	int *day;
	int *doy;
	int *hour;
	int *minute;
	float *second;

   {
     struct date_time dp;

     dp.date = *date;
     dp.year = *year;
     dp.month = *month;
     dp.mname[0] = mname[0];
     dp.mname[1] = mname[1];
     dp.mname[2] = mname[2];
     dp.mname[3] = mname[3];
     dp.day = *day;
     dp.doy = *doy;
     dp.hour = *hour;
     dp.minute = *minute;
     dp.second = *second;

     /* htoe needs julian day, so calculate it if necesary */
     if(*date <= 0)
       {
       if(dp.year < 100)
	 dp.year = dp.year + 1900;
       if(*doy <= 0)
         mdtodate(&dp);
       dp.date = (dp.year*1000)+dp.doy;
       }

     htoe(&dp);

     *epoch = dp.epoch;
	*date = dp.date;
	*year = dp.year;
	*month = dp.month;
	mname[0] = dp.mname[0];
	mname[1] = dp.mname[1];
	mname[2] = dp.mname[2];
	mname[3] = dp.mname[3];
	*day = dp.day;
	*doy = dp.doy;
	*hour = dp.hour;
	*minute = dp.minute;
	*second = dp.second;
     return(1);
   }
/****************************************************************************************************************************
***
***                   [source file "SRC/LIB/LIBUTILS/fetoh_.c"]
***                             int function fetoh_
***
***    Call this function from FORTRAN to convert an "epochal" time into a "human" time.
***
***    Written by Wilmer Rivers on June 26, 1985
***
***    This function calls function "etoh", which resides in /usr/CSS/libtime.a   (link it with the object file by specifying
***    "-ltime" to the loader).
***
****************************************************************************************************************************/

int fetoh_ (epoch, date, year, month, mname, day, doy, hour, minute, second)
   double *epoch;
   long *date;
   int *year, *month, *day, *doy, *hour, *minute;
   char mname[4];
   float *second;
{
   struct date_time
      { double epoch;
        long date;
        int year;
        int month;
        char mname[4];
        int day;
        int doy;
        int hour;
        int minute;
        float second; } dp;

   dp.epoch = *epoch;

   etoh(&dp);

   *date = dp.date;
   *year = dp.year;
   *month = dp.month;
   mname[0] = dp.mname[0];
   mname[1] = dp.mname[1];
   mname[2] = dp.mname[2];
   mname[3] = dp.mname[3];
   *day = dp.day;
   *doy = dp.doy;
   *hour = dp.hour;
   *minute = dp.minute;
   *second = dp.second;

   return (1);
}

#define ISLEAP(yr)      (!(yr % 4) && yr % 100 || !(yr % 400))
#define isleap(yr)      (!(yr % 4) && yr % 100 || !(yr % 400))

static int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31,31};

htoe(dt)
register struct date_time *dt;
{
	double dtoepoch();
	dt->epoch = 
	dtoepoch(dt->date) + 
	dt->hour * 3600. + 
	dt->minute * 60. +
	dt->second;
	return(0);
}
mdtodate(dt)
register struct date_time *dt;
{
	int i,dim;
	dt->doy = 0;
	for( i = 0 ; i < dt->month - 1 ; i++ ){
		dim = days_in_month[i];
		if( i == 1 && ISLEAP(dt->year) ) dim++;
		dt->doy += dim;
	}
	dt->doy += dt->day;
	dt->date = 1000 * dt->year + dt->doy;
	return(0);
}

#define ISLEAP(yr)      (!(yr % 4) && yr % 100 || !(yr % 400))
#define isleap(yr)      (!(yr % 4) && yr % 100 || !(yr % 400))

char *strcpy();

static char *month_name[] =
{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

#define mod(a,b)	(a) - ((int)((a)/(b))) * (b)
etoh(dt)
register struct date_time *dt;
{
	int diy;
	double secleft;

	dt->doy = dt->epoch / 86400.;
	secleft = mod(dt->epoch,86400.0);
	dt->hour = dt->minute = dt->second = 0;

	if(secleft) {			/* compute hours minutes seconds */
		if(secleft < 0) {	/* before 1970 */
			dt->doy--;		/* subtract a day */
			secleft += 86400;	/* add a day */
		}
		dt->hour = secleft/3600;
		secleft = mod(secleft,3600.0);
		dt->minute = secleft/60;
		dt->second = mod(secleft,60.0);
	}

	if(dt->doy >= 0){
		for( dt->year = 1970 ; ; dt->year++ ){
			diy = ISLEAP(dt->year) ? 366:365;
			if( dt->doy < diy ) break;
			dt->doy -= diy;
		}
	}
	else{
		for( dt->year = 1969 ; ; dt->year-- ){
			diy = ISLEAP(dt->year) ? 366:365;
			dt->doy += diy;
			if( dt->doy >= 0 ) break;
		}
	}
	dt->doy++;
	dt->date = dt->year * 1000 + dt->doy;
	month_day(dt);
	return;
}
month_day(dt)
register struct date_time *dt;
{
	int i,dim,leap;

	leap = ISLEAP(dt->year);
	dt->day = dt->doy;
	for( i = 0 ; i < 12 ; i ++ ){
		dim = days_in_month[i];
		if( leap && i == 1 ) dim++;
		if( dt->day <= dim ) break;
		dt->day -= dim;
	}
	dt->month = i + 1;
	strcpy(dt->mname,month_name[i]);
}

/*
 * collection of time conversion utility subroutines
 */
#include <stdio.h>
#define TIMENULL -99999999999.999
 
/*
 * convert julian date to epoch time
 */
double
dtoepoch(date)
register long	date;
{
	register int	cnt;
	double	days;

	cnt = (int)(date / 1000);
	days = 0;
	if (cnt > 1970)
		while (--cnt >= 1970)
			days += ISLEAP(cnt) ? 366 : 365;
	else if (cnt < 1970)
		while (cnt < 1970) {
			days -= ISLEAP(cnt) ? 366 : 365;
			cnt++;
		}
	return((days + (date - 1) % 1000) * 86400.);
}

/*
 * timecon --
 *	convert hh[:mm[:ss[.sss]]] to epoch time
 */
double
timecon(timstr)
register char	*timstr;
{                                       /* seconds conversion factors */
	static float	cinst[3]={3600.,60.,1.};
	register char	*C;		/* traveling pointer */
	register short	colon,		/* colons in string */
			dot;		/* decimal points in string */
	double	tnum,			/* total number of seconds */
		atof();
	char	*s;			/* save input timstr pointer */

	s = timstr;
	dot = colon = 0;
	tnum = 0.0;
	for (C = timstr;*C;++C)
		if (!isdigit(*C))
			if (*C == ':') {
				if (dot || colon > 1) {
					fprintf(stderr,"timecon: bad time string %s\n",s);
					return(TIMENULL);
				}
				tnum += atof(timstr) * cinst[colon];
				timstr = C + 1;
				++colon;
			}
			else if (*C == '.' && dot++)  {
				fprintf(stderr,"timecon: bad time string %s\n",s);
				return(TIMENULL);
			}
	return(tnum + atof(timstr) * cinst[colon]);
}
