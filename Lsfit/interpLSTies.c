#include "stdio.h"
#include "string.h"
#include "math.h"
#include <stdlib.h>
#include "geotiff/xtiffio.h"   /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "landsatSource64/clib/standard.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"

static void lltoxy1(double alat, double alon, double *x, double *y, double dlam, double slat);
static int32_t interpTiePoint(matchResult *matches, lsTiepoints *tiePoints, int32_t j);

/*
  Read tiepoint file for tiepoints.
*/
void interpLSTies(lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints)
{
	int32_t i, j;
	double minX, maxX, minY, maxY;
	double psScale;
	/*
	  Cycle through ties
	*/
	minX = matches->x0;
	maxX = matches->x0 + (matches->nx - 1) * matches->dx * (double)matches->stepX;
	minY = matches->y0;
	maxY = matches->y0 + (matches->ny - 1) * matches->dy * (double)matches->stepY;
	j = 0;
	fprintf(stderr, "%lf %lf %lf %lf\n", minX, maxX, minY, maxY);

	for (i = 0; i < tiePoints->npts; i++)
	{
		if (fitDat->proj == NSIDCNORTH)
		{
			lltoxy1(tiePoints->lat[i], tiePoints->lon[i], &(tiePoints->x[i]), &(tiePoints->y[i]), ROTNORTH, SLATNORTH);
		}
		else if (fitDat->proj == NSIDCSOUTH)
		{
			lltoxy1(tiePoints->lat[i], tiePoints->lon[i], &(tiePoints->x[i]), &(tiePoints->y[i]), ROTSOUTH, SLATSOUTH);
		}
		else
			error("Invalid Projection");
		/*
		  Reject points out of box
		*/
		tiePoints->x[i] *= KMTOMS;
		tiePoints->y[i] *= KMTOMS;
		if ((tiePoints->x[i] >= minX) && (tiePoints->x[i] <= maxX) && (tiePoints->y[i] >= minY) && (tiePoints->y[i] <= maxY))
		{
			/* Move good points forward in arrays */
			tiePoints->x[j] = tiePoints->x[i];
			tiePoints->y[j] = tiePoints->y[i];
			tiePoints->z[j] = tiePoints->z[i];
			tiePoints->lat[j] = tiePoints->lat[i];
			tiePoints->lon[j] = tiePoints->lon[i];
			tiePoints->vx[j] = tiePoints->vx[i];
			tiePoints->vy[j] = tiePoints->vy[i];
			/*
			   Interpolate specific point
			   Note: if not point interpolated, j is not updated, so that element will overwritten on next loop though.
			   This means final array is only good points.
			*/
			if (interpTiePoint(matches, tiePoints, j) == TRUE)
			{
				/*
				   Convert to real meters from PS meters
				*/
				psScale = xyscale(tiePoints->lat[j], fitDat->proj);
				tiePoints->offXT[j] *= psScale;
				tiePoints->offYT[j] *= psScale;
				tiePoints->xyScale[j] = psScale;
				/* if(tiePoints->vy[i] !=0)     fprintf(stderr,"%lf %lf  ....",tiePoints->offXT[j],tiePoints->offYT[j]) ;*/
				/*
				  correct by subtracting displacement computed from tiepoint, note everything is pixels
				*/
				tiePoints->offXT[j] -= ((tiePoints->vx[i] / 365.25) * fitDat->deltaT) / matches->dx;
				tiePoints->offYT[j] -= ((tiePoints->vy[i] / 365.25) * fitDat->deltaT) / matches->dy;
				/* if(tiePoints->vy[i] !=0)  fprintf(stderr,"%lf %lf %lf %lf\n",tiePoints->offXT[j],tiePoints->offYT[j],
				   (tiePoints->vx[i]/365.25)*fitDat->deltaT, (tiePoints->vy[i]/365.25)*fitDat->deltaT) ;		*/
				j++;
			}
		}
	}
	tiePoints->npts = j;
}

static int32_t interpTiePoint(matchResult *matches, lsTiepoints *tiePoints, int32_t j)
{
	double xi, yi, t, u;
	double p1, p2, p3, p4;
	int32_t jm, im;
	/*
	  Image coordinates
	*/
	xi = (tiePoints->x[j] - matches->x0) / (matches->dx * (double)matches->stepX);
	yi = (tiePoints->y[j] - matches->y0) / (matches->dy * (double)matches->stepY);
	/* Integer values, rounded down */
	jm = (int32_t)xi;
	im = (int32_t)yi;
	/* Return if out of bounds, use -2 since  im+1,jm+1 index array,means extreme border points neglected (good thing) */
	if (jm < 0 || im < 0 || jm > (matches->nx - 2) || im > (matches->ny - 2))
		return (FALSE);
	/* Fraction of pixel for interpolation */
	t = (float)(xi - (double)jm);
	u = (float)(yi - (double)im);
	/*
	  X point
	*/
	p1 = matches->X[im][jm];
	p2 = matches->X[im][jm + 1];
	p3 = matches->X[im + 1][jm + 1];
	p4 = matches->X[im + 1][jm];
	/* Don't use if all 4pts aren't good - should ensure better quality data */
	if (p1 <= (NODATA + 1) || p2 <= (NODATA + 1) || p3 <= (NODATA + 1) || p4 <= (NODATA + 1))
		return (FALSE);
	tiePoints->offXT[j] = (double)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);
	/*
	  Y point
	*/
	p1 = matches->Y[im][jm];
	p2 = matches->Y[im][jm + 1];
	p3 = matches->Y[im + 1][jm + 1];
	p4 = matches->Y[im + 1][jm];
	/* Don't use if all 4pts aren't good - should ensure better quality data */
	if (p1 <= (NODATA + 1) || p2 <= (NODATA + 1) || p3 <= (NODATA + 1) || p4 <= (NODATA + 1))
		return (FALSE);
	tiePoints->offYT[j] = (double)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);

	return (TRUE); /* Return true since interpolation appears to have worked */
}

static void lltoxy1(double alat, double alon, double *x, double *y, double dlam, double slat)
{
	/*
		Converts from geodetic latitude and longitude to Polar
		Stereographic (x,y) coordinates for the Polar regions.
		The standard parallel (lat with no distortion) is 70 deg.
		Equations are from Snyder, J.P. 1982, Map Projections Used
		by the U.S. Geological Survey, Geological Survey Bulletin
		1532, U.S. Gov. Printing Office.  See JPL Tech. Memo.
		3349-85-101 for details.  Sub written by C.S. Morris-
		Apr 29, 1985 (Goddard?).
		converted to MASSCOMP 7/7/87 DRT
		on SUN      11/89  EAF
		to C        03/94 IRJ

		alat in	latitude, degrees, +N, -S
		alon in     longitude, degrees 0-360(East-West)
		x	 out	x-coordinate in km
		y	 out	y-coordinate in km
		dlam in	rotation, degrees, 0-360
	*/
	int i1;
	double t[2];
	double e, e2, re, sn, cm, rho, rlat, tmp;
	/*
	  Radius of earth (km) -Hughes Ellipsoid
	  re=6378.273
	*/
	re = 6378.137;
	/*
	  Eccentricity of earth -Hughes Ellipsoid	  e2=0.006693883  changed to WGS 84 10/14/05
	*/
	e2 = 0.0066943801;
	e = sqrt(e2);
	/*
	  Standard parallel  slat=70;

	  Test for N or S hemi, set constants as necessary
	  For SSM/I grid, Northern Hemisphere  sn=1.
	*/
	if (alat < 0)
		sn = -1.;
	else
		sn = 1.;
	/*
	  Compute x, y
	*/
	alat = sn * alat;
	alon = sn * alon;
	if (alat >= 89.995)
	{
		*x = 0.;
		*y = 0.;
	}
	else
	{
		rlat = alat;
		for (i1 = 0; i1 < 2; i1++)
		{
			if (i1 == 1)
				rlat = slat;
			t[i1] = tan((PI / 4.) - (rlat / (2. * RTOD)));
			tmp = pow((1. - e * sin(DTOR * rlat)) / (1. + e * sin(DTOR * rlat)), (e / 2.0));
			t[i1] = t[i1] / tmp;
		}
		cm = cos(DTOR * slat) / sqrt(1. - e2 * pow(sin(DTOR * slat), 2.0));
		rho = re * cm * t[0] / t[1];
		*x = rho * sn * sin(DTOR * (alon + dlam));
		*y = -rho * sn * cos(DTOR * (alon + dlam));
	}
	return;
}
