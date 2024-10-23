#include "stdio.h"
#include "math.h"
#include "string.h"
#include <stdlib.h>
/*#include "geotiff/xtiffio.h"  for TIFF
#include "geotiff/geotiffio.h"  for GeoTIFF */
#include "landsatSource64/clib/standard.h"
#include "landsatSource64/cRecipes/cRecipes.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "landsatSource64/cRecipes/nrutil.h"
#define MAXTIEPOINTS 500000
#define REJECTTHRESHOLD 100 // Increased from 50 to 100 7/25/2024
#define MAXTIEITERATIONS 3

typedef struct
{
	double x;
	double y;
} svdData;

static void planeFit(void *x, int32_t i, double *afunc, int32_t ma);
static void constantFit(void *x, int32_t i, double *afunc, int32_t ma);
static void printCovariance(double **CmatrixX, int32_t ma, FILE *fpOut, char *coord);
static void printFit(double *aX, double *aY, int32_t ma, FILE *fpOut);
static void getCoeffs(double *aX, double *p1, double *p2, double *p3, int32_t ma);
static void computePointsWidth(svdData *xy, int32_t nPts, double *xWidth, double *yWidth);
static void computeStats(lsFit *fitDat, matchResult *matches,
						 svdData *xy, double *dx, double *dy, double *aX, double *aY,
						 double *meanX, double *meanY, double *sigmaX, double *sigmaY, double *varX, double *varY,
						 int32_t nPts, int32_t ma);
static void discardOutliers(lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints,
							svdData *xy, double *dx, double *dy, double *aX, double *aY,
							double sigmaX, double sigmaY, int32_t *nReject, int32_t ma);

static void printFitResult(lsFit *fitDat, double *aX, double *aY, double **CmatrixX, double **CmatrixY, int32_t ma,
						   double meanX, double meanY, double sigmaX, double sigmaY,
						   int32_t nReject, int32_t nRejectForFit, int32_t nptsOrig, int32_t tieCountForFit);
/*
  Read tiepoint file for tiepoints.
*/
void computeLSFit(lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints)
{
	void (*fitFunction)(void *, int i, double *, int);
	double **v, **u, **CmatrixX, **CmatrixY, chisq, *sigB, *sigX, *sigY, *w; /* Svd fit variabiles */
	double *aX, *aY;														 /* Solution for params */
	double *dx, *dy, *dxOrig, *dyOrig;										 /* Offsets data moved to 1...n vector */
	double resX, resY;
	double varX, varY, sigmaX, sigmaY, meanX, meanY;
	double xWidth, yWidth;
	int32_t nData, ma;
	int32_t nptsOrig;
	svdData *xy, *xyOrig;
	int32_t nIterations, notDone; /* Used to control iterations loop */
	int32_t nReject, nRejectForFit, tieCountForFit;
	int32_t i, j;

	/*
	   Init arrays for LS fit
	*/
	aX = (double *)malloc(sizeof(double) * 20); /* 20 just incase parameters increase later */
	aY = (double *)malloc(sizeof(double) * 20); /* 20 just incase parameters increase later */
	dx = dvector(1, tiePoints->npts);
	dy = dvector(1, tiePoints->npts);
	dxOrig = dvector(1, tiePoints->npts);
	dyOrig = dvector(1, tiePoints->npts);
	sigX = dvector(1, tiePoints->npts);
	sigY = dvector(1, tiePoints->npts);
	for (i = 1; i <= tiePoints->npts; i++)
	{
		sigX[i] = 1.0;
		sigY[i] = 1.0;
	}
	/*
	  Load data for fit
	*/
	xy = (svdData *)malloc(sizeof(svdData) * (tiePoints->npts + 1));
	xyOrig = (svdData *)malloc(sizeof(svdData) * (tiePoints->npts + 1));
	fprintf(stderr, "number of tie points %i\n", tiePoints->npts);
	nptsOrig = tiePoints->npts;
	for (i = 1; i <= tiePoints->npts; i++)
	{
		dxOrig[i] = tiePoints->offXT[i - 1];
		dyOrig[i] = tiePoints->offYT[i - 1];
		dx[i] = tiePoints->offXT[i - 1];
		dy[i] = tiePoints->offYT[i - 1];
		/* scale to km to help condition fit */
		xy[i].x = tiePoints->x[i - 1] * MTOKM - (matches->x0 + 0.5 * matches->nx * matches->stepX * matches->dx) * MTOKM;
		xy[i].y = tiePoints->y[i - 1] * MTOKM - (matches->y0 + 0.5 * matches->ny * matches->stepY * matches->dy) * MTOKM;
		xyOrig[i].x = tiePoints->x[i - 1] * MTOKM - (matches->x0 + 0.5 * matches->nx * matches->stepX * matches->dx) * MTOKM;
		xyOrig[i].y = tiePoints->y[i - 1] * MTOKM - (matches->y0 + 0.5 * matches->ny * matches->stepY * matches->dy) * MTOKM;
	}
	notDone = TRUE;
	nIterations = 0;
	nRejectForFit = 0;
	nReject = 0;
	// Loop to do fits
	while (notDone == TRUE)
	{ // Require at least 10 tiepoints
		if (tiePoints->npts <= 3)
			error("No solution, only %i tiepoints (must be > 4)\n", tiePoints->npts);
		// Compute width
		computePointsWidth(xy, tiePoints->npts, &xWidth, &yWidth);
		//	Setup fit.
		//	Added 07/18/24
		//	Downshift to constant fit if data don't support plane fit and
		//	moved inside loop to track points as they are culled.
		fprintf(stderr, "**** xWidth, yWidth %f %f\n", xWidth, yWidth);
		if (tiePoints->npts < 30 || xWidth < 75 || yWidth < 75)
		{
			fitDat->fitType = CONSTANTFIT;
			ma = 1;
			fitFunction = &constantFit;
			fprintf(stderr, "Downshifting to constant fit because too few points or the data too clustered\n");
			fprintf(stderr, "** Fitting for a constant ****\n");
		}
		else
		{
			ma = 3;
			fitFunction = &planeFit;
			fprintf(stderr, "** Fitting for a plane ****\n");
		}
		// Inside loop because ma can change
		u = dmatrix(1, tiePoints->npts, 1, ma);
		v = dmatrix(1, ma, 1, ma);
		CmatrixX = dmatrix(1, ma, 1, ma);
		CmatrixY = dmatrix(1, ma, 1, ma);
		w = dvector(1, ma);
		sigB = dvector(1, ma);
		//
		// Do the fits
		svdfit((void *)xy, dx, sigX, tiePoints->npts, aX, ma, u, v, w, &chisq, fitFunction);
		svdvar(v, ma, w, CmatrixX);
		svdfit((void *)xy, dy, sigY, tiePoints->npts, aY, ma, u, v, w, &chisq, fitFunction);
		svdvar(v, ma, w, CmatrixY);
		tieCountForFit = tiePoints->npts;
		nRejectForFit += nReject;
		// Print for debug, not on first fit covariance will be the same since errors are same for both x & y
		printFit(aX, aY, ma, stderr);
		printCovariance(CmatrixX, ma, stderr, "X");
		printCovariance(CmatrixY, ma, stderr, "Y");
		// Evaluate fit
		computeStats(fitDat, matches, xy, dx, dy, aX, aY,
					 &meanX, &meanY, &sigmaX, &sigmaY, &varX, &varY, tiePoints->npts, ma);
		fprintf(stderr, "nOrig %i nRject %i nReject total %i\n ", nptsOrig, nReject, nRejectForFit);
		fprintf(stderr, "Intermediate Res mean %f %f sigma %f %f\n\n", meanX, meanY, sigmaX, sigmaY);
		// Discard outliers and setup to redo if needed
		discardOutliers(fitDat, matches, tiePoints, xy, dx, dy, aX, aY, sigmaX, sigmaY, &nReject, ma);
		nIterations++;
		//
		// Quit if no more rejects or max iterations reached
		if (nReject == 0 || nIterations >= MAXTIEITERATIONS)
		{
			notDone = FALSE;
			break;
		}
		// Use residual sigmas, after scaling back to pixels for error estimate
		for (i = 1; i <= tiePoints->npts; i++)
		{
			// fprintf(stderr, "%f %f\n", dx[i], dy[i]);
			sigX[i] = sigmaX / (365.25 * matches->dx) * fitDat->deltaT;
			sigY[i] = sigmaY / (365.25 * matches->dy) * fitDat->deltaT;
		}
	} // Finished iterating
	//
	// Use non culled stats for error estimate (note will improve after masking)
	computeStats(fitDat, matches, xyOrig, dxOrig, dyOrig, aX, aY,
					 &meanX, &meanY, &sigmaX, &sigmaY, &varX, &varY, nptsOrig, ma);
	fprintf(stderr, "tiePoints->nPts %i nReject %i\n", tieCountForFit, nRejectForFit);
	fprintf(stderr, "Final Res mean %f %f sigmax %f %f\n", meanX, meanY, sigmaX, sigmaY);
	/*
	  write results
	*/
	printFitResult(fitDat, aX, aY, CmatrixX, CmatrixY, ma,
				   meanX, meanY, sigmaX, sigmaY, nReject, nRejectForFit, nptsOrig, tieCountForFit);
}

static void computeStats(lsFit *fitDat, matchResult *matches,
						 svdData *xy, double *dx, double *dy, double *aX, double *aY,
						 double *meanX, double *meanY, double *sigmaX, double *sigmaY, double *varX, double *varY,
						 int32_t nPts, int32_t ma)
{
	double px1, px2, px3, py1, py2, py3;
	double resX, resY;
	int32_t i;
	*meanX = 0.0;
	*meanY = 0.0;
	*varX = 0.0;
	*varY = 0.0;
	*sigmaX = 0.0;
	*sigmaY = 0.0;
	for (i = 1; i <= nPts; i++)
	{
		getCoeffs(aX, &px1, &px2, &px3, ma);
		resX = (dx[i] - (px1 + px2 * xy[i].x + px3 * xy[i].y)) * (matches->dx) / fitDat->deltaT * 365.25;
		getCoeffs(aY, &py1, &py2, &py3, ma);
		resY = (dy[i] - (py1 + py2 * xy[i].x + py3 * xy[i].y)) * (matches->dy) / fitDat->deltaT * 365.25;
		*meanX += resX;
		*varX += resX * resX;
		*meanY += resY;
		*varY += resY * resY;
	}
	*varX = *varX / ((double)nPts);
	*varY = *varY / ((double)nPts);
	*meanX = *meanX / ((double)nPts);
	*meanY = *meanY / ((double)nPts);
	*sigmaX = sqrt(*varX - *meanX * *meanX);
	*sigmaY = sqrt(*varY - *meanY * *meanY);
}

//	Ad hoc metric to make sure points cover a wide enough area for the fit
static void computePointsWidth(svdData *xy, int32_t nPts, double *xWidth, double *yWidth)
{

	int32_t i;
	double meanX, meanY;
	// Compute mean

	meanX = 0.0;
	meanY = 0.0;
	for (i = 1; i <= nPts; i++)
	{
		meanX += xy[i].x;
		meanY += xy[i].y;
	}
	meanX = meanX / nPts;
	meanY = meanY / nPts;
	// Determine max distance from center
	*xWidth = 0.0;
	*yWidth = 0.0;
	for (i = 1; i <= nPts; i++)
	{
		*xWidth = max(*xWidth, fabs(xy[i].x - meanX));
		*yWidth = max(*yWidth, fabs(xy[i].y - meanY));
	}
}

static void discardOutliers(lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints,
							svdData *xy, double *dx, double *dy, double *aX, double *aY,
							double sigmaX, double sigmaY, int32_t *nReject, int32_t ma)
{
	double px1, px2, px3, py1, py2, py3;
	double resX, resY;
	int32_t i;
	*nReject = 0;
	// fprintf(stderr, "%f %f %f %f\n", meanX, meanY, sigmaX, sigmaY);
	for (i = 1; i <= tiePoints->npts; i++)
	{
		// fprintf(stderr, "%i %i\n", i, *nReject);
		getCoeffs(aX, &px1, &px2, &px3, ma);
		resX = (dx[i] - (px1 + px2 * xy[i].x + px3 * xy[i].y)) * (matches->dx) / fitDat->deltaT * 365.25;
		getCoeffs(aY, &py1, &py2, &py3, ma);
		resY = (dy[i] - (py1 + py2 * xy[i].x + py3 * xy[i].y)) * (matches->dy) / fitDat->deltaT * 365.25;
		// fprintf(stderr, "%f %f  %f %f %i\n", fabs(resX), fabs(resY), sigmaX, sigmaY, REJECTTHRESHOLD);
		if (fabs(resX) > min(2 * sigmaX, REJECTTHRESHOLD) || fabs(resY) > min(2 * sigmaY, REJECTTHRESHOLD))
		{
			(*nReject)++;
		}
		else
		{
			/* keep only points that passed test in if statement */
			dx[i - *nReject] = dx[i];
			dy[i - *nReject] = dy[i];
			xy[i - *nReject].x = xy[i].x;
			xy[i - *nReject].y = xy[i].y;
		}
	} /* end for */
	tiePoints->npts -= *nReject;
}

// Set coefficients based on order of fit.
static void getCoeffs(double *aX, double *p1, double *p2, double *p3, int32_t ma)
{
	*p1 = aX[1];
	if (ma == 3)
	{
		*p2 = aX[2];
		*p3 = aX[3];
	}
	else if (ma == 1)
	{
		*p2 = 0.0;
		*p3 = 0.0;
	}
	else
		error("printFit: invalid number of params %i", ma);
}

static void printFit(double *aX, double *aY, int32_t ma, FILE *fpOut)
{
	double px1, py1, px2, py2, px3, py3;
	getCoeffs(aX, &px1, &px2, &px3, ma);
	fprintf(fpOut, "; dx fit, const, x coeff, y coeff\n");
	fprintf(fpOut, "Xfit =  %le %le %le\n", px1, px2 * MTOKM, px3 * MTOKM);
	getCoeffs(aY, &py1, &py2, &py3, ma);
	fprintf(fpOut, "; dy fit, const, x coeff, y coeff\n");
	fprintf(fpOut, "Yfit =  %le %le %le\n", py1, py2 * MTOKM, py3 * MTOKM);
}

static void printCovariance(double **Cmatrix, int32_t ma, FILE *fpOut, char *coord)
{
	fprintf(fpOut, "; %s_covariance_matrix\n", coord);
	if (ma == 3)
	{
		fprintf(fpOut, "C%s_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n", coord, 1, 1, 1,
				Cmatrix[1][1], Cmatrix[1][2] * MTOKM, Cmatrix[1][3] * MTOKM);
		fprintf(fpOut, "C%s_%1i1_%1i2_%1i3 =  %9.6le %9.6le %9.6le \n", coord, 2, 2, 2,
				Cmatrix[2][1] * MTOKM, Cmatrix[2][2] * MTOKM * MTOKM, Cmatrix[2][3] * MTOKM * MTOKM);
		fprintf(fpOut, "C%s_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n", coord, 3, 3, 3,
				Cmatrix[3][1] * MTOKM, Cmatrix[3][2] * MTOKM * MTOKM, Cmatrix[3][3] * MTOKM * MTOKM);
	}
	else if (ma == 1)
	{
		fprintf(fpOut, "C%s_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n", coord, 1, 1, 1,
				Cmatrix[1][1], 0., 0.);
		fprintf(fpOut, "C%s_%1i1_%1i2_%1i3 =  %9.6le %9.6le %9.6le \n", coord, 2, 2, 2,
				0., 0., 0.);
		fprintf(fpOut, "C%s_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n", coord, 3, 3, 3,
				0., 0., 0.);
	}
}

static void printFitResult(lsFit *fitDat, double *aX, double *aY, double **CmatrixX, double **CmatrixY, int32_t ma,
						   double meanX, double meanY, double sigmaX, double sigmaY,
						   int32_t nReject, int32_t nRejectForFit, int32_t nptsOrig, int32_t tieCountForFit)
{
	FILE *fpOut;

	fpOut = fopen(fitDat->fitFile, "w");
	fprintf(fpOut, "; Fit to data from match file(.dx,.dy)  %s\n", fitDat->matchFile);
	fprintf(fpOut, "; Tiefile %s \n", fitDat->tieFile);
	if (fitDat->fitType == PLANEFIT)
		fprintf(fpOut, "; Fit type used, x,y plane \n");
	else if (fitDat->fitType == CONSTANTFIT)
		fprintf(fpOut, "; Fit type used, constant \n");
	fprintf(fpOut, "; n rejected for fit, and n would reject on next iteration, percent rejected \n");
	fprintf(fpOut, "N_ties_rejected =  %i %i %i\n", nRejectForFit, nReject,
			(int32_t)(100.0 * (double)nRejectForFit / (double)(nptsOrig)));
	fprintf(fpOut, ";  Number of ties used in fit\n");
	fprintf(fpOut, "N_ties_used =  %i \n", tieCountForFit);
	printFit(aX, aY, ma, fpOut);
	fprintf(fpOut, "; X residual (all points), meanY, sigmaY (m/yr -- sigmaY (m) \n");
	fprintf(fpOut, "X_residual =  %10.2lf %10.2lf   %10.2lf\n", meanX, sigmaX, sigmaX / 365.25 * fitDat->deltaT);
	fprintf(fpOut, "; Y residual (all points), meanY, sigmaY (m/yr -- sigmaY (m) \n");
	fprintf(fpOut, "Y_residual =  %10.2lf %10.2lf  %10.2lf \n", meanY, sigmaY, sigmaY / 365.25 * fitDat->deltaT);
	printCovariance(CmatrixX, ma, fpOut, "X");
	printCovariance(CmatrixY, ma, fpOut, "Y");
	fprintf(fpOut, "; Y_covariance_matrix\n");
	fprintf(fpOut, "&");
	fclose(fpOut);
}

static void planeFit(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  Plane equation
	*/
	svdData *xx;
	double x1, y1;
	xx = (svdData *)x;
	afunc[1] = 1;
	afunc[2] = xx[i].x;
	afunc[3] = xx[i].y;
	return;
}

static void constantFit(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  Plane equation
	*/
	svdData *xx;
	double x1, y1;
	xx = (svdData *)x;
	afunc[1] = 1;
	return;
}
