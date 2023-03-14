#include "stdio.h"
#include "math.h"
#include "string.h"
#include <stdlib.h>
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "landsatSource64/clib/standard.h"
#include "landsatSource64/cRecipes/cRecipes.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "landsatSource64/cRecipes/nrutil.h"
#define MAXTIEPOINTS 500000
#define REJECTTHRESHOLD 50
#define MAXTIEITERATIONS 3
 
static void planeFit(void *x,int i,double *afunc, int ma);

typedef struct  {
        double x;
        double y;
} svdData;

/*
  Read tiepoint file for tiepoints.
*/
void computeLSFit( lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints)
{
	double **v,**u, **CmatrixX,**CmatrixY, chisq, *sigB,*sigX,*sigY, *w; /* Svd fit variabiles */ 
	double *aX,*aY; /* Solution for params */
	double *dx,*dy,*dxOrig,*dyOrig;   /* Offsets data moved to 1...n vector */
	double resX, resY;
	double varX,varY,sigmaX,sigmaY,meanX,meanY;
	char *fitFile;
	int nData, ma;
	int32 nptsOrig;
	double normX,normY;
	svdData *xy, *xyOrig;
	int32 nIterations,notDone; /* Used to control iterations loop */
	int32 nReject,nRejectForFit, tieCountForFit;
	FILE *fpOut;
	int32 i,j;

	if(fitDat->fitType==CONSTANTFIT) { ma=1; fprintf(stderr,"** Fitting for a constant ****\n");}
	else  if(fitDat->fitType==PLANEFIT) { ma=3; fprintf(stderr,"** Fitting for a plane ****\n");}
	else error("computeLSFit: invalid fit type \n");
	/* 
	   Init arrays for LS fit
	*/
	aX=(double *)malloc(sizeof(double)*20); /* 20 just incase parameters increase later */
	aY=(double *)malloc(sizeof(double)*20); /* 20 just incase parameters increase later */
	dx = dvector(1,tiePoints->npts);
	dy = dvector(1,tiePoints->npts);
	dxOrig = dvector(1,tiePoints->npts);
	dyOrig = dvector(1,tiePoints->npts);
	sigX = dvector(1,tiePoints->npts);	for(i=1; i<=tiePoints->npts; i++) sigX[i]=1.0;
	sigY = dvector(1,tiePoints->npts);	for(i=1; i<=tiePoints->npts; i++) sigY[i]=1.0;
	u = dmatrix(1, tiePoints->npts,1,ma);
	v = dmatrix(1,ma,1,ma);
	CmatrixX = dmatrix(1,ma,1,ma);
	CmatrixY = dmatrix(1,ma,1,ma);
	w = dvector(1,ma); 
	sigB=dvector(1,ma);
	/*
	  Load data for fit
	*/
	xy=(svdData *) malloc(sizeof(svdData)*(tiePoints->npts+1));
	xyOrig=(svdData *) malloc(sizeof(svdData)*(tiePoints->npts+1));
	fprintf(stderr,"number of tie points %i\n", tiePoints->npts);
	normX=matches->stepX*matches->dx *matches->nx;
	normY=matches->stepY*matches->dy *matches->ny;
	nptsOrig= tiePoints->npts;
	for(i=1; i <= tiePoints->npts; i++) {
		dxOrig[i] = tiePoints->offXT[i-1];
		dyOrig[i] = tiePoints->offYT[i-1];
		dx[i]=tiePoints->offXT[i-1];
		dy[i]=tiePoints->offYT[i-1];
		/* scale to km to help condition fit */
		xy[i].x=tiePoints->x[i-1] *MTOKM         -(matches->x0 + 0.5 * matches->nx*matches->stepX * matches->dx)*MTOKM;
		xy[i].y=tiePoints->y[i-1] * MTOKM        -(matches->y0 + 0.5 * matches->ny*matches->stepY * matches->dy)*MTOKM;
		xyOrig[i].x=tiePoints->x[i-1] *MTOKM  -(matches->x0 + 0.5 * matches->nx*matches->stepX * matches->dx)*MTOKM;
		xyOrig[i].y=tiePoints->y[i-1] * MTOKM -(matches->y0 + 0.5 * matches->ny*matches->stepY * matches->dy)*MTOKM;
		/*error("%f %f \n",-matches->x0,-matches->y0);*/
		/*	if(i % 1 == 0) fprintf(stdout," %lf %lf %lf %lf\n",	xy[i].x,	xy[i].y,dx[i],dy[i]);*/
	}
	/*
	  do fits 
	*/
	
	notDone=TRUE;
	nIterations =0; 
	nRejectForFit=0;
	nReject=0;
	while(notDone == TRUE && nIterations < MAXTIEITERATIONS) {
		if(fitDat->fitType == PLANEFIT) {
			svdfit((void *)xy,dx,sigX,tiePoints->npts,aX,ma,u,v,w, &chisq,&planeFit);
			svdvar(v,ma,w,CmatrixX);
			fprintf(stderr,"Xoff fit %le %le %le\n",aX[1],aX[2]*MTOKM,aX[3]*MTOKM);		
			svdfit((void *)xy,dy,sigY,tiePoints->npts,aY,ma,u,v,w, &chisq,&planeFit);
			fprintf(stderr,"Yoff fit %le %le %le\n",aY[1],aY[2]*MTOKM,aY[3]*MTOKM);
			svdvar(v,ma,w,CmatrixY);
			fprintf(stderr,"CX_%1i1_%1i2_%1i3 %9.6le %9.6le %9.6le \n",1,1,1,CmatrixX[1][1],		  CmatrixX[1][2]*MTOKM,	       CmatrixX[1][3]*MTOKM);
       			fprintf(stderr,"CX_%1i1_%1i2_%1i3 %9.6le %9.6le %9.6le \n",2,2,2,CmatrixX[2][1]*MTOKM, CmatrixX[2][2]*MTOKM*MTOKM,CmatrixX[2][3]*MTOKM*MTOKM);
       			fprintf(stderr,"CX_%1i1_%1i2_%1i3 %9.6le %9.6le %9.6le \n",3,3,3,CmatrixX[3][1]*MTOKM, CmatrixX[3][2]*MTOKM*MTOKM,CmatrixX[3][3]*MTOKM*MTOKM);
			fprintf(stderr,"CY_%1i1_%1i2_%1i3 %9.6le %9.6le %9.6le \n",1,1,1,CmatrixY[1][1], 		  CmatrixY[1][2]*MTOKM,		CmatrixY[1][3]*MTOKM);
       			fprintf(stderr,"CY_%1i1_%1i2_%1i3 %9.6le %9.6le %9.6le \n",2,2,2,CmatrixY[2][1]*MTOKM, CmatrixY[2][2]*MTOKM*MTOKM,CmatrixY[2][3]*MTOKM*MTOKM);
       			fprintf(stderr,"CY_%1i1_%1i2_%1i3 %9.6le %9.6le %9.6le \n",3,3,3,CmatrixY[3][1]*MTOKM, CmatrixY[3][2]*MTOKM*MTOKM,CmatrixY[3][3]*MTOKM*MTOKM);
		} else error("computeLSFit: incorrect fit type \n");
		/*
		  Evaluate fit
		*/	
		meanX=0.0; meanY=0.0;	varX=0.0; varY=0.0;	sigmaX=0.0; sigmaY=0.0;
		/*
		  Compute residual in m/yr
		 */
		 
		for(i=1; i <= tiePoints->npts; i++) {
			resX =(dx[i]-(aX[1] +aX[2]*xy[i].x + aX[3]*xy[i].y)) * (matches->dx) /fitDat->deltaT * 365.25;
			resY =(dy[i]-(aY[1] +aY[2]*xy[i].x +  aY[3]*xy[i].y)) * (matches->dy) /fitDat->deltaT * 365.25;
			meanX += resX; varX += resX * resX;
			meanY += resY; varY += resY * resY;
		}
		sigmaX = (varX- meanX*meanX) / ((double)(tiePoints->npts - 1));
		sigmaX=sqrt(sigmaX);
		sigmaY = (varY- meanY*meanY) / ((double)(tiePoints->npts - 1));
		sigmaY=sqrt(sigmaY);
		varX=sqrt(varX/((double)tiePoints->npts));
		varY=sqrt(varY/((double)tiePoints->npts));
		meanX=meanX/((double)tiePoints->npts);
		meanY=meanY/((double)tiePoints->npts);
		/* Use residual, after scaling back to pixels for error estimate */
		for(i=1; i<=tiePoints->npts; i++) { sigX[i]=sigmaX/(365.25*matches->dx) *  fitDat->deltaT ; sigY[i]=sigmaY/(365.25*matches->dy) *  fitDat->deltaT;}
		/*
		  Toss outliers and set up do do over
		*/
		nRejectForFit+=nReject; tieCountForFit=	tiePoints->npts; /* save for output stats, so as not to count the last iteration */
		nReject=0;
		fprintf(stderr,"%f %f %f %f\n",meanX, meanY, sigmaX, sigmaY);
		for(i=1; i <= tiePoints->npts; i++) {
			resX =(dx[i]-(aX[1] +aX[2]*xy[i].x + aX[3]*xy[i].y)) * (matches->dx) /fitDat->deltaT * 365.25;
			resY =(dy[i]-(aY[1] +aY[2]*xy[i].x +  aY[3]*xy[i].y)) * (matches->dy) /fitDat->deltaT * 365.25;
			if(fabs(resX) > min(3*sigmaX ,REJECTTHRESHOLD) || fabs(resY) > min(3*sigmaY ,REJECTTHRESHOLD) ) {
				nReject++;
			} else {
				/* keep only points that passed test in if statement */
				dx[i-nReject]=dx[i];
				dy[i-nReject]=dy[i];
				xy[i-nReject].x=xy[i].x;
				xy[i-nReject].y=   xy[i].y;
			}
		} /* end for */
		tiePoints->npts -= nReject;
		fprintf(stderr, "tiePoints->nPts %i\n", tiePoints->npts);
		nIterations++;
	}
	/*
	  Evaluate final fit, but keep outliers
	*/	
	meanX=0.0; meanY=0.0;		varX=0.0; varY=0.0;		sigmaX=0.0; sigmaY=0.0;
	for(i=1; i <= nptsOrig; i++) {
		resX =(dxOrig[i]-(aX[1] +aX[2]*xyOrig[i].x + aX[3]*xyOrig[i].y)) * (matches->dx) /fitDat->deltaT * 365.25;
		resY =(dyOrig[i]-(aY[1] +aY[2]*xyOrig[i].x +  aY[3]*xyOrig[i].y)) * (matches->dy) /fitDat->deltaT * 365.25;
		meanX += resX; varX += resX * resX;
		meanY += resY; varY += resY * resY;
	}
	fprintf(stderr,"%lf %lf %i\n",meanX,  varX,nptsOrig);
	meanX=meanX/(double)(nptsOrig );
	meanY=meanY/(double)(nptsOrig );
	varX=varX/(double)(nptsOrig - 1);
	varY=varY/(double)(nptsOrig - 1);
	sigmaX = (varX- meanX*meanX); 		sigmaX=sqrt(sigmaX);
	sigmaY = (varY- meanY*meanY) ;		sigmaY=sqrt(sigmaY);
	/*
	  write results
	*/
	/*fitFile = (char *)malloc(strlen(fitDat->matchFile)+6); 	fitFile[0]='\0';
	fitFile=strcpy(fitFile,fitDat->matchFile);
	fitFile=strcat(fitFile,".fit");*/
	fpOut=fopen(fitDat->fitFile,"w");
	fprintf(fpOut,"; Fit to data from match file(.dx,.dy)  %s\n",fitDat->matchFile);
	fprintf(fpOut,"; Tiefile %s \n",fitDat->tieFile);
	if(fitDat->fitType == PLANEFIT) fprintf(fpOut,"; Fit type used, x,y plane \n");
	fprintf(fpOut,"; n rejected for fit, and n would reject on next iteration, percent rejected \n");
	fprintf(fpOut,"N_ties_rejected =  %i %i %i \n",nRejectForFit, nReject, (int32)(100.0*(double)nRejectForFit/(double)(nptsOrig)));

	fprintf(fpOut,";  Number of ties used in fit\n");
	fprintf(fpOut,"N_ties_used =  %i \n",tieCountForFit);

	fprintf(fpOut,"; dx fit, const, x coeff, y coeff\n");
	fprintf(fpOut,"Xfit =  %le %le %le\n",aX[1],aX[2]*MTOKM,aX[3]*MTOKM);		

	fprintf(fpOut,"; dy fit, const, x coeff, y coeff\n");
	fprintf(fpOut,"Yfit =  %le %le %le\n",aY[1],aY[2]*MTOKM,aY[3]*MTOKM);	

	fprintf(fpOut,"; X residual (all points) meanX, sigmaX (m/yr) - sigmaX (m)\n");
	fprintf(fpOut,"X_residual =  %10.2lf %10.2lf   %10.2lf\n",  meanX,sigmaX,sigmaX/365.25 * fitDat->deltaT  );

	fprintf(fpOut,"; Y residual (all points), meanY, sigmaY (m/yr -- sigmaY (m) \n");
	fprintf(fpOut,"Y_residual =  %10.2lf %10.2lf  %10.2lf \n",  meanY,sigmaY,sigmaY/365.25 * fitDat->deltaT  );
	fprintf(fpOut,"; X_covariance_matrix\n");

	fprintf(fpOut,"CX_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n",1,1,1,CmatrixX[1][1],		  CmatrixX[1][2]*MTOKM,	       CmatrixX[1][3]*MTOKM);
	fprintf(fpOut,"CX_%1i1_%1i2_%1i3 =  %9.6le %9.6le %9.6le \n",2,2,2,CmatrixX[2][1]*MTOKM, CmatrixX[2][2]*MTOKM*MTOKM,CmatrixX[2][3]*MTOKM*MTOKM);
	fprintf(fpOut,"CX_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n",3,3,3,CmatrixX[3][1]*MTOKM, CmatrixX[3][2]*MTOKM*MTOKM,CmatrixX[3][3]*MTOKM*MTOKM);
	fprintf(fpOut,"; Y_covariance_matrix\n");
	fprintf(fpOut,"CY_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n",1,1,1,CmatrixY[1][1], 		  CmatrixY[1][2]*MTOKM,		CmatrixY[1][3]*MTOKM);
	fprintf(fpOut,"CY_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n",2,2,2,CmatrixY[2][1]*MTOKM, CmatrixY[2][2]*MTOKM*MTOKM,CmatrixY[2][3]*MTOKM*MTOKM);
	fprintf(fpOut,"CY_%1i1_%1i2_%1i3 = %9.6le %9.6le %9.6le \n",3,3,3,CmatrixY[3][1]*MTOKM, CmatrixY[3][2]*MTOKM*MTOKM,CmatrixY[3][3]*MTOKM*MTOKM);
	fprintf(fpOut,"&");
	fclose(fpOut);
}


static  void planeFit(void *x,int i,double *afunc, int ma) 
{
	/*
	  Plane equation
	*/
	svdData *xx;
	double x1,y1;
	xx = (svdData *)x;
	/*   fprintf(stderr,"%i %f\n",i,xy[0].x);*/
	afunc[1] = 1;
	afunc[2] = xx[i].x;
	afunc[3] = xx[i].y;
	/*	if(i % 100 == 0) fprintf(stderr,"%i %lf %lf %lf\n",afunc[1],afunc[2],afunc[3]);*/
	return;
}
