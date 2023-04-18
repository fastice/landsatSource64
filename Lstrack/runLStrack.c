#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "fftw3.h"
#include"string.h"
#include "clib/standard.h"
#include "time.h"
#include "lstrack.h"

/*
  Helper functions
*/
static  void loadTiff(char *file,LSimage *Image,TIFF *tif,GTIF *gtif); /* Read tiff files in and apply high-pass filter */
static void writeTiffAsRaw(LSimage *Image,char *file);  /* Used for debugging - dump image buffers */
static void mallocMatch(matchResult *matches) ;    /* Memory allocation  */
static int32 getLatePosition(LSimage *earlyImage,LSimage *lateImage,matchParams *matchP,
	uint32 ie,uint32 je,uint32 *il,uint32 *jl,int32 *chipSize,int32 *lx1, int32 *lx2,int32 *ly1,int32 *ly2) ;
static float fastCorr(LSimage *earlyImage,LSimage *lateImage,int32 chipSize,int32 ie,int32 je, int32 il, 
	int32 jl, int32 lx1, int32 lx2,int32 ly1,int32 ly2,float *offX,float *offY);
static unsigned char isMasked(double x, double y, matchParams *matchP);
static unsigned char velGuess(double x, double y, matchParams *matchP, double dx,double dy, int32  *iVel, int32 *jVel,
	 int32 *slow, int32_t *highStrain) ;
static float **lsmatrix(int32 nrl, int32 nrh, int32 ncl, int32 nch)	;
static float overSamplePeakCorr(int32 iMax, int32 jMax, float *offX, float *offY);
static void parseLSDate(char *file,LSimage *Image);
static long julday(int32_t mm, int32_t id, int32_t iyyy);

/*
  Globals for FFTW and corelation
*/
float **corr;
fftwf_complex *corrPatchSpace, *corrPatchFreq, *corrPatchOverFreq, *corrPatchOverSpace;
fftwf_plan forwardCorrPlan, inverseCorrPlan;

/************************************ runLStrack ******************************************
Main routine to do LANDSAT speckle tracking
************************************ runLStrack ******************************************/
matchResult runLStrack(char *earlyFile, char *lateFile, matchParams *matchP )
{
	extern float **corr;
	extern fftwf_complex *corrPatchSpace, *corrPatchFreq, *corrPatchOverFreq,*corrPatchOverSpace;
	extern fftwf_plan forwardCorrPlan, inverseCorrPlan;
	extern int32 nMatch, nAttempt, nTotal;

	TIFF *tifEarly=(TIFF*)0, *tifLate=(TIFF*)0 ; 	 /* TIFF-level descriptor */
	GTIF *gtifEarly=(GTIF*)0, *gtifLate=(GTIF*)0;	 /* GeoKey-level descriptor */
	LSimage earlyImage, lateImage; 			/* tiff images after being read in */
	matchResult matches;						/* match result */
	uint32 idum,idum1;					        /* Dummy variables */
	void *pdum;	      							/* Dummy variables */
	int32 im,jm,ie,je,il,jl,i,j;								/* LCV - counters for main loops */
	double minX, minY;						/* Temp variabiles used to compute bounds */
	double maxX, maxY;						/* Temp variabiles used to compute bounds */
	double nX, nY;							/* Temp variabiles used to compute bounds */
	double x,y;								/* Current position being tracked */
	int32 slow, highStrain; 								/* Set when velocity map indicates chip can be increased and search region reduced */
	int32 chipSize;							/* Chipsize */
	int32 lx1, ly1, lx2, ly2, k, k1, iVel, jVel; 					/* Search region  */   
	float offX, offY;							/* Returned offsets - i.e, final result */
	uint32 edgeBuf;							/* "margin" around edges - i.e., don't use data within edgebuf of image edges */
	int32_t width, widthO;							/* wdith of corrPatch, and oversampled version corrPatchOver */
	float pixelCorr;							/* Returned correlation value */
	time_t currTime, startTime;					/* Used for timer */
	FILE *fpDebug, *fpDebug1, *fpDebug2;     		 /* Debug variabiles */
	float **tmp, *tmp1, **tmpa, *tmp1a, fdum[1];      /* Debug variabiles */
	FILE *fp;
	uint32 countS,countM,countF; 
	/*debug 
	  tmp=(float **)malloc(sizeof(float*)*OSSIZE*2*NOVER);     tmp1=(float *)malloc(sizeof(float)*OSSIZE*2*NOVER*OSSIZE*2*NOVER);
	  for(i=0; i < OSSIZE*2*NOVER; i++) tmp[i]=&(tmp1[i*OSSIZE*2*NOVER]);
	  tmpa=(float **)malloc(sizeof(float*)*OSSIZE*2);   tmp1a=(float *)malloc(sizeof(float)*OSSIZE*2*OSSIZE*2);
	  for(i=0; i < OSSIZE*2; i++) tmpa[i]=&(tmp1a[i*OSSIZE*2]);
	  end debug */
	/*
	  Hardwire parameters now - change later ***
	*/
	countS=0;
	countM=0;
	countF=0;
	matchP->deltaX=100;      matchP->deltaY=100; /* In pixels */
	matchP->nx=0;	   	     matchP->ny=0;
	nMatch=0; nAttempt=0; nTotal=0; /* Zero counters */
	edgeBuf=100;
	/* 
	   Open TIFF and GTIFF descriptor to read GeoTIFF tags 
	*/
	startTime=time(NULL);
	fprintf(stderr,"Early/latefile file %s %s\n",earlyFile,lateFile);
	tifEarly=XTIFFOpen(earlyFile,"r"); if (!tifEarly) error("couldn't open tifEarly");
	tifLate =XTIFFOpen(lateFile,"r");  if (!tifLate) error("couldn't open tifLate");
	gtifEarly = GTIFNew(tifEarly); if (!gtifEarly) error("couldn't open geo tif Early");
	gtifLate = GTIFNew(tifLate); if (!gtifLate) error("couldn't open geo tif Late");
	/*
	  Read in images
	*/
	loadTiff(earlyFile,&earlyImage, tifEarly,gtifEarly); matches.jdEarly=earlyImage.jd;
	loadTiff(lateFile,&lateImage, tifLate,gtifLate); matches.jdLate=lateImage.jd;
	matchP->deltaT=lateImage.jd - earlyImage.jd;
	fprintf(stderr, "Dates and deltaT %f %f %f\n", earlyImage.jd, lateImage.jd, matchP->deltaT);
	/*
	  Check projection
	*/
	if(matchP->proj == NSIDCSOUTH || matchP->proj == NSIDCNORTH) {
		/* Early proj */
		if(earlyImage.proj < 0) earlyImage.proj = matchP->proj; 
		else if(earlyImage.proj != matchP->proj) 
			error("cmd line proj %i inconsistent with %s %i\n",matchP->proj,earlyFile,earlyImage.proj);
		/* late proj*/
		if(lateImage.proj < 0) lateImage.proj = matchP->proj; 
		else if(lateImage.proj != matchP->proj) 
			error("cmd line proj %i inconsistent with %s %i\n",matchP->proj,lateFile,lateImage.proj);
    } else {
			fprintf(stderr,"EPSG 1/2 %i %i\n", earlyImage.proj, lateImage.proj);
			if(earlyImage.proj == NSIDCSOUTH || earlyImage.proj == NSIDCNORTH) matchP->proj = earlyImage.proj;
            else if(lateImage.proj == NSIDCSOUTH || lateImage.proj == NSIDCNORTH) matchP->proj = lateImage.proj;	
	}
 	if(matchP->proj !=NSIDCSOUTH && matchP->proj != NSIDCNORTH && matchP->proj != NULLPROJ) error("Invalid Projection either not set in tiff or cmd line %i \n",matchP->proj);
	fprintf(stderr,"Match EPSG %i\n",matchP->proj);
	fprintf(stderr,"Time To Load Image %f \n",(float)(time(NULL)-startTime)  ); currTime=time(NULL);
	/*
	  Match coordinate range
	*/
	minX=max(earlyImage.x0,lateImage.x0);       minY=max(earlyImage.y0,lateImage.y0);
	/* Ensure origin for stepping is integer number of pixels in first image */
	if(minX != earlyImage.x0) {
		idum=(uint32)((minX-earlyImage.x0)/earlyImage.dx); minX=earlyImage.x0 + idum * earlyImage.dx;	
	}
	if(minY != earlyImage.y0) {
		idum=(uint32)((minY-earlyImage.y0)/earlyImage.dy); minY=earlyImage.y0 + idum * earlyImage.dy;
	}
	matches.x0=minX+matchP->deltaX * earlyImage.dx;  
	matches.y0=minY+matchP->deltaY * earlyImage.dy; 
	matches.dx=earlyImage.dx; matches.dy=earlyImage.dy;
	matches.stepX=matchP->stepX; 	matches.stepY=matchP->stepY;
    if( fabs(earlyImage.dx-lateImage.dx) > 0.0001  || fabs(earlyImage.dy-lateImage.dy) > 0.0001) {
		error ("Image pixel sizes differ %10.5f %10.5f %10.5f %10.5\n", earlyImage.dx, earlyImage.dy,
				lateImage.dx,lateImage.dy);
	}
	/*
	  Compute max dimensions
	*/
	maxX=min(earlyImage.x0 + earlyImage.nx*earlyImage.dx, lateImage.x0 + lateImage.nx*lateImage.dx);
	maxY=min(earlyImage.y0 + earlyImage.ny*earlyImage.dy, lateImage.y0 + lateImage.ny*lateImage.dy);
	/*
	  Compute total number of matches to do in x and y directions
	*/
	nX=(uint32)((maxX-matches.x0)/earlyImage.dx) -e dgeBuf; 
	nX=nX/matchP->stepX; /*  pixels, then steps  */
	nY=(uint32)((maxY-matches.y0)/earlyImage.dy) - edgeBuf; 
	nY=nY/matchP->stepY;
	fprintf(stderr,"nX,nY %lf %lf  nx,ny %i %i	\n",nX,nY,matchP->nx,matchP->ny);
	
	if( (matchP->nx) <=0) {	
		matches.nx=(uint32)nX;    
	} else {   	
		matches.nx=min(nX,matches.nx); 
	}/* Use specified values if less than calcuated */

	if((matchP->ny) <=0) {	
		matches.ny=(uint32)nY;
	} else {
		matches.ny=min(nY,matches.ny); 
	
	fprintf(stderr,"%lf %lf %lf %lf \n \n",minX,maxX,minY,maxY);
	fprintf(stderr,"nX,nY %lf %lf  nx,ny %i %i	\n",nX,nY,matchP->nx,matchP->ny);
	fprintf(stderr,"%f %f\n",(maxX-matches.x0),((maxX-matches.x0)/earlyImage.dx) );
	fprintf(stderr,"match start %lf %lf number of matches in x/y %i %i %i %i\n",matches.x0,matches.y0,matches.nx,matches.ny,matchP->stepX,matchP->stepY);
	/* 
	   Do final mallocs
	*/
	mallocMatch(&matches);   
	/* debug fpDebug=fopen("testcorr","w"); fpDebug1=fopen("testcorr1","w");  fpDebug2=fopen("testcorr2","w"); end debug */
	/* 
	   Loop do correlation (i,j) index match output
	*/
	ie=(uint32)( (matches.y0-earlyImage.y0)/earlyImage.dy + 0.5);	
	for(im=0; im < matches.ny; im++) {
		je = (uint32)( (matches.x0-earlyImage.x0)/earlyImage.dx + 0.5); /* Compute first pixel - note j=x, i=y, always index i,j */
		y = earlyImage.y0 + ie*earlyImage.dy;
		for(jm=0; jm < matches.nx; jm++) {
			matches.Rho[im][jm]=-1.0;    matches.X[im][jm]=NODATA;	matches.Y[im][jm]=NODATA; matches.type[im][jm]=0;
			nTotal++;
			x=earlyImage.x0+je*earlyImage.dx;
			slow=0;
			/*
			  Is masked?
			*/
			if(isMasked(x,y,matchP) == TRUE) { matches.type[im][jm] |= ISMASKED; goto skipmatch;}
			/*
			  Get postion in late image. Skip if not found
			*/
			if(getLatePosition(&earlyImage,&lateImage,matchP,ie, je,&il,&jl,&chipSize,&lx1,&lx2,&ly1,&ly2)  < 0)  { matches.type[im][jm] |= NODATATOUSE; goto skipmatch;}
			/* 
			   velocity guess
			*/
			velGuess( x, y,matchP, earlyImage.dx,earlyImage.dy,&iVel, &jVel,&slow,&highStrain) ;

			il+=iVel; jl+=jVel;
			if(highStrain == TRUE) matches.type[im][jm] |= HIGHSTRAIN;
			if(slow > 0 ) {     /* Slow area with little variation, so increase chip and reduce search */
				if(slow == 1) {
					chipSize *= 1.5;
					if(chipSize % 2 == 0) chipSize++;
					lx1=-MINSEARCH; ly1=-MINSEARCH;
					lx2=MINSEARCH; ly2=MINSEARCH;
					matches.type[im][jm] |= SLOWVEL;
					countS++;
				 } else {
				 /* This is a bit of an afterthought, and it allows mid-ranage speeds (100,1000) to use a smaller chip size */
					lx1=-MINSEARCH*2; ly1=-MINSEARCH*2;
					lx2=MINSEARCH*2; ly2=MINSEARCH*2;				 
					matches.type[im][jm] |= FASTVEL;
					countM++;
				 }
			} else {	
				matches.type[im][jm] |= FASTVEL; 
				countF++; 
		    }
			/*
			  Run correlator.
			*/
		  	/*fprintf(stderr,"%i %i -- %i %i --- %i %i %i %i\n",ie,je,il,jl,lx1,lx2,ly1,ly2);*/
			pixelCorr=fastCorr(&earlyImage,&lateImage,chipSize,ie,je,il,jl, lx1, lx2,ly1,ly2,&offX,&offY) ;
			if(pixelCorr > NOUSEABLEDATA) { nAttempt++; }  /*  means no correlation tried because there was no data to correlate */
			if(pixelCorr==NOCORR) {matches.type[im][jm] |= LOWCORR; goto skipmatch;}
			if(pixelCorr==NOUSEABLEDATA) {matches.type[im][jm] |= NODATATOUSE; goto skipmatch;}
			/*	fprintf(stderr,"%li %li %li %li\n",ie,je,il,jl);*/
			matches.X[im][jm]=offX+(float)jVel;
			matches.Y[im][jm]=offY +(float)iVel;
			matches.Rho[im][jm]= pixelCorr;
			nMatch++;
		skipmatch:
			je += matchP->stepX;
		}
		ie += matchP->stepY;
		if( im % 20 ==  0 && im > 1) { fprintf(stderr,"%i %i %i %i %i %i %% Match %8.1f %% nMatch %i nAttempt %i nPts %i  ",ie,je,im,jm,ly1,ly2, ((float)nMatch/(float)max(nAttempt,1))*100,nMatch,nAttempt,nTotal);
			fprintf(stderr,"Time for this batch %8.1f (s) Total Elapsed  %8.2f  (m)\n",
				(float)(time(NULL)-currTime), (float)(time(NULL)-startTime)/60.); currTime=time(NULL);}
	}
	/*	writeTiffAsRaw(&earlyImage,"test1");
		writeTiffAsRaw(&lateImage,"test2");*/
	fprintf(stderr,"count Slow,Med,Fast  %i %i %i\n",countS,countM,countF);
	return(matches);
}

/******************************* isMasked  ***************************************
Checks if point is masked
*************************************** ****************************************/	
static unsigned char isMasked(double x, double y, matchParams *matchP) 
{
	int32 im,jm;
	/*
	  no mask
	*/
	if(matchP->maskFile == NULL) return(FALSE);  /* no mask, so not masked */

	jm=(int32)((x-matchP->mask.x0)/matchP->mask.dx + 0.5);
	im=(int32)((y-matchP->mask.y0)/matchP->mask.dy + 0.5);
	if(jm < 0 || im < 0 || jm >= matchP->mask.nx || im >= matchP->mask.ny ) return(FALSE); /* outside of bounds -> no masked */
	
	if(matchP->mask.m[im][jm] > 0) return(FALSE); /* Positive mask value indicates not masked */

	return(TRUE);
}

/******************************* velGuess  ***************************************
Guess at offset from velocity
*************************************** ****************************************/	
static unsigned char velGuess(double x, double y, matchParams *matchP, double dx,double dy, int32  *iVel, int32 *jVel,int32 *slow,int32_t *highStrain) 
{
	int32 im,jm,i,j;
	double vxPt,vyPt, t,u,xi,yi;
	double p1,p2,p3,p4, vsq;
	double mean,sigma,var;
	int32 n;
	/*
	  Default slow flag
	*/
	*slow=FALSE;
	*highStrain=FALSE;
	/* 
	   Default to zero so if no data, no harm
	*/
	*iVel=0;
	*jVel=0;
	/*
	  no velocity
	*/
	if(matchP->velFile == NULL) return(FALSE);  /* no velocity, so return */
	/* PS to velocity coordinates */
	xi=((x-matchP->velocity.x0)/matchP->velocity.dx + 0.5);
	yi=((y-matchP->velocity.y0)/matchP->velocity.dy + 0.5);
	/* integer coords */
	jm=(int32)xi;
	im=(int32)yi;
	/*	fprintf(stderr,"%lf %lf %lf %lf %i %i %f %f %f\n",x,y,xi,yi,jm,im,matchP->velocity.x0,matchP->velocity.y0,matchP->velocity.dx); */

	if(jm < 0 || im < 0 || jm >= matchP->velocity.nx || im >= matchP->velocity.ny ) return(FALSE); /* outside of bounds -> no velocity */
	
	if(matchP->velocity.vx[im][jm]  <(NODATA+1) || matchP->velocity.vy[im][jm]  <(NODATA+1)) return(FALSE); /* no velocity value */

	t = (float)(xi - (double)jm);
	u = (float)(yi - (double)im);
	/* 
	   Interpolate vx
	*/
	p1 = matchP->velocity.vx[im][jm];              p2 = matchP->velocity.vx[im][jm+1];
	p3 = matchP->velocity.vx[im+1][jm+1];    p4 = matchP->velocity.vx[im+1][jm];
	if(p1 <= (NODATA+1) || p2 <= (NODATA+1) || p3 <= (NODATA+1) || p4 <= (NODATA+1)) vxPt=max(max(max(p1,p2),p3),p4);
	else vxPt = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 +    (1.0 - t) * u* p4);
	/*	fprintf(stderr,"%lf %f %i %i\ %lf %lf \n",(double)(vxPt/365.), matchP->velocity.vx[im][jm],im,jm,t,u);*/
	/* projected delta in pixels */
       	*jVel=(int32)( (vxPt/365. ) *matchP->deltaT/dx);
	/* 
	   Interpolate vy
	*/
	p1 = matchP->velocity.vy[im][jm];              p2 = matchP->velocity.vy[im][jm+1];
	p3 = matchP->velocity.vy[im+1][jm+1];    p4 = matchP->velocity.vy[im+1][jm];
	if(p1 <= (NODATA+1) || p2 <= (NODATA+1) || p3 <= (NODATA+1) || p4 <= (NODATA+1)) vyPt=max(max(max(p1,p2),p3),p4);
	else vyPt = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 +     (1.0 - t) * u* p4);
	/* projected delta in pixels */
       	*iVel=(int32)( (vyPt/365. ) *matchP->deltaT/dy);   
	/*
	  determin whether to set flow flag
	  set slow flag if a) mean is less than SLOWTHRESH, b) sigma of speed for region around point is < 50 m/yr
	*/
	mean=0.0;
	var=0.0;
	n=2;
	if(matchP->slowFlag == TRUE) {
		for(i=im-n; i <= im+n; i++) {
			for(j=jm-n; j <= jm+n; j++) {
				if(matchP->velocity.vx[i][j] < (NODATA+1) || matchP->velocity.vy[i][j] < (NODATA+1)) return(TRUE); /* Don't bother to near a hole */
				vsq=matchP->velocity.vx[i][j] * matchP->velocity.vx[i][j]  + matchP->velocity.vy[i][j] * matchP->velocity.vy[i][j] ;
				var+=vsq;
				mean+=sqrt(vsq);
			}
		}
		var=var/((2*n+1)*(2*n+1));
		mean=mean/((2*n+1)*(2*n+1));
		sigma=sqrt(var - mean*mean);
		if(sigma < SIGMATHRESH && mean < SLOWTHRESH) *slow=1;
		else if(sigma < 4*SIGMATHRESH &&  mean < MIDTHRESH) *slow=2;
		if(sigma > SIGMATHRESH) *highStrain=TRUE;
	}
	return(TRUE);
}

/******************************* fastCorr  ***************************************
do correlation in space domain on integer pixels, then oversample peak by factor NOVER
returns peak correlation value and *offX and *offY  
*************************************** ****************************************/	
static float fastCorr(LSimage *earlyImage, LSimage *lateImage, int32 chipSize, int32 ie, int32 je, 
					  int32 il, int32 jl, int32 lx1, int32 lx2,int32 ly1, int32 ly2, float *offX, float *offY)	
{
	int32 i,j; /* Outer loop to move search window */
	int32 ii,jj,iiMin,iiMax,jjMin,jjMax;
	double meanL,meanE,varE,varL,Xcorr,sigmaE2,sigmaL2,sigmaMaxE,sigmaMaxL;
	int32 lx1a,lx2a,ly1a,ly2a;
	float maxCorr,maxCorrHi;
	double tmpCorr;
	int32 iMax,jMax,iMax1,jMax1;
	int32 corrReturn;
	int32 width, widthO;
	double npts;
	extern float **corr;
	*offX=NODATA; *offY=NODATA;
	iiMin=-(chipSize-1)/2;  iiMax=(chipSize-1)/2; jjMin=-(chipSize-1)/2;   jjMax=(chipSize-1)/2;
	/*
	 Check bounds and return no data if out of bounds
	*/
        if ( (il+iiMax) > lateImage->ny || (ie+iiMax) > earlyImage->ny || (jl+jjMax) > lateImage->nx || (je+jjMax) > earlyImage->nx)  { return(NOUSEABLEDATA);}
        if ( (il+iiMin) < 0		       || (ie+iiMin) < 0 			|| (jl+jjMin) < 0 		     || (je+jjMin) < 0)  			{ return(NOUSEABLEDATA);}
	/* 
		Proceed if bounds ok
	*/
	npts=(double)(chipSize)*(double)(chipSize);
	/* Skip if any missing data in search window */
	varL=0.0;meanL=0.0;
	for(ii=iiMin; ii <=iiMax; ii++) {
		for(jj=jjMin; jj <=jjMax; jj++) {
			if(lateImage->fimage[il+ii][jl+jj]<(NODATA+1)) { return(NOUSEABLEDATA);}
			varL+=    lateImage->fimage[il+ii][jl+jj] * lateImage->fimage[il+ii][jl+jj]; /* Sum variables for late window, which is fixed in extent */
			meanL+=lateImage->fimage[il+ii][jl+jj];
		}
	}
	varL=varL/npts; meanL=meanL/npts; sigmaL2=varL-meanL*meanL; /* Normalize variables */
	/* 
	   Return if any missing points in early image 
	*/
	for(ii=iiMin+ly1; ii <=iiMax+ly2; ii++) {
		for(jj=jjMin+lx1; jj <=jjMax+lx2; jj++) {
			if(earlyImage->fimage[ie+ii][je+jj] <(NODATA+1)) { return(NOUSEABLEDATA);}
		}
	}
	/* -----------------------------------------------------------------------------------------------------------
	   Main correlation code
	   ----------------------------------------------------------------------------------------------------------- */
	maxCorr=-1.0; jMax=1000; iMax=1000;
	for(i=ly1; i <=ly2; i++) { 
		for(j=lx1; j <=lx2; j++) {
			corr[i][j]=0.0;  Xcorr=0.0; varE=0.0; meanE=0.0; sigmaE2=0.0;
			/* Inner correlation loop */
			for(ii=iiMin; ii <=iiMax; ii++) {
				for(jj=jjMin; jj <=jjMax; jj++) {
					Xcorr += (double)earlyImage->fimage[ie+i+ii][je+j+jj] * 
						(double)lateImage->fimage[il+ii][jl+jj];
					varE += ((double)earlyImage->fimage[ie+i+ii][je+j+jj]) * 
						((double)earlyImage->fimage[ie+i+ii][je+j+jj]);
					meanE += (double)earlyImage->fimage[ie+i+ii][je+j+jj];
				}
			} /* end inner correlation */
			varE=varE/npts; meanE=meanE/npts; sigmaE2=varE-meanE*meanE; /* Normalize variables */
			Xcorr=Xcorr/npts;
			corr[i][j]=(float)( (Xcorr-meanE*meanL)/sqrt(sigmaE2*sigmaL2));
			if( corr[i][j] > maxCorr) {
				maxCorr=corr[i][j]; iMax=i; jMax=j; sigmaMaxE=sigmaE2;  sigmaMaxL=sigmaL2; 
			} /* end if */
		} /* j: end main correlatin loop */
	} /* i: end main correlation loop */
	/*
	  reject rho
	*/
	if(maxCorr < RHOTHRESH) {*offX=NODATA; *offY=NODATA; return(NOCORR);}
	/*-----------------------------------------------------------------------------------------------------------
	  Check if correlation result needs to be extended so that the oversample can center on pean +/- OSSIZE
	  -----------------------------------------------------------------------------------------------------------*/
	if(jMax > lx2 || jMax < lx1 | iMax > ly2 || iMax < ly1) { *offX=NODATA; *offY=NODATA; return(NOCORR); }
	lx1a=lx2+1; lx2a=lx1-1; ly1a=ly2+1; ly2a=ly1-1; /* Will cause to skip loop unless modified below */
	if( (lx2-jMax) <= OSSIZE) lx2a=jMax+OSSIZE;
	if( (ly2-iMax) <= OSSIZE) ly2a=iMax+OSSIZE;
	if( (jMax-lx1) <= OSSIZE) lx1a=jMax-OSSIZE;
	if( (iMax-ly1) <= OSSIZE) ly1a=iMax-OSSIZE;
	/*
	  Extend rows if needed
	*/
	for(i=ly1a; i <=ly2a; i++) { 
		for(j=jMax-OSSIZE; j <=lx2a+OSSIZE; j++) {
			Xcorr=0.0; varE=0.0; meanE=0.0; sigmaE2=0.0;  /* Zero summation variables */
			/* Inner correlation loop */
			for(ii=iiMin; ii <=iiMax; ii++) {
				for(jj=jjMin; jj <=jjMax; jj++) {
					Xcorr+=  (double)earlyImage->fimage[ie+i+ii][je+j+jj] *      (double)lateImage->fimage[il+ii][jl+jj];
					varE+=   ((double)earlyImage->fimage[ie+i+ii][je+j+jj]) * ((double)earlyImage->fimage[ie+i+ii][je+j+jj]);
					meanE+=(double)earlyImage->fimage[ie+i+ii][je+j+jj];
				}
			} /* end inner correlation */
			varE=varE/npts; meanE=meanE/npts; sigmaE2=varE-meanE*meanE; /* Normalize variables */
			Xcorr=Xcorr/npts;
			corr[i][j]=(float)( (Xcorr-meanE*meanL)/sqrt(sigmaE2*sigmaL2));
		} /* j: end main correlation loop */
	} /* i: end main correlation loop */
	/*
	  Extend columns if needed,
	*/
	for(i=iMax-OSSIZE; i <=iMax+OSSIZE; i++) {  /* technically this will excute empty loop if not needed */
		for(j=lx1a; j <=lx2a; j++) {
			Xcorr=0.0; varE=0.0; meanE=0.0; sigmaE2=0.0; /* Zero summation variables */
			/* Inner correlation loop */
			for(ii=iiMin; ii <=iiMax; ii++) {
				for(jj=jjMin; jj <=jjMax; jj++) {
					Xcorr+=  (double)earlyImage->fimage[ie+i+ii][je+j+jj] *      (double)lateImage->fimage[il+ii][jl+jj];
					varE+=   ((double)earlyImage->fimage[ie+i+ii][je+j+jj]) * ((double)earlyImage->fimage[ie+i+ii][je+j+jj]);
					meanE+=(double)earlyImage->fimage[ie+i+ii][je+j+jj];
				}
			} /* end inner correlation */
			varE=varE/npts; meanE=meanE/npts; sigmaE2=varE-meanE*meanE; /* Normalize variables */
			Xcorr=Xcorr/npts;
			corr[i][j]=(float)( (Xcorr-meanE*meanL)/sqrt(sigmaE2*sigmaL2));
		} /* j: end main correlatin loop */
	} /* i: end main correlation loop */
	/*-----------------------------------------------END EXTRA CORRELATION----------------------------------------------*/
	/*
	  Oversample correlation (returns maxCorrHi, and offX, offY)
	*/
	maxCorrHi=overSamplePeakCorr(iMax, jMax, offX, offY);
	/*
	  In theory late-image patch slides over early image, but they way its implemented above, the early image slides under the late patch.
	  As a consequence the signs are flipped. This undoes the sign flip.
	*/
	*offY*=-1;
	*offX*=-1;
	return(maxCorrHi);
}

/********************************** overSamplePeakCorr ********************************************
   Routine to overample correlation peak
   input iMax,jMax; Results are calculated using global array. 
   ******************************************** ********************************************************* */
static float overSamplePeakCorr(int32 iMax, int32 jMax, float *offX, float *offY)
{
	extern float **corr;
	extern fftwf_complex *corrPatchSpace, *corrPatchFreq, *corrPatchOverFreq,*corrPatchOverSpace;
	extern fftwf_plan forwardCorrPlan, inverseCorrPlan;
	int32 i,j,i1,j1,i2,j2,width,widthO, iMax1,jMax1;
	float maxCorrHi;
	float tmpCorr;
	float Crms;
	width=2*OSSIZE;
	widthO=2*OSSIZE*NOVER;
	/*
	  Load array for fft 
	*/
	i1=0; j1=0;
	Crms=0;
	for(i=iMax-OSSIZE; i <(iMax+OSSIZE); i++) {  
		j1=0;
		for(j=jMax-OSSIZE; j < (jMax+OSSIZE); j++) {
			corrPatchSpace[i1*width+j1][0]=corr[i][j]; 
			corrPatchSpace[i1*width+j1][1]=0.0;
			Crms+=corr[i][j]*corr[i][j];
			j1++;
		}
		i1++;
	}
	Crms -= corr[iMax][jMax]*corr[iMax][jMax] ; /* to sample only area around the peak */

	Crms=sqrt((double)Crms/(2*2*OSSIZE*OSSIZE));
	/*
	  Reject for low correlation AND poor peak.
	 */
	if( (sqrt(corr[iMax][jMax]) < RHOTHRESHFINAL) && (sqrt(corr[iMax][jMax])/Crms < 9.0)) {
		*offX=NODATA; *offY=NODATA;
		return(NOCORR);
	}

	/* DEBUG */
	/*
	 *offY=Crms; *offX=sqrt(corr[iMax][jMax])/Crms;
	 return(corr[iMax][jMax]);*/
	/* end DEBUG */
	/*
	  Forward FFT
	*/
	fftwf_execute ( forwardCorrPlan );
	/*
	  Fill matrix for oversample
	*/
	for(i=0; i <OSSIZE; i++ ) {
		i1=widthO - OSSIZE +i;
		i2=OSSIZE + i;
		j1=widthO- OSSIZE;
		j2=OSSIZE;
		for(j=0; j < OSSIZE; j++) {
			corrPatchOverFreq[i*widthO + j][0] = corrPatchFreq[i*width + j][0]; 
			corrPatchOverFreq[i*widthO + j][1] = corrPatchFreq[i*width + j][1];
			corrPatchOverFreq[i1*widthO+j][0] = corrPatchFreq[i2*width + j][0];
			corrPatchOverFreq[i1*widthO+j][1] = corrPatchFreq[i2*width + j][1];
			corrPatchOverFreq[i1*widthO+j1][0] = corrPatchFreq[i2*width + j2][0];
			corrPatchOverFreq[i1*widthO+j1][1] = corrPatchFreq[i2*width + j2][1];
			corrPatchOverFreq[i*widthO + j1][0] = corrPatchFreq[i*width + j2][0];
			corrPatchOverFreq[i*widthO + j1][1] = corrPatchFreq[i*width + j2][1];
			j1++; j2++;
		} /* End for j */
	} /* End for i */
	fftwf_execute ( inverseCorrPlan );
	/* 
	   Find peak of oversampled  - search two original pixels (NOVER*2) around peak
	*/
	iMax1=NODATA;	jMax1=NODATA;	maxCorrHi=-(float) (OSSIZE*2*OSSIZE*2);
	for(i=(OSSIZE*NOVER)-NOVER*2; i <  (OSSIZE*NOVER) + NOVER*2; i++)  {
		for(j=(OSSIZE*NOVER)-NOVER*2; j <  (OSSIZE*NOVER) + NOVER*2; j++) {
			tmpCorr=(float)corrPatchOverSpace[i*widthO+j][0] ;
			if(  tmpCorr > maxCorrHi) {
				{maxCorrHi= tmpCorr; iMax1=i; jMax1=j;}   
			}
		}
	}
	maxCorrHi /=(float)( (OSSIZE*2)*(OSSIZE*2));
	*offX=jMax+(float)jMax1/(float)(NOVER) - OSSIZE;
	*offY=iMax+(float)iMax1/(float)(NOVER) - OSSIZE;
	
	return(maxCorrHi);
}

/***************************************GetLatePosition ************************************************
  Given position in early image, take first cut at location it in second image (largley dead reckoning using images coordinates 
*******************************************************************************************************/
static int32 getLatePosition(LSimage *earlyImage,LSimage *lateImage,matchParams *matchP,uint32 ie,uint32 je,uint32 *il,uint32 *jl,int32 *chipSize,int32 *lx1, int32 *lx2,int32 *ly1,int32 *ly2)	 
{
	double xe,ye;
	double xj,yi;
	/* 
	   Absolute position in early image
	*/
	xe=earlyImage->x0 + je*earlyImage->dx;
	ye=earlyImage->y0 + ie*earlyImage->dy;
	/* 
	   fractional image coordinate in image2
	*/
	xj=(xe-lateImage->x0)/lateImage->dx;
	yi=(ye-lateImage->y0)/lateImage->dy;
	/*
	  Chip size
	*/
	*chipSize=matchP->chipSize; /* propagate default */
	/* ADD FANCIER CHIP SIZE LATER AS NEEDED */
	*lx1=-matchP->searchX;
	*ly1=-matchP->searchY;
	*lx2=matchP->searchX;
	*ly2=matchP->searchY;
	/*
	  Compute return values	
	*/
	if(xj < 0 || yi < 0) { *il=-1; *jl=-1;return(-1);}

	*il=(int32)(yi + 0.5);
	*jl=(int32)(xj + 0.5);
	return(1);
  
}
/*
*************************** loadTiff - read a tiff file in ***************************************
*/
void loadTiff(char *file,LSimage *Image,TIFF *tif,GTIF *gtif)
{
	uint16 *buf, *tileBuf, **tile, **Im,*stripBuf;     /* Buffers to read in image */
	float *fbuf;
	uint32 i,j,i1,j1;  		           /* LCV */		
	size_t bufSize;		           /* Used to do buffer size math */
	uint32 tileWidth, tileLength;  /* Tile size info */
	uint32 x, y;				  /* Used to index into tiles */
	int32 ik,jk,nk;			/* loop indices, for filter convolution, nk=filter halfwidth width=2*nk+1 */
	uint32 idum;			 	 /* Dummy variable */
	double fdum;			 	 /* Dummy variable */
	float **fk;				/* High-pass filter kernerl */
	float tmpF;				/* Debug */
	double xc,yc,xc1,yc1,dx,dy;   /* Variables used to compute corners location and pixel size */
	uint32 *bc, *sc;
	int32 rowsperstrip,stripoffsets,stripbytecounts,strip,idum1;
	/*
	  Filter kernel - uses a straightforward edge filter from Ahn and Howat - could be made more customizable
	*/
	nk=1;
	fk=lsmatrix(-nk,nk,-nk,nk );
	/* Edge filter as default */
	for(ik=-nk; ik <=nk; ik++)	for(jk=-nk; jk <=nk; jk++) fk[ik][jk]=-1.0;
	fk[0][0]=8.0;
	fprintf(stderr,"--------------- Filter Kernel ----------\n");
	for(ik=-nk; ik <=nk; ik++)	{ for(jk=-nk; jk <=nk; jk++) { fprintf(stderr,"%8.3f ",fk[ik][jk]); } fprintf(stderr,"\n");}
	fprintf(stderr,"--------------------------------- ----------\n");
	/* 
	   Tiff data
	*/
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &idum); Image->nx=idum;
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &idum); Image->ny=idum;
	fprintf(stderr,"Image %s width/length %i %i\n",file,Image->nx, Image->ny);
	/*
	  Get lower left corner
	*/
	xc=0.5;yc=0.0+Image->ny -0.5;
	GTIFImageToPCS( gtif, &xc, &yc );
	Image->x0=xc; Image->y0=yc;
	/*
	  Maybe there is an easier way, but brute force compute pixel size
	*/
	xc=0.0; yc=0.0; xc1=1; yc1=1;     /* Coordinates from origin and one pixel from origin */
	GTIFImageToPCS( gtif, &xc, &yc);
	GTIFImageToPCS( gtif, &xc1, &yc1);
	dx=xc1-xc; if(dx < 0) dx *=-1.0;  /* Compute dx */
	dy=yc1-yc; if(dy < 0) dy *=-1.0;   /* Compute dx */
	Image->dx=dx; Image->dy=dy;
	fprintf(stderr,"x0/y0/dx/dy %lf %lf %lf %lf\n",Image->x0,Image->y0,Image->dx,Image->dy);
	idum1=0;
	GTIFKeyGet(gtif,ProjectedCSTypeGeoKey,&idum1,0,1);
	fprintf(stderr,"EPSG %i\n",idum1); Image->proj=idum1;
        if(Image->proj != NSIDCNORTH && Image->proj != NSIDCSOUTH) Image->proj=-1;
	/* 
	   Malloc image buffer
	*/
	bufSize=Image->ny * Image->nx * sizeof(uint16);
	buf=malloc(bufSize);
	bufSize=Image->ny * sizeof(uint16 *);
	Image->image=(uint16 **)malloc(bufSize);
	for(i=0; i < Image->ny; i++) Image->image[i] = &(buf[i*Image->nx]);
	/* 
	   Malloc floating point  image buffer
	*/
	bufSize=Image->ny * Image->nx * sizeof(float);
	fbuf=(float *)malloc(bufSize);
	bufSize=Image->ny * sizeof(float *);
	Image->fimage=(float **)malloc(bufSize);
	for(i=0; i < Image->ny; i++) Image->fimage[i] = &(fbuf[i*Image->nx]);
	/*
	  Strip info
	*/
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
	TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &bc);

	fprintf(stderr,"rows per strip, byte count %i %i %li %i\n",rowsperstrip,bc[0],TIFFStripSize(tif) ,TIFFIsTiled(tif));
	if(TIFFIsTiled(tif) ==0) {  /* Striped */
		fprintf(stderr,"Strip input \n");
		stripBuf = _TIFFmalloc(TIFFStripSize(tif));
		for (strip = 0; strip < TIFFNumberOfStrips(tif); strip++) {
			TIFFReadEncodedStrip(tif, strip, stripBuf, (tsize_t) -1);
			for(j=0; j < Image->nx; j++) { 
				Image->image[Image->ny-1-strip][j]=stripBuf[j]; 
			}
		} fprintf(stderr,"strips %i\n",strip);
	} else { /* Tiled */
		fprintf(stderr,"Tile input \n");

		TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength); 		/*		  Get tile size		*/
		TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
		fprintf(stderr,"Tile Size %li %i %i\n", TIFFTileSize(tif),tileLength,tileWidth);

		/*
		  Malloc tile buffer
		*/ 
		tileBuf = _TIFFmalloc(TIFFTileSize(tif));
		bufSize=tileLength * sizeof(uint16 *);
		tile=(uint16 **)malloc(bufSize);
		for(i=0; i < tileLength; i++) tile[i]=&(tileBuf[i*tileWidth]);
		/*
		  read tiles
		*/
		for (y = 0; y < Image->ny; y += tileLength)
			for (x = 0; x < Image->nx; x += tileWidth) {
				TIFFReadTile(tif, tileBuf, x, y, 0,0);
				i1=0; 
				for(i=y; i < min((y+tileLength),Image->ny); i++ ) {
					j1=0;
					for(j=x; j < min((x+tileWidth),Image->nx); j++) { 
						Image->image[Image->ny-1-i][j]=tile[i1][j1]; 
						j1++;
					}
					i1++;
				}
			}
		free(tile); free(tileBuf);
	}
	/*
	  High pass filter and load into floating point array
	*/
	for(i=nk; i < Image->ny-nk; i++) {
		for(j=nk; j < Image->nx-nk; j++)  {
			/*	Image->fimage[i][j]=0.0;*/
			Image->fimage[i][j]=fk[0][0]* (float)Image->image[i][j]; /* Center element */
			for(ik=1; ik <= nk; ik++)  {
				Image->fimage[i][j] +=(float)(((float)fk[0][ik]  * (float)Image->image[i][j+ik]) + (float)((float)fk[0][-ik]  * (float)Image->image[i][j-ik])+
							      (float)((float)fk[ik][0]  * (float)Image->image[i+ik][j]) + (float)((float)fk[-ik][0]  * (float)Image->image[i-ik][j])); /* cross "+" elements */
				for(jk=1; jk <= nk; jk++)  { /* non-zero elements */
					Image->fimage[i][j] +=(float)(((float)fk[ik][jk]  * (float)Image->image[i+ik][j+jk])+ ((float)fk[-ik][-jk]* (float)Image->image[i-ik][j-jk])  +
								      ((float)fk[-ik][jk] * (float)Image->image[i-ik][j+jk])  + ((float)fk[ik][-jk]* (float)Image->image[i+ik][j-jk]) ); 
				}
			}
			
			/*	Image->fimage[i][j]= 8.0*(float)Image->image[i][j]
				tmpF= 8.0*(float)Image->image[i][j]
				-(float) (  Image->image[i-1][j-1]+Image->image[i-1][j] + Image->image[i-1][j+1]  
				+   Image->image[i][j-1] +					      Image->image[i][j+1]  +  
				Image->image[i+1][j-1]+Image->image[i+1][j]+Image->image[i+1][j+1]  ); 
			
				if(i > 3000 && i < 3100 && j==6000) fprintf(stderr,"%f %f \n",tmpF,	Image->fimage[i][j] );
			*/
			if(Image->image[i][j] ==0) {Image->fimage[i][j] =(float) NODATA;} /* Use our NODATA */
		}
	}
	/* Edges */	     
	for(i=0; i < Image->ny; i++)  for(j=0; j < nk; j++)  { 
		Image->fimage[i][j]=(float) NODATA; Image->fimage[i][Image->nx-j-1]=(float) NODATA;
	}
	for(j=0; j < Image->nx; j++) {}
		for(i=0; i < nk; i++) { 
			Image->fimage[i][j]= (float) NODATA;Image->fimage[Image->ny-i-1][j]=(float) NODATA;
		}
	}
	free( Image->image[0]);
	/*
	  Parse date and path/row from file name (may change means to get this information later 
	*/
	parseLSDate(file,Image);
}
/**************************** parseLSDate-parse LS filename to get date, path, row ***************************************

*************************************** *************************************** ****************************************/
static void parseLSDate(char *file, LSimage *Image)
{
	int32 i,j, slash,idum,leadZero,path,row,year,doy,jd;
	char tmp[1500],*tmp1;
	i=0; 
	slash=0;
        while(file[i] != '\0') {
		if(file[i]=='/') slash=i;
		i++;
		if(i > 1500) error("Error parsing landsat file %s\n",file);
	}      
	tmp1=&(file[slash+1]);
        fprintf(stderr,"%s %s \n",file,tmp1);
	fprintf(stderr,"%i \n",slash);
	Image->sat=(char *)malloc(4*sizeof(char)); for(i=0; i<3; i++) Image->sat[i]=tmp1[i]; Image->sat[3]='\0';
        if( tmp1[2] == '0') tmp1 +=4; else tmp1+=3;
        leadZero=-1;
	for(i=0; i < 3; i++) { tmp[i]=tmp1[i]; if(leadZero < 0 && tmp[i]=='0') tmp[i]=' '; else leadZero=1; } tmp[3]='\0';
	sscanf(&tmp[0],"%i",&path); 
	tmp1+=3;
        leadZero=-1;
	for(i=0; i < 3; i++) { tmp[i]=tmp1[i]; if(leadZero < 0 && tmp[i]=='0') tmp[i]=' '; else leadZero=1; } tmp[3]='\0';
	sscanf(&tmp[0],"%i",&row); 
	tmp1+=3;
        leadZero=-1;
	for(i=0; i < 4; i++) { tmp[i]=tmp1[i]; if(leadZero < 0 && tmp[i]=='0') tmp[i]=' '; else leadZero=1; } tmp[4]='\0';
	sscanf(tmp,"%i",&year); 

	tmp1+=4;
        leadZero=-1;
	for(i=0; i < 3; i++) { tmp[i]=tmp1[i]; if(leadZero < 0 && tmp[i]=='0') tmp[i]=' '; else leadZero=1; } tmp[3]='\0';
	sscanf(tmp,"%i",&doy); 

        Image->path=path;
        Image->row=path;
	Image->year=year;
	Image->doy=doy;
	jd=julday(1,1,year)+doy-1;
	Image->jd=(double)jd;
	fprintf(stderr,"path,row,year,doy %7i %7i %7i %7i %10i\n",path,row,year,doy,jd);      
}
	    
/****************************mallocMatch******************************************************************
-> Malloc space for fft
-> setup fft plans
-> malloc output buffers
************************************************************************************************************/
static void mallocMatch(matchResult *matches) 
{
	float *buf;				   /* Pointer to malloced memory for correlation buffer */
	uint8 *bufu;				   /* Pointer to malloced memory for correlation buffer */
	size_t bufSize;		           /* Used to do buffer size math */
	int32_t width,widthO;
	uint32 i,j;				 /* LCV */
	extern float **corr;
	extern fftwf_complex *corrPatchSpace, *corrPatchFreq, *corrPatchOverFreq,*corrPatchOverSpace;
	extern fftwf_plan forwardCorrPlan, inverseCorrPlan;
	/*
	  setup fft
	*/
	widthO=NOVER*(OSSIZE*2);
	width=OSSIZE*2;
	corrPatchSpace= fftwf_malloc (sizeof(fftwf_complex) * width*width);
	corrPatchFreq= fftwf_malloc (sizeof(fftwf_complex) * width*width);
	corrPatchOverSpace= fftwf_malloc (sizeof(fftwf_complex) * widthO*widthO);
	corrPatchOverFreq  =  fftwf_malloc (sizeof(fftwf_complex) *  widthO*widthO);
	fprintf(stderr,"Start plan\n");
	forwardCorrPlan=fftwf_plan_dft_2d(width, width, corrPatchSpace, corrPatchFreq, FFTW_FORWARD, FFTW_ESTIMATE);
	inverseCorrPlan= fftwf_plan_dft_2d(widthO, widthO, corrPatchOverFreq, corrPatchOverSpace, FFTW_BACKWARD, FFTW_ESTIMATE);
	for(i=0; i < widthO; i++) {
		for(j=0; j < widthO; j++) {
			corrPatchOverFreq[i*widthO+j][0]=0.0;corrPatchOverFreq[i*widthO+j][1]=0.0;
		}
	}
	fprintf(stderr,"End plan\n");
	/*
	  Malloc buffer for correlation - Note setting up to allow negative indices
	*/
	corr = lsmatrix(-MAXW, MAXW, -MAXW, MAXW);
	/*
	  Malloc Output Buffers
	*/
	bufSize= sizeof(float *) * matches->ny;
	fprintf(stderr,"bufSize %li\n", bufSize);
	matches->X = (float **)malloc(bufSize);
	matches->Y = (float **)malloc(bufSize);
	matches->Rho = (float **)malloc(bufSize);

	bufSize=sizeof(uint8 *)*matches->ny;
	matches->type = (uint8 **)malloc(bufSize);

	bufSize=sizeof(float)*matches->ny*matches->nx;
	fprintf(stderr,"bufSize %li\n",bufSize);

	buf = (float *)malloc(bufSize);
	for( i=0; i < matches->ny; i++) matches->Rho[i]=&(buf[i*matches->nx]);
    
	buf = (float *)malloc(bufSize);
	for( i=0; i < matches->ny; i++) matches->X[i]=&(buf[i*matches->nx]);

	buf = (float *)malloc(bufSize);
	for(i=0; i < matches->ny; i++) matches->Y[i]=&(buf[i*matches->nx]);

	bufSize=sizeof(uint8)*matches->ny*matches->nx;
	fprintf(stderr,"bufSize %li\n",bufSize);
	bufu = (uint8 *)malloc(bufSize);
	for(i=0; i < matches->ny; i++) matches->type[i]=&(bufu[i*matches->nx]);
	fprintf(stderr,"Done with Malloc\n");
}


/**************************** lsmatrix ******************************************************
Malloc a floating point matrix - in this case mallocs the corr matrix to allow negative indices
*/
float **lsmatrix(int32 nrl, int32 nrh, int32 ncl, int32 nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int32 i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow)*sizeof(float*)));
	if (!m) error("allocation failure 1 in matrix()");
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol)*sizeof(float)));

	if (!m[nrl]) error("allocation failure 2 in matrix()");
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
static long julday(int32_t mm, int32_t id, int32_t iyyy)
{
        long jul;
        int32_t ja,jy=iyyy,jm;

        if (jy < 0) ++jy;

        if (mm > 2) {
                jm=mm+1;
        } else {
                --jy;
                jm=mm+13;
        }
	/*(15+31*(10+12*1582))*/
        jul = (long) (floor(365.25*jy)+floor(30.6001*jm)+id+1720995);

        if ( id+31L*(mm+12*iyyy) >=  (15+31*(10+12*1582))) {
                ja=(int)(0.01 * jy);
                jul += 2-ja+ (int)(0.25*ja);
        }
        return jul;
}
/*
*************************** writeTiffAsRaw - DEBUG routine ***************************************
*/
static void writeTiffAsRaw(LSimage *Image,char *file)
{
	uint32 i;
	FILE *fp;
	fp=fopen(file,"w");
	/*      fwriteBS(Image->image[0],sizeof(uint16),(size_t)(Image->nx*Image->ny),fp,INT16FLAG); */
	fwriteBS(Image->fimage[0],sizeof(float),(size_t)(Image->nx*Image->ny),fp,FLOAT32FLAG); 
}