#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "landsatSource64/clib/standard.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"

#define MAXTIEPOINTS 500000

/*
  Read tiepoint file for tiepoints.
*/
void readLSTiePoints(char *tieFile, lsTiepoints *tiePoints)
{
	int linelength;   /* Input line length */
	int notdone;      /* Loop flag */
	double z=0;
	double lat,lon;
	double vx=0.0, vy=0.0,vz=0.0;
	char lineBuffer[LINEMAX+1];
	char *line;       /* Input line buffer */
	int lineCount=0;
	FILE *fp;
	/*
	  Open Input file
	*/
	fp=openInputFile(tieFile);
	/*
	  Malloc space
	*/
	line = lineBuffer;  /* Allocate line buffer */
	notdone = TRUE;
	tiePoints->lat = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->lon = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->x = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->y = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->z = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->vx = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->vy = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->offXT = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->offYT = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->xyScale = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	/*	tiePoints->vz = (double *)malloc(sizeof(double)*MAXTIEPOINTS);*/
	/*
	  Loop to read tiepoints
	*/
	tiePoints->npts=0;
	while( notdone == TRUE ) {                 /* Loop to read lines */
		linelength = fgetline(fp,line,LINEMAX); /* Read line */
		lineCount++;
		if( strchr(line,ENDDATA) != NULL ) 
			notdone = FALSE;                    /* End of data, set exit flag */
		else if( strchr(line,COMMENT) == NULL ) {   /* If not comment, parse */
			
			if(sscanf(line,"%lf%lf%lf%lf%lf%lf",
				  &lat,&lon,&z,&vx,&vy,&vz) != 6) {/*motion ties*/
				fprintf(stderr,"%f %f %f %f %f %f",lat,lon,z,vx,vy,vz);
				error("%s  %i",
				      "readTiePoints: -- Invalid # of parameters at line:",
				      lineCount);
			}	if(sscanf(line,"%lf%lf%lf%lf%lf",	  &lat,&lon,&z,&vx,&vy) != 5) {/*motion ties*/
				fprintf(stderr,"%lf %lf %lf %lf %lf %lf",lat,lon,z,vx,vy,vz);
				error("%s  %i",
				      "readTiePoints: -- Invalid # of parameters at line:",
				      lineCount);
			}
			/* 
			   Assign tiepoints and update counter
			*/
			tiePoints->lat[tiePoints->npts] = lat;
			tiePoints->lon[tiePoints->npts] = lon;
			tiePoints->z[tiePoints->npts] = z;    
			tiePoints->vx[tiePoints->npts] = vx;     
			tiePoints->vy[tiePoints->npts] = vy;     
			/*				fprintf(stderr,"+++ %f %f %f %f %f\n",lat,lon,z,vx,vy);*/
			/*			tiePoints->vz[tiePoints->npts] = vz;       */
			(tiePoints->npts)++;
			if(tiePoints->npts >= MAXTIEPOINTS) /* Too many tiepts ? */
				error("readTiePoints -- MAXTIEPOINTS=%i exceeded",
				      MAXTIEPOINTS);
		} /* End else */
	} /* End while */
}
