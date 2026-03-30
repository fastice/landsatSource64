#include "stdio.h"
#include "stdlib.h"
#include "geotiff/xtiffio.h"   /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "string.h"
#include "landsatSource64/clib/standard.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "unistd.h"

static void readArgs(int argc, char *argv[], lsFit *fitDat, char **maskFile);
static void dataWidth(matchResult *matches);
static void LSFitusage();
int32_t nMatch, nAttempt, nTotal;
int main(int argc, char *argv[])
{
	char *tieFile, *matchFile, *fitFile;
	char *offsetDatFile;
	matchParams matchP;
	matchResult matches;
	lsTiepoints tiePoints;
	double xImageWidth, yImageWidth;
	lsFit fitDat;
	char *maskFile;
	/*
		Read Args
	*/
	readArgs(argc, argv, &fitDat, &maskFile);

	fprintf(stderr, "Tiefile File                            %s\n", fitDat.tieFile);
	fprintf(stderr, "Match File                            %s\n", fitDat.matchFile);
	fprintf(stderr, "fitFile File                            %s\n", fitDat.fitFile);
	/*
	   Input offsets
	*/
	fprintf(stderr, "Reading Offsets \n");
	readLSOffsets(&fitDat, &matches, TRUE, maskFile);
	// Compute Nominal Width
	dataWidth(&matches);
	/*
		Read Tie points
	*/
	fprintf(stderr, "Reading Tiepoints \n");
	readLSTiePoints(fitDat.tieFile, &tiePoints);
	/*
		Interp offsets
	*/
	fprintf(stderr, "Interpolating offset data \n");
	interpLSTies(&fitDat, &matches, &tiePoints);
	/*
		 Do fit for LS parameters
	*/
	computeLSFit(&fitDat, &matches, &tiePoints);
}

// Compute max width of the data in the x and y directions.
static void dataWidth(matchResult *matches)	
{
	int minj=2000000, maxj=0, mini=2000000, maxi=0;
	for(int i=0; i < matches->ny; i++)
	{
		for(int j=0; j < matches->nx; j++)
		{
			if(matches->X[i][j] > (NODATA + 1))
			{
				mini = min(i, mini);
				maxi = max(i, maxi);
				minj = min(j, minj);
				maxj = max(j, maxj);
			}
		}
	}
	fprintf(stderr, "min/max i, j %i %i, %i, %i\n", mini, maxi, minj, maxj);
	matches->yDataWidth = (maxi - mini + 1) * matches->dx * (double)matches->stepX * MTOKM;
	matches->xDataWidth = (maxj - minj + 1) * matches->dy * (double)matches->stepY * MTOKM;
	fprintf(stderr, "data X/Y width %lf %lf\n", matches->xDataWidth, matches->yDataWidth);
}


static void readArgs(int argc, char *argv[], lsFit *fitDat, char **maskFile)
{
	int32_t n, i, idum;
	char *argString;
	char *dum;

	*maskFile = NULL;
	if (argc < 3 || argc > (4 + 13))
		LSFitusage(); /* Check number of args */
	fitDat->tieFile = argv[argc - 3];
	fitDat->matchFile = argv[argc - 2];
	fitDat->fitFile = argv[argc - 1];
	fitDat->fitFile = malloc(strlen(argv[argc - 1]) + 1);
	dum = strcpy(fitDat->fitFile, argv[argc - 1]);
	fitDat->fitType = PLANEFIT;
	n = argc - 4;
	fprintf(stderr, "%i %i\n", argc, n);
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "mask") != NULL)
		{
			*maskFile = argv[i + 1];
			i++;
			fprintf(stderr, "Using mask : %s\n ", *maskFile);
			if (access(*maskFile, F_OK) == -1)
			{
				error("Mask file not found %s\n", *maskFile);
			}
		}
		else
			LSFitusage();
	}
	return;
}

static void LSFitusage()
{
	error("\n\n%s\n %s \n  %s \n   %s \n   %s \n  ",
		  "lsfit -mask maskFile tieFile matchFile fitFile \n",
		  "where \n",
		  "tieFile =  tiepoints file ",
		  "matchFile = root for matches (.offX,.offY) ",
		  "fitFile  = file with fit parameters  ");
}
