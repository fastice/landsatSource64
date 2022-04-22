#include "stdio.h"
#include "stdlib.h"
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "lstrack.h"
#include "string.h"
#include "landsatSource64/clib/standard.h"
static void readArgs(int argc,char *argv[],char **earlyFile,char **lateFile,matchParams *matchP);
static void LSusage();
static void writeLSOffsets(matchResult *matches,matchParams *matchP);
matchResult runLStrack(char *earlyFile, char *lateFile, matchParams *matchP );
static void readLSMaskFile(matchParams *matchP);
static void readVelMosaic(matchParams *matchP);
int32 nMatch, nAttempt, nTotal;

void main(int argc, char *argv[])
{   
	char *parFile;
 
	int noComplex;
	int floatFlag;
	char *earlyFile, *lateFile, *outfile;
	matchParams matchP;
	matchResult matches;						/* match result */

	TIFF *tif=(TIFF*)0;  /* TIFF-level descriptor */
	GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */
	/*
	  Read input files
	*/
	readArgs(argc,argv,&earlyFile,&lateFile,&matchP);
	fprintf(stderr,"------------ MATCH PARAMS --------- \n");
	fprintf(stderr,"slowFlag = %9i\n",matchP.slowFlag);
	fprintf(stderr,"Chip Size                                            %5i\n",matchP.chipSize);
	fprintf(stderr,"Search Radius (searchX,searchY) %5i  %5i\n",matchP.searchX,matchP.searchY);
	fprintf(stderr,"Step X and Step Y                            %5i  %5i\n",matchP.stepX,matchP.stepY);
        if(matchP.maskFile != NULL) fprintf(stderr,"maskFile = %s\n",matchP.maskFile); else fprintf(stderr,"maskFile = NULL\n");
        if(matchP.velFile != NULL) fprintf(stderr,"velFile = %s\n",matchP.velFile); else fprintf(stderr,"velFile = NULL\n");
	fprintf(stderr,"Early File                            %s\n", earlyFile);
	fprintf(stderr,"Late File                            %s\n", lateFile);
	fprintf(stderr,"Output File                            %s\n",matchP.outputFile);
	fprintf(stderr,"------------ END MATCH PARAMS --------- \n");
	/*
	  Read mask file
	*/
	if(matchP.maskFile !=NULL) readLSMaskFile(&matchP);

	/*
	  Read mask file
	*/
	if(matchP.velFile !=NULL) readVelMosaic(&matchP);
	/*
	  Call matcher
	*/
	matches=runLStrack(earlyFile, lateFile, &matchP);
	matches.fileEarly=earlyFile;
	matches.fileLate=lateFile;
	fprintf(stderr,"%s \n",matchP.outputFile);
	writeLSOffsets(&matches,&matchP);
	fprintf(stderr,"--- %i %i %i \n",nMatch,nAttempt,nTotal);
} 

/**************************** writeLSOffsets ******************************************************
  Output match results
****************************************************************************************************/
static void writeLSOffsets(matchResult *matches,matchParams *matchP)
{
	uint32 i;
	FILE *fp;
	float successRate;
	char *file1;
	size_t  sl;
	extern int32 nMatch, nAttempt, nTotal;
	sl=1500; 
	file1=(char *)malloc(sl); 

	fprintf(stderr,"File root %s %i %i\n",matchP->outputFile, matches->nx,matches->ny);
	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".rho"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(matches->Rho[0],sizeof(float),(size_t)(matches->nx*matches->ny),fp,FLOAT32FLAG); fclose(fp);

	fprintf(stderr,"C\n");
	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".dx"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(matches->X[0],sizeof(float),(size_t)(matches->nx*matches->ny),fp,FLOAT32FLAG); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".dy"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(matches->Y[0],sizeof(float),(size_t)(matches->nx*matches->ny),fp,FLOAT32FLAG); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".mtype"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwrite(matches->type[0],sizeof(uint8),(size_t)(matches->nx*matches->ny),fp); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".dat"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w"); 
	fprintf(fp,"fileEarly = %s\n",matches->fileEarly);
	fprintf(fp,"fileLate = %s\n",matches->fileLate);
	fprintf(fp,"x0 = %11.2lf\n",matches->x0);
	fprintf(fp,"y0 = %11.2lf\n",matches->y0);
	fprintf(fp,"dx = %10.5lf\n",matches->dx);
	fprintf(fp,"dy = %10.5lf\n",matches->dy);
	fprintf(fp,"stepY = %u\n",matches->stepX);
	fprintf(fp,"stepX = %u\n",matches->stepY);
	fprintf(fp,"nx = %9i\n",matches->nx);
	fprintf(fp,"ny = %9i\n",matches->ny);
	fprintf(fp,"slowFlag = %9i\n",matchP->slowFlag);
	fprintf(fp,"EPSG = %9i\n",matchP->proj);
	fprintf(fp,"earlyImageJD = %10lf\n",matches->jdEarly);
	fprintf(fp,"lateImageJD = %10lf\n",matches->jdLate);
	fprintf(fp,"IntervalBetweenImages = %5i\n",(int)(matches->jdLate-matches->jdEarly));
	if(nAttempt > 0) successRate = ( (float)nMatch/(float)nAttempt)*100.; else successRate=0.0;
	fprintf(fp,"Success_rate_for_attempted_matches(%%) =  %7.2f \n", successRate);
	fprintf(fp,"& \n");
	fclose(fp);
}

/*
  Read read velocity mosaic used for initial guess
*/
static void readVelMosaic(matchParams *matchP)
{
	FILE *fp;
	char *geodatFile, *vxFile,*vyFile;
	char line[1500];
	int lineCount=0, eod;
	float *tmp,*tmp1;
	float dum1,dum2;
	int i;
	/*
	  geodat file name
	*/
	geodatFile = (char *)malloc(strlen(matchP->velFile)+11); 	geodatFile[0]='\0';
	geodatFile=strcpy(geodatFile,matchP->velFile);
	geodatFile=strcat(geodatFile,".vx.geodat");
	fprintf(stderr,"vel geodat file %s\n",geodatFile);
	/* 
	   velocity file names
	*/
	vxFile = (char *)malloc(strlen(matchP->velFile)+4); 	vxFile[0]='\0';
	vxFile=strcpy(vxFile,matchP->velFile);	vxFile=strcat(vxFile,".vx");
	vyFile = (char *)malloc(strlen(matchP->velFile)+4); 	vyFile[0]='\0';
	vyFile=strcpy(vyFile,matchP->velFile);	vyFile=strcat(vyFile,".vy");
	fprintf(stderr,"vx,vy file %s %s\n",vxFile,vyFile);
	/*
	  Open geodat file
	*/
	fp = openInputFile(geodatFile);
	if(fp == NULL)
		error("*** readLSVelFile: Error opening %s ***\n",geodatFile);
	/*
	  Read parameters
	*/
	lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # 2 line */
	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%f %f\n",&dum1,&dum2); /* read as float in case fp value */
	matchP->velocity.nx = (int)dum1;
	matchP->velocity.ny= (int)dum2;

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(matchP->velocity.dx),&(matchP->velocity.dy));

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(matchP->velocity.x0),&(matchP->velocity.y0));
	matchP->velocity.x0 *= KMTOMS;     matchP->velocity.y0 *= KMTOMS;
	fclose(fp);

	fprintf(stderr,"%i %i \n %f %f \n %f %f \n",  matchP->velocity.nx,matchP->velocity.ny,
		matchP->velocity.dx,matchP->velocity.dy,  matchP->velocity.x0,matchP->velocity.y0);
	/*
	  Malloc array
	*/
	matchP->velocity.vx = (float **)malloc(matchP->velocity.ny * sizeof(float *));
	tmp = (float *)malloc(matchP->velocity.nx *matchP->velocity.ny * sizeof(float));
	/*
	  Open vx file
	*/

	if(vxFile == NULL)  error("*** readVelMosaic: Error opening %s ***\n",vxFile);
	fp = openInputFile(vxFile);
	if(fp == NULL) error("*** readVelMosaic: Error opening %s ***\n",vxFile);
	for(i=0; i < matchP->velocity.ny; i++) {
		tmp1 =&(tmp[i*matchP->velocity.nx]);
		freadBS(tmp1,sizeof(float),matchP->velocity.nx,fp,FLOAT32FLAG);
		matchP->velocity.vx[i]=tmp1;
	}
	fclose(fp);
	/*
	  Malloc array
	*/
	matchP->velocity.vy = (float **)malloc(matchP->velocity.ny * sizeof(float *));
	tmp = (float *)malloc(matchP->velocity.nx *matchP->velocity.ny * sizeof(float));
	/*
	  Open vx file
	*/
	if(vyFile == NULL)  error("*** readShelf: Error opening %s ***\n",vyFile);
	fp = openInputFile(vyFile);
	if(fp == NULL) error("*** readShelf: Error opening %s ***\n",vyFile);
	for(i=0; i < matchP->velocity.ny; i++) {
		tmp1 =&(tmp[i*matchP->velocity.nx]);
		freadBS(tmp1,sizeof(float),matchP->velocity.nx,fp,FLOAT32FLAG);
		matchP->velocity.vy[i]=tmp1;
	}
	fclose(fp);
}

/*
  Read landsat mask
*/
static void readLSMaskFile(matchParams *matchP)
{
	FILE *fp;
	char *geodatFile;
	char line[1500];
	int lineCount=0, eod;
	unsigned char *tmp,*tmp1;
	float dum1,dum2;
	int i;
	/*
	  geodat file name
	*/
	geodatFile = (char *)malloc(strlen(matchP->maskFile)+8);
	geodatFile[0]='\0';
	geodatFile=strcpy(geodatFile,matchP->maskFile);
	geodatFile=strcat(geodatFile,".geodat");
	fprintf(stderr,"mask geodat file %s\n",geodatFile);
	/*
	  Open geodat file
	*/
	fp = openInputFile(geodatFile);
	if(fp == NULL)
		error("*** readLSMaskFile: Error opening %s ***\n",geodatFile);
	/*
	  Read parameters
	*/
	lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # 2 line */
	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%f %f\n",&dum1,&dum2); /* read as float in case fp value */
	matchP->mask.nx = (int)dum1;
	matchP->mask.ny= (int)dum2;

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(matchP->mask.dx),&(matchP->mask.dy));

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(matchP->mask.x0),&(matchP->mask.y0));
	matchP->mask.x0 *= KMTOMS;     matchP->mask.y0 *= KMTOMS;
	fclose(fp);

	fprintf(stderr,"%i %i \n %f %f \n %f %f \n",  matchP->mask.nx,matchP->mask.ny,
		matchP->mask.dx,matchP->mask.dy,  matchP->mask.x0,matchP->mask.y0);
	/*
	  Malloc array
	*/
	matchP->mask.m = (unsigned char **)malloc(matchP->mask.ny * sizeof(unsigned char *));
	tmp = (unsigned char *)malloc(matchP->mask.nx *matchP->mask.ny * sizeof(unsigned char));
	/*
	  Open shelfFile file
	*/
	if(matchP->maskFile == NULL)  error("*** readShelf: Error opening %s ***\n",matchP->maskFile);
	fp = openInputFile(matchP->maskFile);
	if(fp == NULL)
		error("*** readShelf: Error opening %s ***\n",matchP->maskFile);
	for(i=0; i < matchP->mask.ny; i++) {
		tmp1 =&(tmp[i*matchP->mask.nx]);
		freadBS(tmp1,sizeof(unsigned char),matchP->mask.nx,fp,BYTEFLAG);
		matchP->mask.m[i]=tmp1;
	}
	fclose(fp);
}

 
static void readArgs(int argc,char *argv[],char **earlyFile,char **lateFile, matchParams *matchP)
{
	int32 n,i,idum;
	char *argString;

	matchP->chipSize=33;
	matchP->searchX=MINSEARCH;	       matchP->searchY=MINSEARCH;
	matchP->stepX=DEFAULTSTEP;	       matchP->stepY=DEFAULTSTEP;
	matchP->maskFile=NULL; 
	matchP->mask.m=NULL;
	matchP->velFile=NULL;
	matchP->velocity.vx=NULL;
	matchP->velocity.vy=NULL;
	matchP->slowFlag=FALSE;
	matchP->proj=-1;
	if( argc < 3 || argc > (4+15) ) LSusage();        /* Check number of args */
	*lateFile = argv[argc-2];
	*earlyFile = argv[argc-3];
	matchP->outputFile = (char *) malloc(strlen(argv[argc-1])+1);
	strcpy(	matchP->outputFile, argv[argc-1]);

	n = argc - 4;
	fprintf(stderr,"%i %i\n",argc,n);
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');
		if(strstr(argString,"chipSize") != NULL) {
			sscanf(argv[i+1],"%i",&idum);  
			if( (idum % 2) == 0) idum++; /* ensure its odd */
			matchP->chipSize=idum;
			i++;
		} else  if(strstr(argString,"maskfile") != NULL) {
			matchP->maskFile=argv[i+1];
			i++;
		}else  if(strstr(argString,"velfile") != NULL) {
			matchP->velFile=argv[i+1];
			i++;
		} else if(strstr(argString,"epsg") != NULL) {
			sscanf(argv[i+1],"%i",&idum);  
                        if(idum !=NSIDCNORTH && idum != NSIDCSOUTH && idum != NULLPROJ ) error("invalid projection, choices are 3413 and 3031 or 9999 (null)\n");
			matchP->proj=idum;
			i++;
		} else if(strstr(argString,"searchRadius") != NULL) {
			sscanf(argv[i+1],"%i",&idum);  
			matchP->searchX=idum;
			matchP->searchY=idum;
			i++;
		} else if(strstr(argString,"stepX") != NULL) {
			sscanf(argv[i+1],"%i",&idum);  
			matchP->stepX=idum;
			i++;
		} else if(strstr(argString,"stepY") != NULL) {
			sscanf(argv[i+1],"%i",&idum);  
			matchP->stepY=idum;
			i++;
		} else if(strstr(argString,"slow") != NULL) {
			matchP->slowFlag=TRUE;
		} else LSusage();
	}
	return;
}

static void LSusage()
{
	error("\033[1m\n\n %s\n %s \n\t%s \n\t%s \n\t%s \n\t%s \n\t%s \n\t%s \n\t%s \n\t%s\n\t%s \n\t%s \n",
	      "lstrack -searchRadius -slow searchRadius -velfile velfile -maskfile maskfile -epsg epsg  \\\n\t-chipSize chipSize-stepX stepX -stepY stepY fileEarly  fileLate outfile \n",
	      "where \n",
	      "slow =\t\tincrement chipSize (+50%%) and reduce search radius to min",
	      "searchRadius =\tsearch this distance in x,y for a match ",
	      "maskfile  =\tbinarymask file to skip areas of water etc  ",
	      "velfile  =\troot for velfile.vx, velfile.vy used for initial guess  ",
	      "epsg  =\t\tprojection 3013 for south, 3413 for north  ",
	      "chipSize =\tsize for the search chip used in a match ",
	      "stepX,stepY =\tcompute a match every stepX by stepY pixels in image",
	      "fileEarly =\tfirst image in temporal sequence",
	      "fileLate =\tsecond image in temporal sequence",
	      "outfile =\troot name for output (outfile.offX, outfile.offY, outfile.rho,outfile.dat)\033[0m");

}
