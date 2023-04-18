#include <stdio.h>
#include <stdint.h>
/*
  Constants
*/
#define OSSIZE 8            	 /* Half width of oversample window - i.e. oversample correlation with 2*OSSIZE */
#define MINSEARCH 5	  /* Minimum search window size */
#define NOVER 16	  	  /* Oversampling factor */
#define MAXW 1001      	 /* Maxsize for correlation calculation, should be ridiculously large */
#define RHOTHRESH 0.15  /* Don't save result that has correlation below this value */
#define RHOTHRESHFINAL 0.5  /* Don't save result that has correlation below this value */
#define NODATA -2.0e9
#define NOCORR -1.0
#define NOUSEABLEDATA -2.0
#define DEFAULTCHIP 33
#define DEFAULTSTEP 32

#define NSIDCNORTH 3413
#define NULLPROJ 9999
#define ROTNORTH 45.0
#define SLATNORTH 70.0
#define ROTSOUTH 0.0
#define SLATSOUTH 71.0

#define NSIDCSOUTH 3031
#define KMTOMS 1000
#define MTOKM 0.001
#define TRUE 1
#define FALSE 0
/* Match type codes */
#define NODATATOUSE 0x1
#define LOWCORR 0x2
#define ISMASKED 0x3
#define SLOWVEL 0x20
#define MIDVEL 0x40
#define FASTVEL 0x80
#define HIGHSTRAIN 0x10
#define SLOWTHRESH 100
#define MIDTHRESH 1000
#define SIGMATHRESH 50
typedef struct {
	uint32_t nx;    	/*  Used to input size of area, defaults whole image */
	uint32_t ny;
	double x0;
	double y0;
	double dx;
	double dy;
	unsigned char **m;
} matchMask;

typedef struct {
	uint32_t nx;    	/*  Used to input size of area, defaults whole image */
	uint32_t ny;
	double x0;
	double y0;
	double dx;
	double dy;
	float  **vx;
	float  **vy;
} velMap;


typedef struct {
	char *maskFile;
	char *velFile;
	char *outputFile;
	int32_t slowFlag;
	int32_t proj;
	uint32_t nx;    	/*  Used to input size of area, defaults whole image */
	uint32_t ny;
	uint32_t stepX;
	uint32_t stepY;
	uint32_t deltaX;  /* Offset from origin to start matching at */
	uint32_t deltaY;  /* Offset from origin to start matching at */
	uint32_t searchX; /* Default size of area to search */
	uint32_t searchY;
	uint32_t chipSize;
	float deltaT;
	matchMask mask;
	velMap velocity;
} matchParams;

typedef struct {
	char *sat;
	int32_t path;
	int32_t row;
	int32_t year;
	int32_t doy;
	double jd;
	uint32_t nx;
	uint32_t ny;
	double x0;
	double y0;
	double dx;
	double dy;
	int32_t proj;
	uint16_t **image;
	float **fimage;
} LSimage;

typedef struct {
	char *fileEarly;
	char *fileLate; 
	double jdEarly;
	double jdLate;
	uint32_t nx;
	uint32_t ny;
	double x0;
	double y0;
	double dx;
	double dy;
	double meanSigmaX;
	double meanSigmaY;
	int32_t stepX;
	int32_t stepY;
	float **X;
	float **Y;
	float **sigmaX;
	float **sigmaY;
	float **Rho;
	uint8_t **type;
} matchResult;


