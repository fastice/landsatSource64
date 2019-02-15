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
	uint32 nx;    	/*  Used to input size of area, defaults whole image */
	uint32 ny;
	double x0;
	double y0;
	double dx;
	double dy;
	unsigned char **m;
} matchMask;

typedef struct {
	uint32 nx;    	/*  Used to input size of area, defaults whole image */
	uint32 ny;
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
	int32 slowFlag;
	int32 proj;
	uint32 nx;    	/*  Used to input size of area, defaults whole image */
	uint32 ny;
	uint32 stepX;
	uint32 stepY;
	uint32 deltaX;  /* Offset from origin to start matching at */
	uint32 deltaY;  /* Offset from origin to start matching at */
	uint32 searchX; /* Default size of area to search */
	uint32 searchY;
	uint32 chipSize;
	float deltaT;
	matchMask mask;
	velMap velocity;
} matchParams;

typedef struct {
	char *sat;
	int32 path;
	int32 row;
	int32 year;
	int32 doy;
	double jd;
	uint32 nx;
	uint32 ny;
	double x0;
	double y0;
	double dx;
	double dy;
	int32 proj;
	uint16 **image;
	float **fimage;
} LSimage;

typedef struct {
	char *fileEarly;
	char *fileLate; 
	double jdEarly;
	double jdLate;
	uint32 nx;
	uint32 ny;
	double x0;
	double y0;
	double dx;
	double dy;
	double meanSigmaX;
	double meanSigmaY;
	int32 stepX;
	int32 stepY;
	float **X;
	float **Y;
	float **sigmaX;
	float **sigmaY;
	float **Rho;
	uint8 **type;
} matchResult;


