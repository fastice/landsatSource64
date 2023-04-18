/*
  Constants
*/
#define CONSTANTFIT 1
#define PLANEFIT 3
/*
Structures
*/
typedef struct {
	unsigned int fitType;
	double *lat;
	double *lon;
	double *x;
	double *y;
	double *z; /* probably not needed but included in most tie files */
	double *vx;
	double *vy;
	double *offXT;
	double *offYT;
	double *xyScale;
	double deltaT;
	int32_t npts;
} lsTiepoints;
#define MAXFITPARAM 3
typedef struct {
	char *tieFile;
	char *matchFile;
	char *fitFile;
	unsigned int fitType;
	int32_t slowFlag;
	int32_t proj;
	double deltaT;
	double pX[MAXFITPARAM];
	double pY[MAXFITPARAM];
	double sigmaXRes; /* Residual from tie fit in meters*/
	double sigmaYRes;
	double Cx[MAXFITPARAM][MAXFITPARAM];
	double Cy[MAXFITPARAM][MAXFITPARAM];
} lsFit;
/*
   Functions
*/
void readLSTiePoints(char *tieFile, lsTiepoints *tiePoints);
void readLSOffsets(lsFit *fitDat, matchResult *matches,int32_t readData,char *mask);
char *stripWS(char *string);
void interpLSTies( lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints);
void parseKeyValue(char *line,char *keyword,char *value);
double xyscale(double latctr,int32_t proj);
void computeLSFit( lsFit *fitDat, matchResult *matches, lsTiepoints *tiePoints);
