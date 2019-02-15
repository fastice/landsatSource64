
#include "stdio.h"
/* #include "malloc.h"*/

#define TRUE 1
#define FALSE 0
#define LINEMAX 10000
#define MAXCOL 500
#define BEGINDATA '#'
#define ENDDATA '&'
#define COMMENT ';'
#define PI 3.141592653589793238
#define DTOR  0.0174532930056254
#define RTOD 57.295777918682049
#define INT16FLAG 1
#define INT32FLAG 2
#define FLOAT32FLAG 3
#define BYTEFLAG 4

#define DELTA 0.0000001
#define CLIGHT 2.99792458e8
#define E0 8.841941e-12
#define U0 1.256637e-6
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
/*
    Include file for library of standard routines
*/



/* Print error and exit */
   void error( char *fmt, ... );     

/* Print warning and write to log file if *fp !=NULL */
   void warning( FILE *fp,char *fmt, ... );     

/* getline: read a line, return length */
   int getlineobsolete(char *line, int imax);  

/* fgetline: read a line from file, return length */
   int fgetline(FILE *fp, char *line, int imax);  

/* 
   Search through file for next occurence of BEGINNINGOFDATA n, 
   where n = number of columns. Return line count.
*/

   int goToBeginningOfData(FILE *fp, int lineCount, int *ncolumns);


/* Search through file for next occurence of ENDOFDATA */

   int goToEndOfData(FILE *fp, int lineCount, int *ndata);

/* Open file for input, exit if error occurs */
   FILE *openInputFile(char *filename);
/* 
   Go to next occurence of "string" in file, return negative if EOF.
*/
 int goToNextStringOccurence(FILE *fp, int lineCount,char *string,char *line);
/* 
  Read through file returning data strings.
  If end of data eod = TRUE else eod = FALSE.
*/
   int getDataString(FILE *fp, int lineCount, char *line,int *eod);
/*
  Read/write raw data and byteswap if needed
*/ 
  size_t freadBS(void * ptr, size_t nitems, size_t size, FILE *fp,int flags);
  size_t fwriteBS(void * ptr, size_t nitems, size_t size, FILE *fp,int flags);
