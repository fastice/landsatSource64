
#include "stdio.h"
#include <stdint.h>
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
   void warning( FILE *fp, char *fmt, ... );     

/* getline: read a line, return length */
   int32_t getlineobsolete(char *line, int32_t imax);  

/* fgetline: read a line from file, return length */
   int32_t fgetline(FILE *fp, char *line, int32_t imax);  

/* 
   Search through file for next occurence of BEGINNINGOFDATA n, 
   where n = number of columns. Return line count.
*/

   int32_t goToBeginningOfData(FILE *fp, int32_t lineCount, int32_t *ncolumns);


/* Search through file for next occurence of ENDOFDATA */

   int32_t goToEndOfData(FILE *fp, int32_t lineCount, int32_t *ndata);

/* Open file for input, exit if error occurs */
   FILE *openInputFile(char *filename);
/* 
   Go to next occurence of "string" in file, return negative if EOF.
*/
 int32_t goToNextStringOccurence(FILE *fp, int32_t lineCount, char *string, char *line);
/* 
  Read through file returning data strings.
  If end of data eod = TRUE else eod = FALSE.
*/
   int32_t getDataString(FILE *fp, int32_t lineCount, char *line,int32_t *eod);
/*
  Read/write raw data and byteswap if needed
*/ 
  size_t freadBS(void * ptr, size_t nitems, size_t size, FILE *fp, int32_t flags);
  size_t fwriteBS(void * ptr, size_t nitems, size_t size, FILE *fp, int32_t flags);
