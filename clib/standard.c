/*#include "stdio.h"*/
#include "stdarg.h"
#include "string.h"
#include "standard.h"
/*#include "geotiff/xtiffio.h"   for TIFF */
#include <stdlib.h>
/* #include <libkern/OSByteOrder.h> */
#
#define SWAP_SHORT(Var)  Var = *(int16_t*) SwapEndian((void*)&Var, sizeof(int16_t))
#define SWAP_USHORT(Var) Var = *(uint16_t*) SwapEndian((void*)&Var, sizeof(uint16_t))
#define SWAP_INT32(Var)   Var = *(int32_t*) SwapEndian((void*)&Var, sizeof(int32_t))
#define SWAP_ULONG(Var)  Var = *(uint32_t*) SwapEndian((void*)&Var, sizeof(uint32_t))
#define SWAP_RGB(Var)    Var = *(int32_t*) SwapEndian((void*)&Var, 3)
#define SWAP_FLOAT(Var)  Var = *(float*) SwapEndian((void*)&Var, sizeof(float))
#define SWAP_DOUBLE(Var) Var = *(double*) SwapEndian((void*)&Var, sizeof(double))

extern void *SwapEndian(void* Addr, const int32_t Nb);

/* 
    Useful procedures that can be used for all programs
*/

/*
     error -- prints error message with arguements "..." with format "fmt" 
              and then exits program. From K&R.
*/

    void error( char *fmt, ... )
{
    va_list args;

    va_start(args,fmt);
    fprintf(stderr,"error:  ");
    vfprintf(stderr,fmt,args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(1);
}

    void warning(FILE *fp, char *fmt, ... )
{
    va_list args;

    va_start(args,fmt);
    fprintf(stderr,"\nWarning:  ");
    if(fp !=NULL) fprintf(fp,"\nWarning:  ");
    vfprintf(stderr,fmt,args);
    if(fp !=NULL) vfprintf(fp,fmt,args);
    fprintf(stderr, "\n");
    if(fp !=NULL) fprintf(fp, "\n");
    va_end(args);
}


/* getline:  read a line, return length.  From K&R */

    int32_t getlineobsolete(char *line, int32_t imax)
{
    if(fgets(line,imax,stdin) == NULL)
        return 0;
    else
        return strlen(line);
}
 
/* fgetline:  read a line from a file, return length. */

    int32_t fgetline(FILE *fp, char *line, int32_t imax)
{
    if( feof(fp) == EOF ) error("%s","fgetline -- EOF encountered");
    if( ferror(fp) != 0 ) error("%s","fgetline -- file error");

    if(fgets(line,imax,fp) == NULL)
        return 0;
    else
        return strlen(line);
}


/* 
   Search through file for next occurence of BEGINDATA n, 
   where n = number of column. Return line count
*/

   int32_t goToBeginningOfData(FILE *fp, int32_t lineCount, int32_t *ncolumns)
{
    int32_t lineLength;  
    char line[LINEMAX];
    char *string;

    while( TRUE ) {  /* Loop to read lines */
 
        lineLength = fgetline(fp,line,LINEMAX); /* Read line */
        lineCount++;
        if(lineLength > 0) {
            if( strchr(line,ENDDATA) != NULL ) /* Read data */
                error("%s","goToBeginningOfData -- & with no matching # "); 
            else if( (string = strchr(line,BEGINDATA) ) != NULL ) { 
                if( sscanf((string+1),"%i",ncolumns) != 1) /*If BEGINDATA,get */
                    error("%s",
                        "goToBeginningOfData -- # found without n columns ");
                return lineCount;   /* Found so exit */
            }  
        }
    }
}


/* 
  Read through file returning data strings.
  If end of data eod = TRUE else eod = FALSE.
*/

   int32_t getDataString(FILE *fp, int32_t lineCount, char *line,int32_t *eod)
{
    int32_t lineLength;  
    int32_t ncolumns;    
    int32_t count;
    count = 0;
    while( TRUE ) {  /* Loop to read lines */
 
        lineLength = fgetline(fp,line,LINEMAX); /* Read line */
        lineCount++;
        count++;
        if( strchr(line,COMMENT) == NULL &&   
         strpbrk(line,
          "0123456789.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
            != NULL ) {
             *eod = FALSE;
             return lineCount; /* Data */ 
         }
         
         if( strchr(line,ENDDATA) != NULL &&
             strchr(line,COMMENT) == NULL ) { /* End of data */
             *eod = TRUE;
             return lineCount;
         }
         if( count >= 100000) error("getDataString: Missing end of data");
    }
}


/* Search through file for next occurence of ENDDATA */

   int32_t goToEndOfData(FILE *fp, int32_t lineCount, int32_t *ndata)
{
    int32_t lineLength;  
    int32_t ncolumns;
    char line[LINEMAX];

    while( TRUE ) {  /* Loop to read lines */
 
        lineLength = fgetline(fp,line,LINEMAX); /* Read line */
        lineCount++;
        if( strchr(line,COMMENT) == NULL &&   
            strpbrk(line,"0123456789.EDed") != NULL ) (*ndata)++; /*Count data*/

        if( strchr(line,ENDDATA) != NULL ) return lineCount;
    }
}

/* Search through file for next occurence of ENDDATA */

   int32_t goToNextStringOccurence(FILE *fp, int32_t lineCount, char *string,char *line)
{
    int32_t lineLength;  
    int32_t ncolumns;

    while( TRUE ) {  /* Loop to read lines */
        lineLength = fgetline(fp,line,LINEMAX); /* Read line */
        lineCount++;
        if(feof(fp) != 0) return -1;
        if( strstr(line,string) != NULL) return lineCount;
    }
}


/* Open file for input, exit if error occurs */

   FILE *openInputFile(char *filename)
{
    FILE *fp;

    if( (fp = fopen(filename,"r")) == NULL )  /* Error if unable to open file */
       error("\n\n%s |%s|\n","openInputFile -- Unable to open file: ",filename);
    return fp;
}


/*
  Read raw data and byteswap if needed
*/ 
  size_t freadBS(void *ptr, size_t nitems, size_t size, FILE *fp, int32_t flags)
{
  int16_t  *ptr16;
  long long i;
  size_t retVal;
  uint32_t *ptrLong;
  size_t n;
  /* Use fread to read the data as requested */
   retVal=fread(ptr,nitems,size,fp);
#ifdef i386
  /*
    ++++++++++++++++++++++++++++++start i386++++++++++++++++++++++++++++++++++++++++
  */
 /* 
    No need to swap byte data 
 */
   if(flags == BYTEFLAG) return(retVal);

 /*
    Float 32 data
 */
   if(flags==FLOAT32FLAG) {
      n=(nitems*size)/sizeof(float);
      ptrLong=(uint32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
      return(retVal);
   }

   if(flags==INT32FLAG) {
      n=(nitems*size)/sizeof(long int);
      ptrLong=(uint32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
      return(retVal);
   }
   if(flags==INT16FLAG) {
      n=(nitems*size)/sizeof(short int);
      ptr16=(int16_t *)ptr;
      for(i=0; i < n; i++) ptr16[i]=SWAP_USHORT(ptr16[i]);
      return(retVal);
   }
  /*------------------------------end i386----------------------------------------*/
#else 
#ifdef x86_64

  /*
    ++++++++++++++++++++++++++++++start i386++++++++++++++++++++++++++++++++++++++++
  */
 /* 
    No need to swap byte data 
 */
   if(flags == BYTEFLAG) return(retVal);

 /*
    Float 32 data
 */
   if(flags==FLOAT32FLAG) {

      n=(nitems*size)/sizeof(uint32_t);
      ptrLong=(uint32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
      return(retVal);
   }

   if(flags==INT32FLAG) {
      n=(nitems*size)/sizeof(uint32_t);
      ptrLong=(uint32_t*)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
      return(retVal);
   }
   if(flags==INT16FLAG) {
      n=(nitems*size)/sizeof(int16_t);
      ptr16=(int16_t*)ptr;
      for(i=0; i < n; i++) ptr16[i]=SWAP_USHORT(ptr16[i]);
      return(retVal);
   }
  /*------------------------------end i386----------------------------------------*/
#else
#ifdef powerpc
  /*++++++++++++++++++++++++++++++start powerpc++++++++++++++++++++++++++++++++++++++++*/
  /* This essentially is does nothing since no byteswapping required */  
   return(retVal);
/*------------------------------end powerpc---------------------------------------- */

#else 
  /*++++++++++++++++++++++++++++++start powerpc++++++++++++++++++++++++++++++++++++++++ */
   error("Invalid processor type; compiled with no processor specifications \n");
/*------------------------------end powerpc----------------------------------------*/

#endif 
#endif
#endif
return 0;
}
/*
  Read raw data and byteswap if needed
*/ 
  size_t fwriteBS(void *ptr, size_t nitems, size_t size, FILE *fp,int32_t flags)
{
  uint16_t *ptr16;
  long long i;
  size_t retVal;
  uint32_t *ptrLong;
  size_t n;
/*
    ++++++++++++++++++++++++++++++start i386++++++++++++++++++++++++++++++++++++++++
*/
#ifdef i386
 /*
    Float 32 data
 */
   if(flags==FLOAT32FLAG) {
      n=(nitems*size)/sizeof(float);
      ptrLong=(int32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   }

   if(flags==INT32FLAG) {
      n=(nitems*size)/sizeof(int32);
      ptrLong=(int32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   }

   if(flags==INT16FLAG) {
      n=(nitems*size)/sizeof(int16_t);
      ptr16=(int16_t *)ptr;
      for(i=0; i < n; i++) ptr16[i]=SWAP_USHORT(ptr16[i]);
   }

  /*------------------------------end i386----------------------------------------*/
#else 
#ifdef  x86_64
 /*
    Float 32 data
 */
   if(flags==FLOAT32FLAG) {
      n=(nitems*size)/sizeof(float);
      ptrLong=(uint32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   }

   if(flags==INT32FLAG) {
      n=(nitems*size)/sizeof(int32_t);
      ptrLong=(uint32_t *)ptr;
      for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   }

   if(flags==INT16FLAG) {
      n=(nitems*size)/sizeof(int16_t);
      ptr16=(uint16_t *)ptr;
      for(i=0; i < n; i++) ptr16[i]=SWAP_USHORT(ptr16[i]);
   }

   /*------------------------------end i386 -----------*/
#else
#ifdef powerpc
/*++++++++++++++++++++++++++++++start powerpc++++++++++++++++++++++++++++++++++++++++*/
  /* This essentially is does nothing since no byteswapping required */  
/*------------------------------end powerpc---------------------------------------- */
#else 
/*++++++++++++++++++++++++++++++start powerpc++++++++++++++++++++++++++++++++++++++++ */
   error("Invalid processor type; compiled with no processor specifications \n");
/*------------------------------end powerpc----------------------------------------*/

#endif 
#endif
#endif
   retVal=fwrite(ptr,nitems,size,fp);

#ifdef i386 
   /* Swap back */
   if(flags==FLOAT32FLAG) for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   if(flags==INT32FLAG)   for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   if(flags==INT16FLAG)   for(i=0; i < n; i++) ptr16[i]=SWAP_USHORT(ptr16[i]);
#endif

#ifdef x86_64
   /* Swap back */
   if(flags==FLOAT32FLAG) for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   if(flags==INT32FLAG)   for(i=0; i < n; i++) ptrLong[i]=SWAP_ULONG(ptrLong[i]);
   if(flags==INT16FLAG)   for(i=0; i < n; i++) ptr16[i]=SWAP_USHORT(ptr16[i]);
#endif
   return(retVal);
}





static long _TestEndian=1;

int32_t IsLittleEndian(void) {
	return *(char*)&_TestEndian;
}

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

void *SwapEndian(void* Addr, const int32_t Nb) {
	static char Swapped[16];
	switch (Nb) {
		case 2:	Swapped[0]=*((char*)Addr+1);
				Swapped[1]=*((char*)Addr  );
				break;
		case 3:	// As far as I know, 3 is used only with RGB images
				Swapped[0]=*((char*)Addr+2);
				Swapped[1]=*((char*)Addr+1);
				Swapped[2]=*((char*)Addr  );
				break;
		case 4:	Swapped[0]=*((char*)Addr+3);
				Swapped[1]=*((char*)Addr+2);
				Swapped[2]=*((char*)Addr+1);
				Swapped[3]=*((char*)Addr  );
				break;
		case 8:	Swapped[0]=*((char*)Addr+7);
				Swapped[1]=*((char*)Addr+6);
				Swapped[2]=*((char*)Addr+5);
				Swapped[3]=*((char*)Addr+4);
				Swapped[4]=*((char*)Addr+3);
				Swapped[5]=*((char*)Addr+2);
				Swapped[6]=*((char*)Addr+1);
				Swapped[7]=*((char*)Addr  );
				break;
		case 16:Swapped[0]=*((char*)Addr+15);
				Swapped[1]=*((char*)Addr+14);
				Swapped[2]=*((char*)Addr+13);
				Swapped[3]=*((char*)Addr+12);
				Swapped[4]=*((char*)Addr+11);
				Swapped[5]=*((char*)Addr+10);
				Swapped[6]=*((char*)Addr+9);
				Swapped[7]=*((char*)Addr+8);
				Swapped[8]=*((char*)Addr+7);
				Swapped[9]=*((char*)Addr+6);
				Swapped[10]=*((char*)Addr+5);
				Swapped[11]=*((char*)Addr+4);
				Swapped[12]=*((char*)Addr+3);
				Swapped[13]=*((char*)Addr+2);
				Swapped[14]=*((char*)Addr+1);
				Swapped[15]=*((char*)Addr  );
				break;
	}
	return (void*)Swapped;
}

