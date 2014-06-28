#ifndef __CHECKPATHS__
#define __CHECKPATHS__


//short should be ok for up to 16 digits on this system
//change to #define later
typedef unsigned long bitvec;

//Consider making this conditionally based on ngen
//Or make it a fixed size? (investigate C standard)
short getBinaryDigit(bitvec v, short digit);
bitvec binDigitsToDec(short numidx, short digit[numidx]);
short binDigitSum(bitvec v, short numdigits);
short MSB(bitvec v);

typedef struct pathstats {
  unsigned long maxfpos;
  unsigned long finished;
  unsigned long terminated;
  unsigned long subopt;

} pathstats;


/****************Traverasl Functions*******************/
/*unsigned long traverseAllPaths(const double* grRates, short nmuts, short callfree,  short (*nodeEval) (const void**));*/

unsigned long traverseAllPaths(const double* grRates, short nmuts);

unsigned long printAllPaths(double* grRates, short nmuts, FILE* trajectfile);

pathstats countTraps(const double* grRates, short ngen, float s);


/**************** Node functions *************************/

short returnTrue(const void** node);


/**************** Util **********************************/

void makeOutFileName(char* outfilename, const char* infilename, const char* suffix);
//Seems odd that we can't specify both dimensions of the "matrix" below. 
//Try to understand this.
void readMutFile(const char* infilename, int nrows, int ncols, double grRatesMat[][ncols]);
//__CHECKPATHS__
#endif 


