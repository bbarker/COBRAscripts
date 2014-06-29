/* Brandon Barker - Dec 01, 2011
 * Requires arg1 arg2 arg3
 * Data file is arg1
 * The file must have 2^arg2 columns and arg3 rows. */

/* TODO: 
 Program Issues:
 1) Fix weird problem for opening concatenation of 4 10x1000 files.
 2) Fix parallel seg fault (on e.g. 4+ threads).
  For 1, consider:
  http://stackoverflow.com/questions/2279052/increase-stack-size-in-linux-with-setrlimit
  Alternatively, incease the stack size on the command line by using ulimit -s NUMBKBYTES

  This very likely is also causing the parallel issue (2) since threads may have smaller
  stack sizes than a unithread process. 

 Goals:
 0) !! Need to do subset checking when counting terminations, suboptimals, etc.
  -- can probably do this simply by bitwise comparison
 1) Explicitly count the number of paths cut off by various traps.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "checkPaths.h"
#include "datastruct.h"

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

/* Returns 1 if the binary digit in v is 1, *
 * otherwise, returns 0.                    */
short getBinaryDigit(bitvec v, short digit) {
  short digval;
  bitvec mask = 0;
  mask = ~(~mask << 1) << digit;
  digval =  (short) ((v & mask) >> digit);
  return digval;
}

//Probably could speed this up by directly using shifts
short binDigitSum(bitvec v, short numdigits) {
  short digsum = 0;
  short i;
  for (i = 0; i <= numdigits; i++) {
    digsum += getBinaryDigit(v,i);
  }
  return digsum;
}
bitvec binDigitsToDec(short numidx, short digit[numidx]) {
  bitvec bv = 0;
  int i;
  for (i=0; i<numidx; i++) {
    bitvec mask = 0;
    bv |= ~(~mask << 1) << (digit[i]); //Assumes min is 0
  }
  return bv;
}

short MSB(bitvec v) {
  short r=0;
  while (v >>= 1) //unroll for more speed
  { //Taken from Sean Eron Anderson Bit Twiddling page
    r++; //Faster methods available on same page.
  }
  return r;
}

/* Need to implement a stack now for DFS-style traversal. 
 * We can then logical-OR each bit with 1 and add those
 * that increase in value to the stack.  Allow a function 
 * pointer (?) as an argument that can be used to test
 * each "node", thus deciding whether its children should
 * be added to the stack (a dummy would always return true).
 * In cases where more data needs to be stored, need to call
 * free on data.  Node creation is still a problem; easier to
 * just copy and alter template function.
 */

//Barebones example
unsigned long traverseAllPaths(const double* grRates, short ngen) {

  StackArray* pathstack;
  unsigned long finishedpaths = 0;
  unsigned ncols;
  ncols = pow((unsigned) 2, (unsigned) ngen);
  /* The maximum depth is the number of mutants, and since we
   * don't care if we've visited a node, memory requirements
   * should be quite minimal. */
  pathstack = stackInit((long) (ngen*ngen), (long) ngen);
  stackPush(pathstack, (void*) grRates);
  double* curnode;
  bitvec m, nextidx;
  bitvec groffset;

  while (pathstack->top > -1) {
    curnode = (double *) stackPop(pathstack);
    groffset = (bitvec)(curnode-grRates);
    //Do any analysis on curnode here
    if (groffset == (ncols - 1)) {
      finishedpaths += 1;
    } 
    else { // Add other nodes to stack
      for(m=0; m<ngen; m++) {
	nextidx = groffset | ~(~0 << 1) << m;
        //printf("Nextidx: %d\n",nextidx);
	if (nextidx > groffset) {
	  //Add additional checks here.
	  stackPush(pathstack, (void*) (grRates+nextidx));
          /*if (pathstack->top > maxstacksz) {
	    maxstacksz = pathstack->top;
	    printf("stack top: %ld\n", maxstacksz);
          }*/
	}
      }// end m<ngen
    }//end else (add nodes to stack)
  } //end main while loop
  stackDestroy(pathstack);
  return finishedpaths;
} //End function traverseAllPaths

unsigned long printAllPaths(double* grRates, short ngen, FILE* trajectfile) {
  unsigned i,j;
  StackArray* pathstack;
  unsigned long finishedpaths = 0;
  unsigned ncols;
  ncols = pow((unsigned) 2, (unsigned) ngen);
  /* The maximum depth is the number of mutants, and since we
   * don't care if we've visited a node, memory requirements
   * should be quite minimal. */
  pathstack = stackInit((long) (ngen*ngen), (long) ngen*2);
  short drows = 2;
  double* grRatesDub[ncols][drows];
  for(i=0; i<ncols; i++) {
    grRatesDub[i][0] = &grRates[i];
    grRatesDub[i][1] = NULL;
  }
  stackPush(pathstack, (void*) grRatesDub);
  void **curnode, **nextnode;
  double **traceback;
  double revarray[ngen+1];
  bitvec revmutnum[ngen+1];
  bitvec m, nextidx;
  bitvec groffset;
  short bdig;
  bitvec numrowbytes = (bitvec) drows*sizeof(pathstack->psize);
  while (pathstack->top > -1) {
    curnode = (void **) stackPop(pathstack);
    groffset = ((bitvec) ((double *)curnode-(double *)grRatesDub))/drows;
    //Do any analysis on curnode here
    if (groffset == (ncols - 1)) {
      finishedpaths += 1;
      traceback = (double **) curnode;
      //Reverse this by using temp array.
      for (i=0; i<=ngen; i++) {
        revarray[ngen-i] = **traceback;
        revmutnum[ngen-i] = ((bitvec) ((double *)traceback-(double *)grRatesDub))/drows;
        traceback = (double **) *(traceback+1);
      }
      for (i=0; i<=ngen; i++) {
        if (i == 0) {
          fprintf(trajectfile,"%ld\t",revmutnum[i]); 
	}
        else {
          bdig = MSB(revmutnum[i]^revmutnum[i-1]);
          fprintf(trajectfile,"%d\t", bdig+1);
        }
      }
      fprintf(trajectfile,"\n");
      for (i=0; i<=ngen; i++) {
        fprintf(trajectfile,"%lf\t",revarray[i]/revarray[0]);
      }
      fprintf(trajectfile,"\n");
    } 
    else { // Add other nodes to stack
      for(m=0; m<ngen; m++) {
	nextidx = groffset | ~(~0 << 1) << m;
        nextnode = (void**) (((double *)grRatesDub) + (nextidx*drows));
        //long diff2 = 0;
        //diff2 = (long)((void*)nextnode-(void*)grRatesDub);
   
	// printf("NextIDX: %d %d %ld\n", groffset, nextidx, diff2);

        //nextnode = (void*) &(grRatesDub[nextidx][0]);
        //printf("Nextidx: %d\n",nextidx);
        //printf("%d ",m);
	if (nextidx > groffset) {
	  //Add additional checks here.
	  //printf("%p %p\n", grRatesDub[ncols-1][0], curnode);
	  nextnode[1]=&curnode[0];
	  stackPush(pathstack, (void *) nextnode);
	}
      }// end m<ngen
      //break;
    }//end else (add nodes to stack)
  } //end main while loop
  stackDestroy(pathstack);
  return finishedpaths;
} //End function traverseAllPaths
 

pathstats countTraps(const double* grRates, short ngen, float s) {
  /* First, find the first occurence of the maximum value 
   * Then, find all nodes x_i such that: x_i < x_m = max, and
   * x_(i-1) < x_i > x_(i+1) and i < n.  Should make fuzzy:
   * x_(i-1)/x_i < 1-s > x_(i+1)/x_i and x_i/x_m + s  < 1. 
   * The problem here is that we only check for abrupt local 
   * maxima; need to extend this: make previous the last node
   * that was different. 
   * It may be better to use x_i/x_m < 0.99 if s sufficiently large*/

  unsigned ncols;
  ncols = pow((unsigned) 2, (unsigned) ngen);
  pathstats pstats;

  pstats.finished = 0;
  pstats.terminated = 0;
  pstats.subopt = 0;
  double maxf = 0;
  pstats.maxfpos = 0;
  unsigned long i;
  for (i=0; i<ncols; i++) {
    if(grRates[i] > maxf) {
      maxf = grRates[i];
      pstats.maxfpos = i;
    }  
  }
  printf("Maxf: %ld %f\n", pstats.maxfpos, maxf);
  /* Need to create a struct with a ptr to the entry in grRates,
   * and a ptr to the previous entry. */

  /* Actually if we really only need two values, may as well just
   * copy them all in to a new Nx2 array  and save some headaches */

  /* !!!!!!! But does this work? is it possible the same node 
   * could be overwritten? Shouldn't be; DFS style search.*/
  short drows = 2;
  double grRatesDub[pstats.maxfpos+1][drows];
  for(i=0; i<=pstats.maxfpos; i++) {
    grRatesDub[i][0] = grRates[i];
    grRatesDub[i][1] = -1.0;
  }
  StackArray* pathstack;
  pathstack = stackInit((long) (ngen*(ngen-1)), (long) 2*ngen);
  stackPush(pathstack, (void*) grRatesDub);  
  double *curnode, *nextnode;
  bitvec m, nextidx;
  bitvec groffset;

  unsigned long counterw = 0;
  //bitvec numrowbytes = (bitvec) drows*sizeof(pathstack->psize);
  bitvec numrowbytes = (bitvec) drows*sizeof(double);
  while (pathstack->top > -1) {
    //if (counterw > 10) {break;} counterw++;
    curnode = (double *) stackPop(pathstack);
    groffset = (bitvec)(((double*)curnode-(double*)grRatesDub)/drows);
    //groffset = ( (bitvec) ((void*)curnode-(void*)grRatesDub)) / drows;
    printf("groffset %d\n", groffset);
    //Do any analysis on curnode here
    if (groffset == pstats.maxfpos) {
      pstats.finished += 1;
    } 
    else { // Add other nodes to stack
      for(m=0; m<ngen; m++) {
	nextidx = groffset | ~(~0 << 1) << m;
        nextnode = ((double *)grRatesDub)+nextidx*drows;
        printf("nextidx: %d %ld\n", nextidx, (long) (nextnode-(double *)grRatesDub));
	if (nextidx > groffset && nextidx <= pstats.maxfpos) {
	  //Add additional checks here.
          if (curnode[1] >= 0 &&
              curnode[0]/maxf + s < 1 && 
              curnode[1]/curnode[0] < 1-s &&
              nextnode[0]/curnode[0] < 1-s) {
            pstats.terminated += 1; //Save these to a file?
            //printf("%lf %lf %lf %lf\n", curnode[1], curnode[0], nextnode[0], maxf);
            /* Note a terminal path could embed other local optima,
             * so not all are currently traced by this method */
          }
          else {
            //Set the prior value of the next node
	    //to the most recent different value.
 	    if (curnode[0] != nextnode[0]) {
              nextnode[1] = curnode[0];
            }
            else {
	      nextnode[1] = curnode[1];
            }
	    stackPush(pathstack, (void*) nextnode);
          }
	}
        else if (nextidx > pstats.maxfpos && nextnode[0]/curnode[0] > 1+s) {
          pstats.subopt += 1;
        }
      }// end m<ngen
    }//end else (add nodes to stack)
  } //end main while loop
  stackDestroy(pathstack);
  return pstats;
} //End function countTraps


//TODO: Read Function Pointers in the Wild:
//http://www.cprogramming.com/tutorial/function-pointers.html

void makeOutFileName(char* outfilename, const char* infilename, const char* suffix) {
  strcpy(outfilename, infilename);
  char *fname, *tmpfname;
  tmpfname = strtok(outfilename, "\\/");
  while (tmpfname != NULL) {
    fname = tmpfname;
    tmpfname = strtok(NULL, "\\/");
  }
  fname = strtok(fname,".");
  strcpy(outfilename,fname);
  strcat(outfilename,suffix);
  return;
}

void readMutFile(const char* infilename, int nrows, int ncols, double grRatesMat[][ncols]) {
  int i,j;
  FILE *infile;
  infile = fopen(infilename,"r");
  if (infile == NULL) {
    printf("Couldn't open infile.\n");
    exit(-2);
  }
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      //Consider doing the following one row at a time
      //if memory becomes an issue.
      fscanf(infile, "%lf,", &grRatesMat[i][j]);
    } //end j < ncols

  } // end j < nrows
  fclose(infile);
  return;
}

/**************** Node functions *************************/

short returnTrue(const void** node) {
  return 1;
}
