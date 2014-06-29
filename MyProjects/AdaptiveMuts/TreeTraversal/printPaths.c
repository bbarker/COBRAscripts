#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "checkPaths.h"

int main(int argc, const char **argv) {
  int i, j;

  if (argc < 4) {
    printf("Not enough arguments\n");
    exit(-1);
  }
  int ngen = atoi(argv[2]);
  printf("Num Genes: %d\n", ngen);
  unsigned int ncols = (unsigned int) pow((unsigned int) 2, (unsigned int) ngen);
  printf("Num Mutants: %d\n", ncols);
  int nrows = atoi(argv[3]);
  double grRatesMat[nrows][ncols];
  readMutFile(argv[1], nrows, ncols, grRatesMat);

  char outfilename[100]; //Better way?
  makeOutFileName(outfilename, argv[1], "_paths.csv");
  printf("%s\n",outfilename);
  FILE *outfile;
  outfile = fopen(outfilename,"w");
  if (outfile == NULL) {
    printf("Couldn't open outfile.\n");
    exit(-2);
  } 
  unsigned long finishedPaths[nrows];
  for (i=0; i<nrows; i++) {
    //printf("Starting row %d.\n",i);
    fprintf(outfile, "\n*** Replicate %d ***\n", i);
    finishedPaths[i] = printAllPaths(grRatesMat[i], ngen, 
                       outfile);
    //finishedPaths[i] = traverseAllPaths(grRatesMat[i], ngen);
  }
  /*for (i=0; i<nrows; i++) {
    printf("%ld ", finishedPaths[i]);
  }*/

  fclose(outfile);
  //Later on, it may be good to analyze trap position as well.
  //Best to wait until parallel bugs are fixed to be sure. 
  return 0;
}
