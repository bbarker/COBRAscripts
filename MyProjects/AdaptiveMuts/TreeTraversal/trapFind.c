#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "trapFind.h"

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


  FILE *outfile;
  char outfilename[100]; //Better way?
  makeOutFileName(outfilename, argv[1], "_trapAnalysis.csv");
  
  pathstats PathStats[nrows];
  #pragma omp parallel for
  for (i=0; i<nrows; i++) {
    //printf("Starting row %d.\n",i);
    PathStats[i] = countTraps(grRatesMat[i], ngen, 0.01);
  }
  //For small sets print all trajectories, not thread safe  
  /*  if (ngen < 6) { 
    char trajfilename[100];
    makeOutFileName(outfilename, argv[1], "_trajectories.csv");
    }*/ //Save for a different program
  outfile = fopen(outfilename,"w");
  if (outfile == NULL) {
    printf("Couldn't open outfile.\n");
    exit(-2);
  }
  short nmutgen = 0;
  for (i=0; i<nrows; i++) {
    nmutgen = binDigitSum((bitvec) PathStats[i].maxfpos, ngen); 
    fprintf(outfile, "%ld,%ld,%ld,%ld\n", (unsigned long) nmutgen,
            PathStats[i].finished, PathStats[i].terminated,
	    PathStats[i].subopt);
  }
  fclose(outfile);
  //Later on, it may be good to analyze trap position as well.
  //Best to wait until parallel bugs are fixed to be sure. 
  return 0;
}
