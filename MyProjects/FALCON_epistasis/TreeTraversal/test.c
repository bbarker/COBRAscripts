#include <stdio.h>

short MSB(unsigned v) {
  short r=0;
  while (v >>= 1) //unroll for more speed
  { //Taken from Sean Eron Anderson Bit Twiddling page
    r++; //Faster methods available on same page.
    printf("%d ", v);
  }
  printf("\n");
  return r;
}

int main(void){
  unsigned i;
  for(i=0; i <20; i++) {
    printf("%d: %d\n", i, MSB(i));
  }
}
