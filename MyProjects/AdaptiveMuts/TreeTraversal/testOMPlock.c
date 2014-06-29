#include <stdio.h>
#include <omp.h>

int main(void) {

  int parasum = 0;
  printf("%d\n",parasum);

  omp_lock_t writelock;
  omp_init_lock(&writelock);

  int i;
  #pragma omp parallel for
  for ( i = 0; i < 10; i++ )
  {
    int tmp = 2*i;
    omp_set_lock(&writelock);
    // one thread at a time stuff
    parasum += tmp;
    omp_unset_lock(&writelock);
    // some stuff
  }

  omp_destroy_lock(&writelock);

  printf("%d\n",parasum);

}
