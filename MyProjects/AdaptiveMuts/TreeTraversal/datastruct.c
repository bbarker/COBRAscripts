#include <stdlib.h>
#include <stdio.h>
#include "datastruct.h"

//Example and test case: prints 1:10 backwards
/*
int main(void) {
  StackArray* stack;
  stack = stackInit(5,5);
  //set RESIZEINC to 5 while testing
  int AnArray[20];
  int i,*x;
  for (i=0; i < 11; i++) {
    AnArray[i]=i;
    stackPush(stack,&AnArray[i]);
  }
  while(stack->top > 0) {
    x = (int*) stackPop(stack);
    printf("%d\n",*x);
  }
  stackDestroy(stack);
  return; 
}
*/

StackArray* stackInit(long initsize, long incsize) {
  void *tmpvoid;
  tmpvoid = (void*)&initsize;
  int psz = sizeof(tmpvoid);  //Better way to get this?
  StackArray* stack;
  stack = malloc(sizeof(StackArray));
  stack->top = -1;
  stack->size = initsize;
  stack->psize = psz;
  stack->incsize = incsize;
  stack->darray = NULL;
  //Next line generates valgrind warning, probably 
  //because we haven't done anything with darray.
  //Shouldn't be an issue. Except it is in OpenMP.
  //So initilize to NULL above
  stack->darray = realloc(stack->darray, psz*initsize);
  return stack;  
}


//Keep track of currentsize externally
void stackChangeSize(StackArray* stack, long sizechg) {
  printf("STACK psize: %d, size: %ld, sizechg: %ld\n",
         stack->psize,stack->size,sizechg);
  long stackbytes = (stack->size+sizechg) *
                    (long)stack->psize;
  printf("stack bytes: %ld\n",stackbytes);
  void** tempstack;
  tempstack = realloc(stack->darray, stackbytes);
  if (tempstack == NULL) {
    stackDestroy(stack);
    printf("OUT OF MEMORY: tried to allocate %ld pointers.\n", stackbytes);
  } 
  else {
    stack->darray = tempstack;
    stack->size = stack->size + sizechg;
  }
  return;
}

void* stackPop(StackArray* stack) {
  void* tmp;
  tmp = stack->darray[stack->top];
  //stack->darray[stack->top] = 0;
  stack->top--;
  //Don't worry about decreasing size for now
  return tmp; 
}

void stackPush(StackArray* stack, void* item) {
  stack->top++;
  if (stack->top == stack->size) {
    stackChangeSize(stack, stack->incsize);
  }
  stack->darray[stack->top] = item; 
  return;
}

void stackDestroy(StackArray* stack) {
  free(stack->darray);
  free(stack);
}





