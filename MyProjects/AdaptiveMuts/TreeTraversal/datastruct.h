#ifndef __DATASTRUCT__
#define __DATASTRUCT__

typedef struct StackArray {
  void** darray;
  int psize;
  long top;
  long size;
  long incsize;
} StackArray;

StackArray* stackInit(long initsize, long incsize);
void stackChangeSize(StackArray* stack, long sizechg);
void* stackPop(StackArray* stack);
void stackPush(StackArray* stack, void* item);
void stackDestroy(StackArray* stack);


//__DATASTRUCT__
#endif 


