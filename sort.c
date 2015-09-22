#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "sort.h"

static unsigned int *global_arr_qsort; // used for quick sort
static FLOAT *global_arr_qsortF;

// used by qsort (in sortarray)
// to sort indexes of an array (no the array itself)
// -> ascending order
int CompareUINT_a(const void *ap, const void *bp)
{
  unsigned int *d=global_arr_qsort;
  unsigned int  a=*((unsigned int*)ap);
  unsigned int  b=*((unsigned int*)bp);
	
  if (d[a]<d[b]) return -1;
  else if (d[a]>d[b]) return 1;
  else return 0;
}

// used by qsort (in sortarray)
// same as CompareUINT_a, descending order
int CompareUINT_d(const void *ap, const void *bp)
{
  unsigned int *d=global_arr_qsort;
  unsigned int  a=*((unsigned int*)ap);
  unsigned int  b=*((unsigned int*)bp);
	
  if (d[a]<d[b]) return 1;
  else if (d[a]>d[b]) return -1;
  else return 0;
}

// returns the array of index that sort array arr
// n is the number of elments in arr
// set reverse!=0 for descending order
unsigned int *SortArray(unsigned int *arr, INT n, int stride_p, int reverse)
{
  unsigned int *index;
  unsigned int i;
  int stride=stride_p;

  if (stride<=1) stride=1;
  
  index = malloc(sizeof(unsigned int)*n);
  for (i=0;i<n;i++)
      index[i]=i*stride;
	
  global_arr_qsort=arr;
	
  if (reverse)
    qsort(index,n,sizeof(unsigned int),CompareUINT_d);
  else
    qsort(index,n,sizeof(unsigned int),CompareUINT_a);
	
  global_arr_qsort=NULL;

  for (i=0;i<n;i++)
      index[i]/=stride;

  return index;
}

// used by qsort (in sortarray)
// to sort indexes of an array (no the array itself)
// -> ascending order
int CompareFLOAT_a(const void *ap, const void *bp)
{
  FLOAT *d=global_arr_qsortF;
  INT a=*((INT*)ap);
  INT b=*((INT*)bp);
	
  if (d[a]<d[b]) return -1;
  else if (d[a]>d[b]) return 1;
  else return 0;
}

// used by qsort (in sortarray)
// same as CompareFLOAT_a, descending order
int CompareFLOAT_d(const void *ap, const void *bp)
{
  FLOAT *d=global_arr_qsortF;
  INT a=*((INT*)ap);
  INT b=*((INT*)bp);
	
  if (d[a]<d[b]) return 1;
  else if (d[a]>d[b]) return -1;
  else return 0;
}

// returns the array of index that sort array arr
// n is the number of elments in arr
// set reverse!=0 for descending order
INT *SortArrayF(FLOAT *arr, INT n, int reverse)
{
  INT *index;
  INT i;
	
  index = malloc(sizeof(INT)*n);
  for (i=0;i<n;i++)
    index[i]=i;
	
  global_arr_qsortF=arr;
	
  if (reverse)
    qsort(index,n,sizeof(INT),CompareFLOAT_d);
  else
    qsort(index,n,sizeof(INT),CompareFLOAT_a);
	
  global_arr_qsortF=NULL;
  	
  return index;
}
