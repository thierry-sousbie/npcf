#ifndef __MY_SORT_INTERFACE_H__
#define __MY_SORT_INTERFACE_H__

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

    unsigned int *SortArray(unsigned int *arr, INT n, int stride,int reverse);
    INT *SortArrayF(FLOAT *arr, INT n, int reverse);

#ifdef __cplusplus
}
#endif

#endif
