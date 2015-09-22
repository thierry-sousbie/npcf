#ifndef _CORREL_THREADS_HEADER_H__
#define _CORREL_THREADS_HEADER_H__

#include "bbtree.h"
#include "correl_parms.h"

typedef struct threaded_subDivide_parms_str {
    subDivide_parms parms;
    BBTreeNode **node;
    int ret;
} threaded_subDivide_parms;

int threaded_subDivide(subDivide_parms *parms);

#endif
