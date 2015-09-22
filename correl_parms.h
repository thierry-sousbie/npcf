#ifndef __CORREL_PARAMS_H__
#define __CORREL_PARAMS_H__

#include "bbtree.h"

typedef struct subDivide_parms_str {
    int npoints;
    int npairs;

    BBTree **tree;
    BBTree *tmp;
    int *index;
    int tmpi;
    
    double sizetol;
    int *nbins;
    int nbins_tot;
    int nbinsp1_tot;
    double max_size;
    double max_size2;
    int periodicity;

    double *val;
    double *wei;
    double **wei_shuffle;

    double **x;
    double **x2;

    float ***wei_perm;
    int nperm;

    double prev_d2;
    double mindist;

} subDivide_parms;


int alloc_parms(subDivide_parms *parms, int npoints, int ndims, int *nbins, int nperm);
int free_parms(subDivide_parms *parms);
int alloc_copy_parms(subDivide_parms *parms, subDivide_parms *parmsi);


#endif
