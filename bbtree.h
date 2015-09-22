#ifndef __BBTREE_HEADER_H__
#define __BBTREE_HEADER_H__

#include "types.h"



typedef struct BBTreeNode_str
{
    struct BBTreeNode_str *next1;
    struct BBTreeNode_str *next2;
    //int prev;
    
    int npart;
    int start_index;
    
    float *com;
    FLOAT size;
    FLOAT weight;

#ifdef TREE_SHUFFLE
    FLOAT *wei_shuffle;
#endif
    
} BBTreeNode;


typedef struct 
{
    int ndims;
    double *x0;
    double *delta;
    int periodicity;
    INT npart;
    INT nnodes;

    float *pos;
    float *weight;
    
    float *val;
#ifdef TREE_SHUFFLE
    int nshuffle;
    float **wei_shuffle;
#endif
    BBTreeNode *node;
} BBTree;

void free_BBtree(BBTree **tree);

BBTree *ComputeBBTree(float *pos,int npart,int ndims,float *weight,float *val,
		      double *x0,double *delta,int Opt_periodic, int inplace, int shuffle);

#endif
