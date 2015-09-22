#ifndef __CORREL_TOOLS_H__
#define __CORREL_TOOLS_H__

#include "bbtree.h"
#include "types.h"
#include "NDfield.h"

int countNPairs(int npoints, int ndims);

BBTree *GenerateRandomRealTree(BBTree *tree, int sizeFactor, NDfield *mask);

double Variance(double **val, int nreal, int index);

void decimateDist(float **pos, INT *npart, int ndims, double fraction);

void printPeriodicity(int Opt_periodic, int ndims);

int loadField(char *fname,float **pos, INT *npart, int *ndims, double **x0,double **delta, double decimate);

double SphereVolume(int ndims);

float **GeneratePermutation(INT N,int Nperm, float *val, int verbose);

INT **GeneratePermutationIndex(INT N,int Nperm);

void *rotateArray(void *array, int datasize, int nval, int delta);

int coord2index(int *coords, int *size, int ndims);

int index2coords(int index, int *coord,int *dims, int ndims);

double nPointsNorm(int npairs, double density,double volume, int ndims, int index,double **bin, int *nbins);

double *collectData(double *arr, int *size, int ndims);

char *OutputResults(int npoints,int ndims, int *nbins,double **bin_val, 
		    double *val, double *val_w, double **val_ws,double **LS,
		    int shuffle,int index,char *out);

#endif
