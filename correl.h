#ifndef _CORREL_HEADER_H__
#define _CORREL_HEADER_H__

#include "NDfield.h"
#include "bbtree.h"
#include "types.h"

void NPointCor(int npoints, BBTree **tree, double sizetol, 
	       double *lmin, double *lmax, int *nbins, int linearbins, 
	       double **h,double **hw, double ***hw_s, double ***x, 
	       int shuffle, int periodicity);

double **NPointCorLS(int npoints,BBTree *tree_D,BBTree *tree_R, double sizetol, 
		     double *lmin, double *lmax, int *nbins, int linearbins, 
		     double **h,double **hw,double ***hws, double ***x, 
		     int shuffle, int periodicity);

/*
void TwoPointCor(BBTree *tree,BBTree *tree2, double sizetol, double lmin, double lmax, int nbins,
		 int linearbins, double **h,double **hw,double ***hws, double **x, int shuffle, int periodicity);

double **TwoPointCorLS(BBTree *tree_D,BBTree *tree_R, double sizetol, double lmin, double lmax, int nbins, 
		       int linearbins, double **h,double **hw,double ***hws, double **x, int shuffle, int periodicity);
*/

#endif
