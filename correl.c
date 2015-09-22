#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "NDfield.h"
#include "types.h"
#include "bbtree.h"
#include "correl.h"
#include "tools.h"
#include "2points.h"
#include "3points.h"
#include "npoints.h"

#ifdef USE_OPENMP
#include "omp.h"
#endif
#ifdef USE_THREADS
#include "correl_threads.h"
#endif

int subDivide(subDivide_parms *parms)
{
    int i;
    BBTreeNode *node[parms->npoints];
  
    for (i=0;i<parms->npoints;i++)
	node[i]=&parms->tree[i]->node[parms->tree[i]->npart];
    
    if (parms->npoints==2)
	SubDivide2(parms,node);
    else if (parms->npoints==3)
	SubDivide3(parms,node);
/*
    else
	SubDivideN(parms,node);
*/  
    return 0;
}

// this computes 1+Cn where Cn is 'npoints'-points correlation function
// h is 1+Cn
// hw is the weighted couterpart of h
// hws is the weighted and shuffled coutnerpart of h (with 'shuffle' different shuffled realizations)
// critangle is not really a critical angle :) If a group of particles has maximal radius R (distance to furthest particle from COM),
// then it is allowed to be counted in a bin of size S if R<=S*critangle. The value of critangle should be <=1, critangle=0.5 being good
// for most purposes, critangle =1 being ok in some case ...
void NPointCor(int npoints, BBTree **tree, double sizetol, double *lmin, double *lmax, int *nbins, int linearbins,
	       double **h,double **hw, double ***hw_s, double ***x, int shuffle, int periodicity)
{
    double size;
    double tmpd;
    int i,j,k;
    double min[20];
    double max[20];
    double cval[20];
    clock_t time;
    struct timeval wc_time1,wc_time2;
    double volume;
    double density;
    int nperm=shuffle;
    int npairs;
    int nbins_tot;
    float **perm[npoints];
    subDivide_parms parms;
    double *norm;
    
    for (i=0;i<npoints-1;i++)
	if (tree[i] != tree[i+1]) break;
   
    if (i==npoints-1)
	printf("Computing %d-points auto-correlation function (Tol=%.2f) ... ",npoints,sizetol);
    else
	printf("Computing %d-points inter-correlation function (Tol=%.2f) ... ",npoints,sizetol);

    for (i=0;i<npoints-1;i++)
	if (tree[i]->ndims != tree[i+1]->ndims) break;
    if (i!=npoints-1)
    {
	fprintf(stderr,"ERROR: in NPointCor, the samples to correlate must have the same dimensionality.\n");
	exit(-1);
    }

    for (j=0;j<npoints-1;j++)
	for (i=0;i<tree[j]->ndims;i++)
	{
	    if ((tree[j]->x0[i]!=tree[j]->x0[i])||
		(tree[j]->delta[i]!=tree[j]->delta[i]))
	    {
		fprintf(stderr,"ERROR: in NPointCor, the samples to correlate must have the same bounding box.\n");
		exit(-1);
	    }
	}
    printf("\n");
    time=clock();
    gettimeofday(&wc_time1, NULL);
    size=tree[0]->delta[0];
    for (i=1;i<tree[0]->ndims;i++)
    {
	if (tree[0]->delta[i]<size) size = tree[0]->delta[i];
    }
    
    size*=0.5;
  
    max[0]=size;
    min[0]=0.01*((double)2*size)/pow(tree[0]->npart,1./tree[0]->ndims);
    
    if (lmin[0]>0) min[0]=lmin[0];
    if (lmax[0]>0) max[0]=lmax[0];
    
    
    nperm=0;
#ifndef TREE_SHUFFLE
    for (i=0;i<npoints;i++)
	if (tree[i]->weight!=NULL)
	{
	    perm[i]=GeneratePermutation(tree[i]->npart,shuffle,tree[i]->weight,0);
	    nperm = shuffle;
	    while ((i<npoints)&&(tree[i+1]==tree[i]))
	    {
		perm[i+1]=perm[i];
		i++;
	    }
	}
	else perm[i]=NULL;

    if (nperm)
    {
	parms.wei_perm[0]=perm1;
	parms.wei_perm[1]=perm2;
    }
#else
    for (i=0;i<npoints;i++)
	if (tree[i]->nshuffle)
	    nperm=tree[i]->nshuffle;
#endif

    alloc_parms(&parms,npoints,tree[0]->ndims,nbins,nperm);
    npairs = parms.npairs;
    nbins_tot=parms.nbins_tot;

    for (i=1;i<npairs;i++)
    {
	min[i]=lmin[i];
	max[i]=lmax[i];
    }

    if (nperm) printf("(%d permutations) ",nperm);fflush(0);

    if (nperm)
    {
	for (i=0;i<npoints;i++)
	    parms.wei_perm[i]=perm[i];
    }
    
    *h=parms.val;
    *hw=parms.wei;
    *hw_s=parms.wei_shuffle;
    *x=parms.x;
    
    for (i=0;i<npoints;i++) 
	parms.tree[i]=tree[i];
    parms.periodicity=periodicity;
    parms.sizetol = sizetol/nbins[0]*log(max[0]/min[0]); // distance 
    cval[0] = parms.sizetol;
    j=0;

    for (i=1;i<npairs;i++)
    {
	tmpd = tan(sizetol*(max[i]-min[i])/180.*PI/nbins[i]); //angle
	cval[i] = tmpd;
	if (tmpd < parms.sizetol) {parms.sizetol = tmpd;j=i;}
    }

    parms.sizetol = parms.sizetol*parms.sizetol;
    parms.max_size =size;
    parms.max_size2 =parms.max_size*parms.max_size;

    printf ("  Angular tolerance : %.2f degrees\n",atan(sqrt(parms.sizetol))/PI*180.);

    if (npoints>2)
	{
	    if (j==0)
		printf ("  Maximal constraint comes from the %d size bins.\n",nbins[j]);
	    else
		printf ("  Maximal constraint comes from the %d angle%d bins.\n",nbins[j],j);

	    printf ("  I suggest one of the following choices:\n");
	    for (i=0;i<npairs;i++)
	    {
		printf ("    nbins=[%d,%d,%d]\n",
			(int)(sizetol*log(max[0]/min[0])/cval[i]),
			(int)(sizetol*(max[1]-min[1])/180.*PI/atan(cval[i])),
			(int)(sizetol*(max[1]-min[1])/180.*PI/atan(cval[i])));
	    }

	}

    if (linearbins)
    {
	printf ("  I am just a computer so I do whatever I'm told.\n");
	printf ("  Still, it is quite stupid to use linear bins with a tree code ...\n");
	
	for (i=0;i<=nbins[0];i++)
	    (*x)[0][i]=min[0]+ ((double)i/nbins[0])*(max[0]-min[0]);
    }
    else
	for (i=0;i<=nbins[0];i++)
	    (*x)[0][i]=min[0]*pow(10,((double)i/nbins[0])*log10(max[0]/min[0]));
    
    
    //for (k=1,j=nbins[0]+1;k<npairs;j+=nbins[k++]+1)
    for (k=1;k<npairs;k++)
	for (i=0;i<=nbins[k];i++)
	{
	    (*x)[k][i]=(min[k]+((double)i/nbins[k])*(max[k]-min[k])) /180.*PI;
	}

    for (i=0;i<=nbins[0];i++) 
	parms.x2[0][i] = parms.x[0][i]*parms.x[0][i];
    
    //for (k=1,j=nbins[0]+1;k<npairs;j+=nbins[k++]+1)
    for (k=1;k<npairs;k++)
	for (i=0;i<=nbins[k];i++)
	{
	    parms.x2[k][i]=cos(parms.x[k][i]);
	}
   

#ifndef USE_THREADS
    subDivide(&parms);
#else
    threaded_subDivide(&parms);
#endif

    
    
    for (j=0,volume=1;j<tree[0]->ndims;j++)
	volume*=tree[0]->delta[j];

    for (i=0,density=1;i<npoints;i++)
	density*=tree[i]->npart/volume;
    
    
    
    norm=malloc(nbins_tot*sizeof(double));
    for (i=0;i<nbins_tot;i++)
	norm[i]=nPointsNorm(npairs,density,volume, tree[0]->ndims,i,*x,nbins);

    collectData(*h,nbins,npairs);
    collectData(*hw,nbins,npairs);
    if (nperm) 
	for (j=0;j<nperm;j++)
	    collectData((*hw_s)[j],nbins,npairs);

    for (i=0;i<nbins_tot;i++) 
	if (norm[i]>0)
	    (*h)[i]/=norm[i];
	else (*h)[i]=0;

    //for (i=0;i<nbins_tot;i++) (*h)[i]-=1;
       
    for (i=0;i<nbins_tot;i++) 
	for (j=0;j<npoints;j++)
	    norm[i]*=(double)tree[j]->node[tree[j]->npart].weight / tree[j]->npart;
              
    for (i=0;i<nbins_tot;i++) 
	if (norm[i]>0)
	    (*hw)[i]/=norm[i];
	else 
	    (*hw)[i]=0;
    //for (i=0;i<nbins_tot;i++) (*hw)[i]-=1;
    
/*
      for (i=0;i<parms.nbins_tot;i++) 
      printf (" %d %d %d %e \n"
      ,(int)(i/parms.nbins/parms.nbins)
      ,(int)(i/parms.nbins)%parms.nbins
      ,i%parms.nbins
      ,(*h)[i]);
      
      exit(0);
*/

    if (nperm)
    {
	for (j=0;j<nperm;j++)
	{
	    for (i=0;i<nbins_tot;i++) 
		if (norm[i]>0)
		    (*hw_s)[j][i]/=norm[i];
		else
		    (*hw_s)[j][i]=0;
	    //for (i=0;i<nbins_tot;i++) (*hw_s)[j][i]-=1;
	}
    }

    //for (k=0,j=0;k<npairs;j+=nbins[k++]+1)
    for (k=0;k<npairs;k++)
	for (i=0;i<nbins[k];i++)
	    (*x)[k][i]=0.5*((*x)[k][i]+(*x)[k][i+1]);

    for (k=1;k<npairs;k++)
	for (i=0;i<nbins[k];i++)
	    (*x)[k][i]=(*x)[k][i]/PI*180.;
    
    for (k=0;k<npairs;k++)
	free(parms.x2[k]);
    free(parms.x2);
/*
    if (npoints==3)
	for (i=0;i<nbins;i++)
	{
	    (*h)[i] = (*h)[i+nbins*i+nbins*nbins*i];
	    (*hw)[i] = (*hw)[i+nbins*i+nbins*nbins*i];
	    if (nperm)
		for (j=0;j<nperm;j++)
		    (*hw_s)[j][i]=(*hw_s)[j][i+nbins*i+nbins*nbins*i];
	}
*/

#ifndef TREE_SHUFFLE
    for (j=0;j<npoints;j++)
    {
	if (tree[j]->weight!=NULL)
	{
	    for (i=0;i<shuffle;i++) 
		free(parms.wei_perm[j][i]);
	    free(parms.wei_perm[j]);
	    while ((j<npoints)&&(tree[j+1]==tree[j])) j++;
	}
    }
#endif

    gettimeofday(&wc_time2, NULL);

#ifdef PARALLEL
    {
	double delta;
	delta = (double)wc_time2.tv_sec + ((double)wc_time2.tv_usec)*1.E-6;
	delta -= (double)wc_time1.tv_sec + ((double)wc_time1.tv_usec)*1.E-6;
	printf("All done. (%.2f s,g=%.2f)\n",delta,((double)(clock()-time)/CLOCKS_PER_SEC)/delta);fflush(0);
    }
#else
    printf("All done. (%.2f s)\n",(double)(clock()-time)/CLOCKS_PER_SEC);fflush(0);
#endif

}

// this computes 1+Cn where Cn is 'npoints'-points correlation function using the generalized LS fromula.
// h is 1+Cn
// hw is the weighted couterpart of h
// hws is the weighted and shuffled coutnerpart of h (with 'shuffle' different shuffled realizations)
// critangle is not really a critical angle :) If a group of particles has maximal radius R (distance to furthest particle from COM),
// then it is allowed to be counted in a bin of size S if R<=S*critangle. The value of critangle should be <=1, critangle=0.5 being good
// for most purposes, critangle =1 being ok in some case ...
double **NPointCorLS(int npoints, BBTree *tree_D,BBTree *tree_R, double critangle, double *lmin, double *lmax, int *nbins, int linearbins,
		       double **h,double **hw, double ***hws, double ***x, int shuffle, int periodicity)
{
    double *h_n[npoints+1];
    double *h_w[npoints+1];
    double **h_ws[npoints+1];
    double factor[npoints+1];
    BBTree *tree[npoints];
    
    double min[20],max[20],size;
    int i,j,k;
    double (**cor);
    char str[npoints+1];
    int npairs = countNPairs(npoints, tree_D->ndims);
    int nbins_tot = 1;

    for (i=0;i<npairs;i++)
	nbins_tot*=nbins[i];

    size=tree_D->delta[0];
    for (i=1;i<tree_D->ndims;i++)
    {
	if (tree_D->delta[i]<size) size = tree_D->delta[i];
    }
    if (lmax[0]<=0) max[0]=size; 
    else max[0]=lmax[0];

    if (lmin[0]<=0) min[0]=0.01*((double)2*size)/pow(tree_D->npart,1./tree_D->ndims); 
    else min[0]=lmin[0];

    for (i=1;i<npairs;i++)
    {
	min[i]=lmin[i];
	max[i]=lmax[i];
    }

    
    for (i=0;i<=npoints;i++)
    {
	for (j=i;j>=0;j--)
	    if ((j==i)||(j==0)) 
		factor[j]=1;
	    else
		factor[j]=factor[j]+factor[j-1];
    }
    for (i=0;i<=npoints;i++) 
	if (i&1) factor[i]=-factor[i];

    //for (j=0;j<=npoints;j++) printf ("%f ",factor[j]);printf("\n");
    for (j=0;j<=npoints;j++)
    {
	strcpy(str,"");
	for (i=0;i<npoints-j;i++)
	{
	    sprintf(str,"%sD",str);
	    tree[i]=tree_D;
	}
	for (i=npoints-j;i<npoints;i++)
	{
	    sprintf(str,"%sR",str);
	    tree[i]=tree_R;
	}

	printf ("Computing %s:\n",str);
	NPointCor(npoints,tree,critangle,min,max,nbins,linearbins,&h_n[j],&h_w[j],&h_ws[j],x,shuffle,periodicity);

	if (j!=npoints)
	{
	    for (i=0;i<npairs;i++)
		free((*x)[i]);
	    free(*x);
	}
    }
    
    (*h)=realloc(*h,sizeof(double)*nbins_tot);
    (*hw)=realloc(*hw,sizeof(double)*nbins_tot);
    memset((*h),0,sizeof(double)*nbins_tot);
    memset((*hw),0,sizeof(double)*nbins_tot);

    for (i=0;i<nbins_tot;i++)
    {

	//printf("bin[%d]=",i);
	for (j=0;j<=npoints;j++)
	{
	    //printf ("+%f*%f ",factor[j],h_n[j][i]);
	    (*h)[i]+=factor[j]*h_n[j][i];
	}
	//printf ("\n");

	if (h_n[npoints][i]!=0) 
	    (*h)[i]/=(h_n[npoints][i]);
	else (*h)[i]=0;
    }

    for (i=0;i<nbins_tot;i++)
    {
	for (j=0;j<=npoints;j++)
	    (*hw)[i]+=factor[j]*h_w[j][i];

	if (h_w[npoints][i]!=0) 
	    (*hw)[i]/=(h_w[npoints][i]);
	else (*hw)[i]=0;
    }
       
    if (h_ws[0]!=NULL)
    {
	for (i=0;i<shuffle;i++)
	{
	    for (j=0;j<nbins_tot;j++)
	    {
		for (k=1;k<npoints;k++)
		    (h_ws)[0][i][j] += factor[k]*(h_ws)[k][i][j];
		(h_ws)[0][i][j] += factor[k]*(h_w)[k][j];

		if (h_w[npoints][j]!=0)
		    (h_ws)[0][i][j]/=((h_w)[npoints][j]);
		else (h_ws)[0][i][j] = 0;
		    
	    }
	}
	*hws = h_ws[0];
	for (k=1;k<npoints;k++)
	{
	    for (i=0;i<shuffle;i++) free(h_ws[k][i]);
	    free(h_ws[k]);
	}
    }

    cor = malloc(sizeof(double*)*(npoints+1)*2);

    for (i=0;i<npoints+1;i++)
	cor[i]=h_n[i];
    for (i=0;i<npoints+1;i++)
	cor[i+npoints+1]=h_w[i];
  
    return cor;

}

/*

// range for x is computed automatically if lmin<=0 or lmax<=0
// the returned correlation function is expressed in units of the density of tree2.
void TwoPointCor(BBTree *tree,BBTree *tree2, double sizetol, double lmin, double lmax, int nbins, int linearbins,
		 double **h,double **hw, double ***hw_s, double **x, int shuffle, int periodicity)
{
    double size;
    int i,j;
    double min,max;
    double factor=SphereVolume(tree->ndims);
    clock_t time;
    struct timeval wc_time1,wc_time2;
    double norm=0;
    double volume;
    int nperm=shuffle;

    float **perm1=NULL;
    float **perm2=NULL;

    subDivide_parms parms;
    

    if (tree==tree2) 
	printf("Computing %d-points auto-correlation function (Tol=%.2f) ... ",2,sizetol);
    else
	printf("Computing %d-points inter-correlation function (Tol=%.2f) ... ",2,sizetol);
    fflush(0);

    if (tree->ndims!=tree2->ndims)
    {
	fprintf(stderr,"ERROR: in TwoPointCor, the samples to correlate must have the same dimensionality.\n");
	exit(-1);
    }
    for (i=0;i<tree->ndims;i++)
    {
	if ((tree->x0[i]!=tree2->x0[i])||
	    (tree->delta[i]!=tree2->delta[i]))
	{
	    fprintf(stderr,"ERROR: in TwoPointCor, the samples to correlate must have the same bounding box.\n");
	    exit(-1);
	}
    }

    time=clock();
    gettimeofday(&wc_time1, NULL);
    size=tree->delta[0];
    for (i=1;i<tree->ndims;i++)
    {
	if (tree->delta[i]<size) size = tree->delta[i];
    }
  
    size*=0.5;
  
    max=size;
    min=0.01*((double)2*size)/pow(tree->npart,1./tree->ndims);
    
    if (lmin>0) min=lmin;
    if (lmax>0) max=lmax;
    
    

#ifndef TREE_SHUFFLE
    if (tree->weight!=NULL)
   	perm1=GeneratePermutation(tree->npart,shuffle,tree->weight,0);
    else
	perm1=NULL;
    
    if ((tree!=tree2)&&(tree2->weight!=NULL))
   	perm2=GeneratePermutation(tree2->npart,shuffle,tree2->weight,0);
    else
	perm2=parms.wei_perm[0];
    
    if ((tree->weight==NULL)&&(tree2->weight!=NULL))
	perm1=parms.wei_perm[2];

    if (perm1!=NULL)
	nperm = shuffle;
    else
	nperm = 0;

#else
    if (tree->nshuffle) nperm=tree->nshuffle;
    else if (tree2->nshuffle) nperm=tree2->nshuffle;
    else nperm=0;

    parms.wei_perm=NULL;
#endif

    alloc_parms(&parms,2,tree->ndims,&nbins,nperm);

    if (nperm) printf("(%d permutations) ",nperm);fflush(0);

    if (nperm)
    {
	parms.wei_perm[0]=perm1;
	parms.wei_perm[1]=perm2;
    }
    
    *h=parms.val;
    *hw=parms.wei;
    *hw_s=parms.wei_shuffle;
    *x=parms.x;
    parms.tree[0]=tree;
    parms.tree[1]=tree2;
    
    parms.periodicity=periodicity;
    parms.sizetol = sizetol;
    parms.max_size =size;
    parms.max_size2 =parms.max_size*parms.max_size;
        
    if (linearbins)
	for (i=0;i<=nbins;i++)
	    (*x)[i]=min+ ((double)i/nbins)*(max-min);
    else
	for (i=0;i<=nbins;i++)
	    (*x)[i]=min*pow(10,((double)i/nbins)*log10(max/min));
    
    for (i=0;i<=nbins;i++) parms.x2[i] = parms.x[i]*parms.x[i];

#ifndef USE_THREADS
    subDivide(&parms);
#else
    threaded_subDivide(&parms);
#endif

    volume=1.;
    for (i=0;i<tree->ndims;i++) volume*=tree->delta[i];
    
    norm=volume/((double)tree->node[tree->npart].npart*(double)tree2->node[tree2->npart].npart);
    
    for (i=0;i<nbins;i++) (*h)[i]*=norm/(factor*(pow((*x)[i+1],tree->ndims)-pow((*x)[i],tree->ndims)));
    for (i=0;i<nbins;i++) (*h)[i]-=1;
    
    norm=volume/((double)tree->node[tree->npart].weight*(double)tree2->node[tree2->npart].weight);
          
    for (i=0;i<nbins;i++) (*hw)[i]*=norm/(factor*(pow((*x)[i+1],tree->ndims)-pow((*x)[i],tree->ndims)));
    for (i=0;i<nbins;i++) (*hw)[i]-=1;

    if (nperm)
    {
	for (j=0;j<nperm;j++)
	{
	    for (i=0;i<nbins;i++) (*hw_s)[j][i]*=norm/(factor*(pow((*x)[i+1],tree->ndims)-pow((*x)[i],tree->ndims)));
	    for (i=0;i<nbins;i++) (*hw_s)[j][i]-=1;
	}
    }

    for (i=0;i<nbins;i++) (*x)[i]=0.5*((*x)[i]+(*x)[i+1]);
  
    free(parms.x2);
#ifndef TREE_SHUFFLE
    if (tree->weight!=NULL)
    {
	for (i=0;i<shuffle;i++) 
	    free(parms.wei_perm[0][i]);
	free(parms.wei_perm[0]);
    }
    if ((tree!=tree2)&&(tree2->weight!=NULL))
    {
	for (i=0;i<shuffle;i++) 
	    free(parms.wei_perm[1][i]);
	free(parms.wei_perm[1]);
    }
#endif

    gettimeofday(&wc_time2, NULL);
#ifdef PARALLEL
    {
	double delta;
	delta = (double)wc_time2.tv_sec + ((double)wc_time2.tv_usec)*1.E-6;
	delta -= (double)wc_time1.tv_sec + ((double)wc_time1.tv_usec)*1.E-6;
	printf("done. (%.2f s,g=%.2f)\n",delta,((double)(clock()-time)/CLOCKS_PER_SEC)/delta);fflush(0);
    }
#else
    printf("done. (%.2f s)\n",(double)(clock()-time)/CLOCKS_PER_SEC);fflush(0);
#endif
}

// Computes (DD-2DR+RR)/RR
double **TwoPointCorLS(BBTree *tree_D,BBTree *tree_R, double critangle, double lmin, double lmax, int nbins, int linearbins,
		       double **h,double **hw, double ***hw_s, double **x, int shuffle, int periodicity)
{
    double *hdd=NULL;
    double *hdr=NULL;
    double *hrr=NULL;

    double *hdd_w=NULL;
    double *hdr_w=NULL;
    double *hrr_w=NULL;

    double **hdd_ws=NULL;
    double **hdr_ws=NULL;
    double **hrr_ws=NULL;

    double min,max,size;
    int i,j;
    double (**cor);

    size=tree_D->delta[0];
    for (i=1;i<tree_D->ndims;i++)
    {
	if (tree_D->delta[i]<size) size = tree_D->delta[i];
    }
    if (lmax<=0) max=size; else max=lmax;
    if (lmin<=0) min=0.01*((double)2*size)/pow(tree_D->npart,1./tree_D->ndims); else min=lmin;
    
    printf ("Computing DD:\n");
    TwoPointCor(tree_D,tree_D,critangle,min,max,nbins,linearbins,&hdd,&hdd_w,&hdd_ws,x,shuffle,periodicity);
    free(*x);
    
    printf ("Computing DR:\n");
    TwoPointCor(tree_D,tree_R,critangle,min,max,nbins,linearbins,&hdr,&hdr_w,&hdr_ws,x,shuffle,periodicity);
    free(*x);
    
    printf ("Computing RR:\n");
    TwoPointCor(tree_R,tree_R,critangle,min,max,nbins,linearbins,&hrr,&hrr_w,&hrr_ws,x,shuffle,periodicity);
    
    (*h)=realloc(*h,sizeof(double)*nbins);
    (*hw)=realloc(*hw,sizeof(double)*nbins);
    
    for (i=0;i<nbins;i++)
    	(*h)[i]=hdd[i]+hrr[i]-2*hdr[i];
    for (i=0;i<nbins;i++)
    	(*hw)[i]=hdd_w[i]+hrr_w[i]-2*hdr_w[i];
  
    if (hdd_ws!=NULL)
    {
	for (i=0;i<shuffle;i++)
	{
	    for (j=0;j<nbins;j++)
		if (hrr_w[j]!=-1.)
		{
		    (hdd_ws)[i][j] = (hdd_ws[i][j]+hrr_w[j]-2*hdr_ws[i][j])/(1.+hrr_w[j]);
		}
		else
		    (hdd_ws)[i][j] = 0;
	}
	*hw_s = hdd_ws;
	for (i=0;i<shuffle;i++) free(hdr_ws[i]);
	free(hdr_ws);
    }

    for (i=0;i<nbins;i++)
    {
    	if (hrr[i]!=-1.) 
	    (*h)[i]/=(1.+hrr[i]);
	else (*h)[i]=0;
    }

    for (i=0;i<nbins;i++)
    {
    	if (hrr_w[i]!=-1) 
	    (*hw)[i]/=(1.+hrr_w[i]);
	else (*hw)[i]=0;
    }

    cor = malloc(sizeof(double*)*6);

    cor[0]=(hdd);
    cor[1]=(hdr);
    cor[2]=(hrr);

    cor[3]=(hdd_w);
    cor[4]=(hdr_w);
    cor[5]=(hrr_w);
  
    return cor;

}

*/
