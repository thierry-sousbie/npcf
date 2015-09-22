#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "tools.h"
#include "correl_parms.h"


int alloc_parms(subDivide_parms *parms, int npoints, int ndims, int *nbins, int nperm)
{
    int i,j;
    int npairs;
    int nbins_tot;
    int nbinsp1_tot;
    
    parms->npoints=npoints;
    parms->npairs=countNPairs(npoints,ndims);

    npairs = parms->npairs;

    nbins_tot=1;
    for (i=0;i<npairs;i++)
	nbins_tot*=nbins[i];

    nbinsp1_tot=1;
    for (i=0;i<npairs;i++)
	nbinsp1_tot*=(nbins[i]+1);

    parms->nbins = calloc(npairs,sizeof(int));
    for (i=0;i<npairs;i++) parms->nbins[i]=nbins[i];

    parms->nbins_tot = nbins_tot;
    parms->nbinsp1_tot = nbinsp1_tot;

    parms->x=calloc(npairs,sizeof(double*));
    parms->x2=calloc(npairs,sizeof(double*));
    for (i=0;i<npairs;i++)
    {
	parms->x[i]=calloc(nbins[i]+1,sizeof(double));
	parms->x2[i]=calloc(nbins[i]+1,sizeof(double));
    }
    

    parms->val=calloc(nbins_tot,sizeof(double));
    parms->wei=calloc(nbins_tot,sizeof(double));
    parms->tree = calloc(parms->npoints,sizeof(BBTree*));
    parms->index = calloc(parms->npoints,sizeof(int));

    
    parms->nperm=nperm;

    if (nperm)
    {
	parms->wei_shuffle = calloc(nperm,sizeof(double*));
	for (i=0;i<nperm;i++)
	    parms->wei_shuffle[i] = calloc(nbins_tot,sizeof(double));
    }
    else parms->wei_shuffle = NULL;
	
    if (nperm)
    {
	parms->wei_perm=calloc(nperm,sizeof(float**));
    }
  

    return 1;
}

int alloc_copy_parms(subDivide_parms *parms, subDivide_parms *parmsi)
{
    int *nbins = parmsi->nbins;
    int i,j;
    int nbins_tot;
    int nbinsp1_tot;
    int npairs = parmsi->npairs;

    nbins_tot = parmsi->nbins_tot;
    nbinsp1_tot = parmsi->nbinsp1_tot;
    memcpy(parms,parmsi,sizeof(subDivide_parms));
    
    parms->tree = calloc(parmsi->npoints,sizeof(BBTree*));
    parms->index = calloc(parms->npoints,sizeof(int));

    parms->val=calloc(nbins_tot,sizeof(double));
    
    parms->wei=calloc(nbins_tot,sizeof(double));

    parms->x=calloc(npairs,sizeof(double*));
    parms->x2=calloc(npairs,sizeof(double*));
    for (i=0;i<npairs;i++)
    {
	parms->x[i]=calloc(nbins[i]+1,sizeof(double));
	parms->x2[i]=calloc(nbins[i]+1,sizeof(double));
    }
            
    parms->nbins = calloc(npairs,sizeof(int));
    for (i=0;i<npairs;i++) parms->nbins[i]=parmsi->nbins[i];

    memcpy(parms->val,parmsi->val,(nbins_tot)*sizeof(double));
    memcpy(parms->wei,parmsi->wei,(nbins_tot)*sizeof(double));

    memcpy(parms->tree,parmsi->tree,parmsi->npoints*sizeof(BBTree*));
    memcpy(parms->index,parmsi->index,parmsi->npoints*sizeof(int));

    for (i=0;i<npairs;i++)
    {
	memcpy(parms->x[i],parmsi->x[i],(nbins[i]+1)*sizeof(double));
	memcpy(parms->x2[i],parmsi->x2[i],(nbins[i]+1)*sizeof(double));
    }

    parms->nperm=parmsi->nperm;
    if (parmsi->nperm)
    {
	parms->wei_shuffle = calloc (parmsi->nperm,sizeof(double*));
	for (i=0;i<parmsi->nperm;i++)
	{
	    parms->wei_shuffle[i]=calloc(nbins_tot,sizeof(double));
	    memcpy(parms->wei_shuffle[i],parmsi->wei_shuffle[i],nbins_tot*sizeof(double));
	}
    }

    if (parmsi->nperm)
    {
	parms->wei_perm=calloc(parmsi->nperm,sizeof(float**));
	memcpy(parms->wei_perm,parmsi->wei_perm,sizeof(float**)*parmsi->nperm);
    }
    
    return 1;
}


int free_parms(subDivide_parms *parms)
{
    int i;
   
    free(parms->val);
    free(parms->wei);
    free(parms->tree);
    free(parms->index);
    
    free(parms->nbins);
    for (i=0;i<parms->npairs;i++)
    {
	free(parms->x[i]);
	free(parms->x2[i]);
    }
    free(parms->x);
    free(parms->x2);

    if (parms->nperm)
    {
	for (i=0;i<parms->nperm;i++) 
	    free(parms->wei_shuffle[i]);
    }
    free(parms->wei_shuffle);
    memset(parms,0,sizeof(subDivide_parms));

    if (parms->nperm)
	free(parms->wei_perm);

    return 1;
}
