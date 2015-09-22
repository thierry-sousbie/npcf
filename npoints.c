#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "tools.h"
#include "npoints.h"

#define TREE(n) (parms->tree[ n ])
#define NODE(n) (parms->tree[ n ]->node[index[ n ]])

void AutoInterSubDivideN(subDivide_parms *parms,int *index, int ncall, int *binindex);

void SubDivideN(subDivide_parms *parms, int *index)
{
    int bin_index[parms->npairs];
    int i,j;
    double factor;

    for (i=0;i<parms->npairs;i++) bin_index[i]=0;
	
    memset(parms->val,0,sizeof(double)*parms->nbins_tot);
    memset(parms->wei,0,sizeof(double)*parms->nbins_tot);
    for (i=0;i<parms->nperm;i++)
	memset(parms->wei_shuffle[i],0,sizeof(double)*parms->nbins_tot);


    AutoInterSubDivideN(parms,index,0,bin_index);
    
    for (i=0,j=2,factor=1;i<parms->npoints-1;i++)
	if (parms->tree[i]==parms->tree[i+1]) 
	    factor*=(j++);
	    
    for (i=0;i<parms->nbins_tot;i++) parms->val[i]*=factor;
    for (i=0;i<parms->nbins_tot;i++) parms->wei[i]*=factor;

    if (parms->nperm)
	for (j=0;j<parms->nperm;j++)
	    for (i=0;i<parms->nbins_tot;i++)
		parms->wei_shuffle[j][i]*=factor;
}

void AutoInterSubDivideN(subDivide_parms *parms,int *index, int ncall, int *binindex)
{
    int i,j;
    double d2;
    int old_index;
    //printf ("%d %d %d\n",index[0],index[1],index[2]);
    /*
    for (i=0;i<parms->npoints;i++) printf("%d ",index[i]);
    printf(" (%d calls)\n",ncall);
    */
    for (i=0,j=0,d2=1;i<parms->npoints;i++)
    {
	if ((i!=parms->npoints-1)&&
	    (TREE(i)==TREE(i+1))&&
	    (index[i]==index[i+1]))
	    j++;
	else if (j>0)
	{
	    int k,l;
	    int newindex[parms->npoints]; 
	    
	    for (k=0,d2=0;k<=j+1;k++)
	    {
		if ((k>1)&&
		    (NODE(i).next1<TREE(i)->npart)) 
		    continue;
		
		if ((j>k)&&
		    (NODE(i).next2<TREE(i)->npart)) 
		    continue;

		memcpy(newindex,index,parms->npoints*sizeof(int));
		for (l=i-j;l<i-j+k;l++)
		    newindex[l]=NODE(i).next1;
		for (l=i-j+k;l<i+1;l++)
		    newindex[l]=NODE(i).next2;

		
		AutoInterSubDivideN(parms,newindex,ncall,binindex);
	    }
	    j=0;
	}
    }

    if (d2==0) return;

    for (i=0,d2=0;i<TREE(0)->ndims;i++)
	d2 += ((NODE(0).com)[i]-(NODE(1).com)[i])*
	    ((NODE(0).com)[i]-(NODE(1).com)[i]);
    
    if ((parms->periodicity)&&(d2>parms->max_size2))
    {
	double e;
	for (i=0,d2=0;i<TREE(0)->ndims;i++)
	{
	    e=fabs((NODE(0).com)[i]-(NODE(1).com)[i]);
	    if ((parms->periodicity&(1<<i))&&(e>TREE(0)->delta[i]*0.5))
		e=TREE(0)->delta[i]-e;
	    d2+=e*e;
	}
    }

    if (d2<parms->x2[0])
    {
	if (sqrt(d2)+NODE(0).size+NODE(1).size
	    <parms->x[0]+parms->sizetol*(parms->x[1]-parms->x[0]))
	    return ;
    }
    
    if (d2>parms->x2[parms->nbins])
    {
	if (sqrt(d2)-NODE(0).size-NODE(1).size
	    >parms->x[parms->nbins]-parms->sizetol*(parms->x[parms->nbins]-parms->x[parms->nbins-1]))
	    return ;
    }
    
    old_index=binindex[0];
    i=binindex[0];

    if (d2<parms->x2[1]) i=0;
    else if (d2<parms->x2[i])
	do {i--;} while ((i>0)&&(d2<parms->x2[i]));

    if (d2>parms->x2[parms->nbins-1]) i=parms->nbins-1;
    else if (d2>parms->x2[i+1])
	do {i++;} while ((i<parms->nbins)&&(d2>parms->x2[i+1]));
    
   
    binindex[0]=i;
    
    if ((index[0]<TREE(0)->npart)||(NODE(0).size < parms->sizetol*(parms->x[i+1]-parms->x[i])))
    {
	if (ncall==parms->npairs-1)
	{
	    double w1;
	    int cur_index;
	    
	    if ((d2<parms->x2[0])||(d2>parms->x2[parms->nbins])) return;
	    
	   
	    cur_index=coord2index(binindex,parms->nbins,parms->npairs);
	    
	    
	    for (w1=1,i=0;i<parms->npoints;i++) w1*=(double)NODE(i).npart;
	    parms->val[cur_index] += w1;
	    for (w1=1,i=0;i<parms->npoints;i++) w1*=(double)NODE(i).weight;
	    parms->wei[cur_index] += w1;
	    
	    for (i=0;i<parms->npoints;i++) 
		if (TREE(i)->weight) break;

	    if (i<parms->npoints)
	    {
#ifndef TREE_SHUFFLE	
		int k;
		double w;
		
		
		for (i=0;i<parms->nperm;i++)
		{
		    w1=1;
		    for (k=0;k<parms->npoints;k++)
		    {
			if (TREE(k)->weight!=NULL)
			    for (j=NODE(k).start_index,w=0;j<NODE(k).start_index+NODE(k).npart;j++)
				w+=parms->wei_perm[k][i][j];
			else w=NODE(k).npart;
			w1*=w;
		    }
		    parms->wei_shuffle[i][cur_index]+=w1;
		}
#else
		for (i=0;i<parms->nperm;i++)
		{
		    for (w1=1,j=0;j<parms->npoints;j++) 
		    {
			if (TREE(j)->weight!=NULL)
			    w1*=(double)NODE(j).wei_shuffle[i];
			else
			    w1*=(double)NODE(j).npart;
		    }
		    parms->wei_shuffle[i][cur_index]+=w1;
		}
#endif

	    }
	
	}
	else
	{
	    rotateArray(index,sizeof(*index),parms->npoints,1);
	    rotateArray(binindex,sizeof(*binindex),parms->npoints,1);
	    rotateArray(parms->tree,sizeof(*parms->tree),parms->npoints,1);
	    AutoInterSubDivideN(parms,index,ncall+1,binindex);
	    rotateArray(binindex,sizeof(*binindex),parms->npoints,-1);
	    rotateArray(index,sizeof(*index),parms->npoints,-1);
	    rotateArray(parms->tree,sizeof(*parms->tree),parms->npoints,-1);
	}
	   
    }
    else
    {
	int newindex[parms->npoints]; 

	memcpy(newindex,index,parms->npoints*sizeof(int));
	newindex[0]=NODE(0).next1;
	AutoInterSubDivideN(parms,newindex,0,binindex);
	newindex[0]=NODE(0).next2;
	AutoInterSubDivideN(parms,newindex,0,binindex);
    }

    binindex[0]=old_index;
    
    return ;
}

#undef TREE
#undef NODE
