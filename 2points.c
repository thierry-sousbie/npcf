#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "tools.h"
#include "2points.h"

#define TREE(n) (parms->tree[ n ])

/*
void AutoSubDivide2(subDivide_parms *parms, int node_index1, int node_index2, int ncall, int bin_index);
int InterSubDivide2(subDivide_parms *parms, int node_index1, int node_index2, int ncall, int bin_index);
*/
void AutoInterSubDivide2(subDivide_parms *parms, BBTreeNode * node1,BBTreeNode * node2, int ncall);

void SubDivide2(subDivide_parms *parms, BBTreeNode **node)
{
    int i,j;

    memset(parms->val,0,sizeof(double)*parms->nbins_tot);
    memset(parms->wei,0,sizeof(double)*parms->nbins_tot);
    for (i=0;i<parms->nperm;i++)
	memset(parms->wei_shuffle[i],0,sizeof(double)*parms->nbins_tot);


    AutoInterSubDivide2(parms, node[0],node[1],0);

    if (parms->tree[0]==parms->tree[1])
    {
	for (i=0;i<parms->nbins_tot;i++) parms->val[i]*=2;
	for (i=0;i<parms->nbins_tot;i++) parms->wei[i]*=2;

	if (parms->nperm)
	    for (j=0;j<parms->nperm;j++)
		for (i=0;i<parms->nbins_tot;i++)
		    parms->wei_shuffle[j][i]*=2;
    }

}


void AutoInterSubDivide2(subDivide_parms *parms, BBTreeNode * node1,BBTreeNode * node2, int ncall)
{
    int i;
    double d2;
        
    if (node1==node2)
    {

	if (node1->next2->npart>1) 
	    AutoInterSubDivide2(parms, node1->next2,node1->next2,0);

	if (node1->next1->npart>1)
	    AutoInterSubDivide2(parms, node1->next1,node1->next1,0);
	
	AutoInterSubDivide2(parms, node1->next1,node1->next2,0);

	return;
    }
    
    for (i=0,d2=0;i<TREE(0)->ndims;i++)
	d2 += ((node1->com)[i]-(node2->com)[i])*
	    ((node1->com)[i]-(node2->com)[i]);
    
   
    // periodicity and out of range checks
    if ((parms->periodicity)&&(d2>parms->max_size2))
    {
	double e;
	for (i=0,d2=0;i<TREE(0)->ndims;i++)
	{
	    e=fabs((node1->com)[i]-(node2->com)[i]);
	    if ((parms->periodicity&(1<<i))&&(e>TREE(0)->delta[i]*0.5))
		e=TREE(0)->delta[i]-e;
	    d2+=e*e;
	}
    }

    if (d2<(*parms->x2)[0])
    {
	if (sqrt(d2)+node1->size+node2->size<(*parms->x)[0])
	    return ;
    }
    
    if (d2>(*parms->x2)[*parms->nbins])
    {
	if (sqrt(d2)-node1->size-node2->size>(*parms->x)[*parms->nbins])
	    return ;
    }
    
    if ((node1->npart<2)||(node1->size*node1->size < d2*parms->sizetol))
    {
	if (ncall==1)
	{
	    
	    if ((d2<(*parms->x2)[0])||(d2>(*parms->x2)[*parms->nbins])) return;

	    i=*parms->nbins>>1;

	    if (d2<(*parms->x2)[1]) i=0;
	    else if (d2<(*parms->x2)[i])
		do {i--;} while ((i>0)&&(d2<(*parms->x2)[i]));
	    
	    if (d2>(*parms->x2)[*parms->nbins-1]) i=*parms->nbins-1;
	    else if (d2>(*parms->x2)[i+1])
		do {i++;} while ((i<*parms->nbins)&&(d2>(*parms->x2)[i+1]));
	    
	    parms->val[i]+= (double)node1->npart*(double)node2->npart;
	    parms->wei[i]+= (double)node1->weight*(double)node2->weight;
	    if ((TREE(0)->weight!=NULL)||(TREE(1)->weight!=NULL))
	    {
		int j;
#ifndef TREE_SHUFFLE	
		double w1=0;
		double w2=0;
		int k;
		
		for (j=0;j<parms->nperm;j++)
		{
		    w1=w2=0;
		    
		    if (TREE(0)->weight!=NULL)
			for (k=node1->start_index;k<node1->start_index+node1->npart;k++)
			    w1+=parms->wei_perm[0][j][k];
		    else w1=(double)node1->npart;
		    
		    if (TREE(1)->weight!=NULL)
			for (k=node2->start_index;k<node2->start_index+node2->npart;k++)
			    w2+=parms->wei_perm[1][j][k];
		    else w2=(double)node2->npart;

		    parms->wei_shuffle[j][i]+=w1*w2;
		}
#else
		for (j=0;j<parms->nperm;j++)
		{
		    double w1,w2;
		    if (node1->wei_shuffle!=NULL)
			w1=(double)node1->wei_shuffle[j];
		    else
			w1=(double)node1->npart;
		    
		    if (node2->wei_shuffle!=NULL)
			w2=(double)node2->wei_shuffle[j];
		    else
			w2=(double)node2->npart;
		    
		    parms->wei_shuffle[j][i]+= w1*w2;
		}
#endif

	    }
	    
	    return;
	}
	else 
	{
	    //parms->tmp=TREE(1);
	    //TREE(1)=TREE(0);
	    //TREE(0)=parms->tmp;
	    AutoInterSubDivide2(parms, node2, node1, ncall+1);
	    //parms->tmp=TREE(1);
	    //TREE(1)=TREE(0);
	    //TREE(0)=parms->tmp;
	}
    }
    else
    {
	AutoInterSubDivide2(parms, node1->next2, node2, 0);
	AutoInterSubDivide2(parms, node1->next1, node2, 0);
    }

    return;

}

/*
void AutoSubDivide2(subDivide_parms *parms, int node_index1, int node_index2, int ncall, int bin_index)
{
    int i;
    double d2;
    BBTree *tree = (*parms->tree);
    
    if (ncall==2) return;

    if (node_index1 == node_index2)
    {

	if (tree->node[node_index1].next2>=tree->npart) 
	    AutoSubDivide2(parms, tree->node[node_index1].next2,tree->node[node_index1].next2,0,bin_index);

	if (tree->node[node_index1].next1>=tree->npart) 
	    AutoSubDivide2(parms, tree->node[node_index1].next1,tree->node[node_index1].next1,0,bin_index);
	
	AutoSubDivide2(parms, tree->node[node_index1].next1,tree->node[node_index1].next2,0,bin_index);

	return;
    }
    
    for (i=0,d2=0;i<tree->ndims;i++)
	d2 += ((tree->node[node_index1].com)[i]-(tree->node[node_index2].com)[i])*
	    ((tree->node[node_index1].com)[i]-(tree->node[node_index2].com)[i]);
    
   
    // periodicity and out of range checks
    if ((parms->periodicity)&&(d2>parms->max_size2))
    {
	double e;
	for (i=0,d2=0;i<tree->ndims;i++)
	{
	    e=fabs((tree->node[node_index1].com)[i]-(tree->node[node_index2].com)[i]);
	    if ((parms->periodicity&(1<<i))&&(e>tree->delta[i]*0.5))
		e=tree->delta[i]-e;
	    d2+=e*e;
	}
    }
    
    if (d2<parms->x2[0])
    {
	if (sqrt(d2)+tree->node[node_index1].size+tree->node[node_index2].size
	    <parms->x[0]+parms->sizetol*(parms->x[1]-parms->x[0]))
	    return ;
    }
    
    if (d2>parms->x2[parms->nbins])
    {
	if (sqrt(d2)-tree->node[node_index1].size-tree->node[node_index2].size
	    >parms->x[parms->nbins]-parms->sizetol*(parms->x[parms->nbins]-parms->x[parms->nbins-1]))
	    return ;
    }
    
    i=bin_index;

    if (d2<parms->x2[1]) i=0;
    else if (d2<parms->x2[i])
	do {i--;} while ((i>0)&&(d2<parms->x2[i]));

    if (d2>parms->x2[parms->nbins-1]) i=parms->nbins-1;
    else if (d2>parms->x2[i+1])
	do {i++;} while ((i<parms->nbins)&&(d2>parms->x2[i+1]));

    
    if ((node_index1<tree->npart)||(tree->node[node_index1].size < parms->sizetol*(parms->x[i+1]-parms->x[i])))
    {
	if (ncall==1)
	{
	    
	    if ((d2<parms->x2[0])||(d2>parms->x2[parms->nbins])) return;
	    
	    parms->val[i]+= (double)tree->node[node_index1].npart*(double)tree->node[node_index2].npart;
	    parms->wei[i]+= (double)tree->node[node_index1].weight*(double)tree->node[node_index2].weight;
	    if (tree->weight!=NULL)
	    {
		int j;
#ifndef TREE_SHUFFLE	
		double w1=0;
		double w2=0;
		int k;
		
		for (j=0;j<parms->nperm;j++)
		{
		    w1=w2=0;
		    
		    for (k=tree->node[node_index1].start_index;k<tree->node[node_index1].start_index+tree->node[node_index1].npart;k++)
			w1+=parms->wei_perm[0][j][k];
		    for (k=tree->node[node_index2].start_index;k<tree->node[node_index2].start_index+tree->node[node_index2].npart;k++)
			w2+=parms->wei_perm[0][j][k];
		    
		    parms->wei_shuffle[j][i]+=w1*w2;
		}
#else
		for (j=0;j<parms->nperm;j++)
		    parms->wei_shuffle[j][i]+=(double)tree->node[node_index1].wei_shuffle[j]*(double)tree->node[node_index2].wei_shuffle[j];
#endif

	    }
	    
	    return;
	}
	else AutoSubDivide2(parms, node_index2, node_index1, ncall+1,i);
    }
    else
    {
	AutoSubDivide2(parms, tree->node[node_index1].next2, node_index2, 0,i);
	AutoSubDivide2(parms, tree->node[node_index1].next1, node_index2, 0,i);
    }

    return;

}

int InterSubDivide2(subDivide_parms *parms, int node_index1, int node_index2, int ncall, int bin_index)
{
    int i;
    int k1,l1;
    int k2,l2;
    double d2;
    int ntot=0;
    BBTree* tree=parms->tree[0];
    BBTree* tree2=parms->tree[1];

    if (ncall==2) return ntot;

    for (i=0,d2=0;i<tree->ndims;i++)
	d2 += ((tree->node[node_index1].com)[i]-(tree2->node[node_index2].com)[i])*
	    ((tree->node[node_index1].com)[i]-(tree2->node[node_index2].com)[i]);
    
    // periodicity and out of range checks
    if ((parms->periodicity)&&(d2>parms->max_size2))
    {
	//if (parms->periodicity)
	{
	    double e;
	    for (i=0,d2=0;i<tree->ndims;i++)
	    {
		e=fabs((tree->node[node_index1].com)[i]-(tree2->node[node_index2].com)[i]);
		if ((parms->periodicity&(1<<i))&&(e>tree->delta[i]*0.5))
		    e=tree->delta[i]-e;
		d2+=e*e;
	    }
	}
	
//	if (d-tree->node[node_index1].size-tree2->node[node_index2].size > parms->x[parms->nbins])
//	    return 1;
    }

//    if (d<parms->x[0])
//    {
//	if (d+tree->node[node_index1].size+tree2->node[node_index2].size < parms->x[0])
//	    return 1;
//    }

    i=bin_index;

    if (d2<parms->x2[1]) i=0;
    else if (d2<parms->x2[i])
	do {i--;} while ((i>0)&&(d2<parms->x2[i]));

    if (d2>parms->x2[parms->nbins-1]) i=parms->nbins-1;
    else if (d2>parms->x2[i+1])
	do {i++;} while ((i<parms->nbins)&&(d2>parms->x2[i+1]));

    if ((node_index1<tree->npart)||(tree->node[node_index1].size < parms->sizetol*(parms->x[i+1]-parms->x[i])))
    {
	if (ncall==1)
	{
	    if ((d2<parms->x2[0])||(d2>parms->x2[parms->nbins])) return 1;

	    parms->val[i]+= (double)tree->node[node_index1].npart*(double)tree2->node[node_index2].npart;
	    parms->wei[i]+= (double)tree->node[node_index1].weight*(double)tree2->node[node_index2].weight;
	
	    if ((tree->weight!=NULL)||(tree2->weight!=NULL))
	    {
		double w1=0;
		double w2=0;
		int j;
#ifndef TREE_SHUFFLE
		int k;

		for (j=0;j<parms->nperm;j++)
		{
		    w1=w2=0;

		    if (tree->weight!=NULL)
			for (k=tree->node[node_index1].start_index;k<tree->node[node_index1].start_index+tree->node[node_index1].npart;k++)
			    w1+=parms->wei_perm[0][j][k];
		    else
			w1 = (double)tree->node[node_index1].weight;

		    if (tree2->weight!=NULL)
			for (k=tree2->node[node_index2].start_index;k<tree2->node[node_index2].start_index+tree2->node[node_index2].npart;k++)
			    w2+=parms->wei_perm[1][j][k];
		    else
			w2 = (double)tree2->node[node_index2].weight;

		    parms->wei_shuffle[j][i]+=w1*w2;
		}
#else
		    for (j=0;j<parms->nperm;j++)
		    {
			w1=w2=0;

			if (tree->weight!=NULL)
			    w1 = (double)tree->node[node_index1].wei_shuffle[j];
			else
			    w1= (double)tree->node[node_index1].npart;

			if (tree2->weight!=NULL)
			    w2 = (double)tree2->node[node_index2].wei_shuffle[j];
			else
			    w2= (double)tree2->node[node_index2].npart;

			parms->wei_shuffle[j][i]+=w1*w2;
		    }
#endif

	    }
	   
	    return 1;
	}
	else 
	{
	    parms->tmp=parms->tree[1];
	    parms->tree[1]=parms->tree[0];
	    parms->tree[0]=parms->tmp;
	    ntot+=InterSubDivide2(parms, node_index2, node_index1, ncall+1,bin_index);
	    parms->tmp=parms->tree[1];
	    parms->tree[1]=parms->tree[0];
	    parms->tree[0]=parms->tmp;
	}
    }
    else
    {


	k1=tree->node[node_index1].next1;
	k2=tree->node[node_index1].next2;
	
	//if (node_index2>=tree2->npart)
	if (tree2->node[node_index2].size > parms->sizetol*(parms->x[i+1]-parms->x[i]))
	{
	    l1=tree2->node[node_index2].next1;
	    l2=tree2->node[node_index2].next2;
	    
	    ntot+=InterSubDivide2(parms, k1,l1,ncall,bin_index);
	    ntot+=InterSubDivide2(parms, k2,l2,ncall,bin_index);
	    ntot+=InterSubDivide2(parms, k1,l2,ncall,bin_index);
	    ntot+=InterSubDivide2(parms, k2,l1,ncall,bin_index);
	}
	else
	{
	    ntot+=InterSubDivide2(parms, k1,node_index2,ncall,bin_index);
	    ntot+=InterSubDivide2(parms, k2,node_index2,ncall,bin_index);
	}
    }
    
    
    return ntot;

}
*/

#undef TREE
