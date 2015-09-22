#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "tools.h"
#include "3points.h"

#define TREE(n) (parms->tree[ n ])
#define PX2(n) (parms->x2[ n ])
#define PX(n)  (parms->x[ n ])


void AutoInterSubDivide3(subDivide_parms *parms, BBTreeNode *node1,BBTreeNode *node2,BBTreeNode *node3, int ncall);

void SubDivide3(subDivide_parms *parms, BBTreeNode **node)
{
    int i,j;
    int factor;
    
    memset(parms->val,0,sizeof(double)*parms->nbins_tot);
    memset(parms->wei,0,sizeof(double)*parms->nbins_tot);
    for (i=0;i<parms->nperm;i++)
	memset(parms->wei_shuffle[i],0,sizeof(double)*parms->nbins_tot);
    for (i=0;i<parms->npoints;i++) parms->index[i]=i;
    
    AutoInterSubDivide3(parms, node[0],node[1],node[2],0);
    
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


void AutoInterSubDivide3(subDivide_parms *parms, BBTreeNode *node1,BBTreeNode *node2,BBTreeNode *node3, int ncall)
{
    int i;
    double d2;
       
    if (node1==node2)
    {

	if (node2==node3)
	{
	    if (node1->next2->npart>1)
		AutoInterSubDivide3(parms, 
				    node1->next2,
				    node1->next2,
				    node1->next2,
				    0);

	    if (node1->next2->npart>1) 
		AutoInterSubDivide3(parms, 
				    node1->next1,
				    node1->next2,
				    node1->next2,
				    0);

	    if (node1->next1->npart>1) 
		AutoInterSubDivide3(parms, 
				    node1->next1,
				    node1->next1,
				    node1->next2,
				    0);
	    
	    if (node1->next1->npart>1) 
		AutoInterSubDivide3(parms, 
				    node1->next1,
				    node1->next1,
				    node1->next1,
				    0);
	}
	else
	{
	    if (node1->next2->npart>1) 
		AutoInterSubDivide3(parms, 
				    node1->next2,
				    node1->next2,
				    node3,
				    0);
	   
	    AutoInterSubDivide3(parms, 
				node1->next1,
				node1->next2,
				node3,
				0);
	    
	    if (node1->next1->npart>1) 
		AutoInterSubDivide3(parms, 
				    node1->next1,
				    node1->next1,
				    node3,
				    0);
	}
	return ;
    }
    else if (node2==node3)
    {
	
	if (node2->next2->npart>1) 
	    AutoInterSubDivide3(parms, 
				node1,
				node2->next2,
				node2->next2,
				0);
	
	AutoInterSubDivide3(parms, 
			    node1,
			    node2->next1,
			    node2->next2,
			    0);
	
	if (node2->next1->npart>1) 
	    AutoInterSubDivide3(parms, 
				node1,
				node2->next1,
				node2->next1,
				0);
	
	return ;
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
    /*
    {
	double v[3][NDIMS_MAX];
	double u[3][NDIMS_MAX];
	int id[3];
	double size2[3];
	
	size2[0]=size2[1]=size2[2]=0;
	// here check periodicity 
	for (i=0;i<TREE(0)->ndims;i++)
	{
	    v[0][i]=(node2->com)[i]-(node3->com)[i];
	    size2[0]+=v[0][i]*v[0][i];
	    v[1][i]=(node1->com)[i]-(node3->com)[i];
	    size2[1]+=v[1][i]*v[1][i];
	    v[2][i]=(node1->com)[i]-(node2->com)[i];
	    size2[2]+=v[2][i]*v[2][i];
	}
	
	if (size2[0]<size2[1])
	    if (size2[1]<size2[2]) {id[0]=0;id[1]=1;id[2]=2;}
	    else if (size2[2]>size2[0]) {id[0]=0;id[1]=2;id[2]=1;}
	    else {id[0]=2;id[1]=0;id[2]=1;}
	else // (size2[1]<size2[0])
	    if (size2[0]<size2[2]) {id[0]=1;id[1]=0;id[2]=2;}
	    else if (size2[2]<size2[1]) {id[0]=2;id[1]=1;id[2]=0;}
	    else {id[0]=1;id[1]=2;id[2]=0;}
	
	size2[0]=size2[1]=size2[2]=0;
	
	for (i=0;i<TREE(0)->ndims;i++)
	{
	    u[0][i]=v[id[0]][i];
	    size2[0]+=u[0][i]*u[0][i];
	    
	    u[1][i]=v[id[1]][i];
	    size2[1]+=u[1][i]*u[1][i];
	    if ((id[0]!=1)&&(id[1]!=1))
		u[0][i] = -u[0][i];
	    
	    u[2][i]=v[id[2]][i];
	    size2[2]+=u[2][i]*u[2][i];
	    if ((id[2]!=1)&&(id[1]!=1))
		u[2][i]=-u[2][i];
	}
	
	a1=(u[0][0]*u[1][0]+u[0][1]*u[1][1]+u[0][2]*u[1][2])/
	    sqrt(size2[0]*size2[1]); //cos(a1)
	a2=(u[2][0]*u[1][0]+u[2][1]*u[1][1]+u[2][2]*u[1][2])/
	    sqrt(size2[1]*size2[2]); //cos(a2)

    }
    */


/*
    if (*parms->index == 0)
    {
	double size = d2+(node1->size+node2->size)*(node1->size+node2->size);
	//double a = 2*d2*(node1->size+node2->size);

	if (size<PX2(0)[0])
	{
	    if (sqrt(d2)+(node1->size+node2->size)<PX(0)[0])
		return ;
	}
	
	if (size>PX2(0)[*parms->nbins])
	{
	    if (sqrt(d2)-(node1->size+node2->size)>PX(0)[*parms->nbins])
		return ;
	}
    }
    else
    {
	double size=0;
	
	for (i=0;i<TREE(0)->ndims;i++)
	    size+=((node1->com)[i]-(node2->com)[i])*(node1->com)[i]-(node2->com)[i];

	size*=size/(d2*parms->prev_d2); //cos(a)^2

    }
*/

/*
    if (d2<PX2(0)[0])
    {
	if (sqrt(d2)+(node1->size+node2->size)<PX(0)[0])
	    return ;
    }
    
    if (d2>PX2(0)[*parms->nbins])
    {
	if (sqrt(d2)-(node1->size+node2->size)>PX(0)[*parms->nbins])
	    return ;
    }
*/
      
    if ((node1->npart<2)||(node1->size*node1->size < d2*parms->sizetol))
    {
	if (ncall==2)
	{
	    int cur_index;
	    int bin_index[3];
	    double v[3][NDIMS_MAX];
	    double u[3][NDIMS_MAX];
	    int id[3];
	    double size2[3];
	    double a1,a2;
	    
	    
	    size2[0]=size2[1]=size2[2]=0;
	    // here check periodicity 
	    for (i=0;i<TREE(0)->ndims;i++)
	    {
		v[0][i]=(node2->com)[i]-(node3->com)[i];
		size2[0]+=v[0][i]*v[0][i];
		v[1][i]=(node1->com)[i]-(node3->com)[i];
		size2[1]+=v[1][i]*v[1][i];
		v[2][i]=(node1->com)[i]-(node2->com)[i];
		size2[2]+=v[2][i]*v[2][i];
	    }

	    if (size2[0]<size2[1])
		if (size2[1]<size2[2]) {id[0]=0;id[1]=1;id[2]=2;}
		else if (size2[2]>size2[0]) {id[0]=0;id[1]=2;id[2]=1;}
		else {id[0]=2;id[1]=0;id[2]=1;}
	    else // (size2[1]<size2[0])
		if (size2[0]<size2[2]) {id[0]=1;id[1]=0;id[2]=2;}
		else if (size2[2]<size2[1]) {id[0]=2;id[1]=1;id[2]=0;}
		else {id[0]=1;id[1]=2;id[2]=0;}

	    if ((size2[id[0]]<PX2(0)[0])||(size2[id[0]]>PX2(0)[parms->nbins[0]])) return;

	    size2[0]=size2[1]=size2[2]=0;

	    for (i=0;i<TREE(0)->ndims;i++)
	    {
		u[0][i]=v[id[0]][i];
		size2[0]+=u[0][i]*u[0][i];

		u[1][i]=v[id[1]][i];
		size2[1]+=u[1][i]*u[1][i];
		if ((id[0]!=1)&&(id[1]!=1))
		    u[0][i] = -u[0][i];

		u[2][i]=v[id[2]][i];
		size2[2]+=u[2][i]*u[2][i];
		if ((id[2]!=1)&&(id[1]!=1))
		    u[2][i]=-u[2][i];
	    }

	    a1=(u[0][0]*u[1][0]+u[0][1]*u[1][1]+u[0][2]*u[1][2])/
		sqrt(size2[0]*size2[1]); //cos(a1)
	    a2=(u[2][0]*u[1][0]+u[2][1]*u[1][1]+u[2][2]*u[1][2])/
		sqrt(size2[1]*size2[2]); //cos(a2)
	    
	    //if (acos(a1)/PI*180.+acos(a2)/PI*180. > 179.9)
	    //	printf ("d=%g %g %g, a1=%g, a2=%g\n",sqrt(size2[0]),sqrt(size2[1]),sqrt(size2[2]),acos(a1)/PI*180.,acos(a2)/PI*180.);
	    
            //printf ("min=%g max=%g\n",acos(PX2(2)[0])/PI*180,acos(PX2(2)[parms->nbins[1]])/PI*180);
	    if ((a1>PX2(1)[0])||(a1<PX2(1)[parms->nbins[1]])) return;
	    if ((a2>PX2(2)[0])||(a2<PX2(2)[parms->nbins[2]])) return;

	    //printf ("d=%g %g %g, a1=%g, a2=%g\n",sqrt(size2[0]),sqrt(size2[1]),sqrt(size2[2]),acos(a1)/PI*180.,acos(a2)/PI*180.);

	    //printf("OK\n");
	    // we bin the size of the smallest segment
	    d2=size2[0];
	    i=parms->nbins[0]>>1;
	    if (d2<PX2(0)[1]) i=0;
	    else if (d2<PX2(0)[i])
		do {i--;} while ((i>0)&&(d2<PX2(0)[i]));
	    
	    if (d2>PX2(0)[parms->nbins[0]-1]) i=parms->nbins[0]-1;
	    else if (d2>PX2(0)[i+1])
		do {i++;} while ((i<parms->nbins[0])&&(d2>PX2(0)[i+1]));
	    bin_index[0]=i;

	    //for (i=0;i<=parms->nbins[1];i++)
	    //printf ("%g\n",PX2(1)[i]);

	    //exit(0);
	    // first angle
	    d2 = a1;
	    
	    i=parms->nbins[1]>>1;
	    if (d2>PX2(1)[1]) i=0;
	    else if (d2>PX2(1)[i])
		do {i--;} while ((i>0)&&(d2>PX2(1)[i]));
	    
	    if (d2<PX2(1)[parms->nbins[1]-1]) i=parms->nbins[1]-1;
	    else if (d2<PX2(1)[i+1])
		do {i++;} while ((i<parms->nbins[1])&&(d2<PX2(1)[i+1]));
	    bin_index[1]=i;

	    // second angle
	    d2 = a2;

	    i=parms->nbins[2]>>1;
	    if (d2>PX2(1)[2]) i=0;
	    else if (d2>PX2(2)[i])
		do {i--;} while ((i>0)&&(d2>PX2(2)[i]));
	    
	    if (d2<PX2(2)[parms->nbins[2]-1]) i=parms->nbins[2]-1;
	    else if (d2<PX2(2)[i+1])
		do {i++;} while ((i<parms->nbins[2])&&(d2<PX2(2)[i+1]));
	    bin_index[2]=i;
	    
	    cur_index=coord2index(bin_index,parms->nbins,parms->npairs);
	    
	    parms->val[cur_index]+= 
		(double)node1->npart*
		(double)node2->npart*
		(double)node3->npart;
	    parms->wei[cur_index]+= 
		(double)node1->weight*
		(double)node2->weight*
		(double)node3->weight;

	    if ((TREE(0)->weight!=NULL)||(TREE(1)->weight!=NULL)||(TREE(2)->weight!=NULL))
	    {
		int j;
#ifndef TREE_SHUFFLE	
		double w1,w2,w3;
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

		    if (TREE(2)->weight!=NULL)
			for (k=node3->start_index;k<node3->start_index+node3->npart;k++)
			    w2+=parms->wei_perm[2][j][k];
		    else w2=(double)node3->npart;

		    parms->wei_shuffle[j][cur_index]+=w1*w2*w3;
		}
#else
		for (j=0;j<parms->nperm;j++)
		{
		    double w1,w2,w3;
		    if (node1->wei_shuffle!=NULL)
			w1=(double)node1->wei_shuffle[j];
		    else
			w1=(double)node1->npart;
		    
		    if (node2->wei_shuffle!=NULL)
			w2=(double)node2->wei_shuffle[j];
		    else
			w2=(double)node2->npart;

		    if (node3->wei_shuffle!=NULL)
			w3=(double)node3->wei_shuffle[j];
		    else
			w3=(double)node3->npart;
		    
		    parms->wei_shuffle[j][cur_index]+= w1*w2*w3;
		}
#endif
	    }
	    
	}
	else 
	{
	    //parms->tmp=TREE(0);TREE(0)=TREE(1);TREE(1)=TREE(2);TREE(2)=parms->tmp;
	    //parms->tmpi=parms->index[0];parms->index[0]=parms->index[1];parms->index[1]=parms->index[2];parms->index[2]=parms->tmpi;
   
	    AutoInterSubDivide3(parms, node2, node3,node1, ncall+1);

	    //parms->tmp=TREE(0);TREE(0)=TREE(2);TREE(2)=TREE(1);TREE(1)=parms->tmp;
	    //parms->tmpi=parms->index[0];parms->index[0]=parms->index[2];parms->index[2]=parms->index[1];parms->index[1]=parms->tmpi;
	}
    }
    else
    {
	AutoInterSubDivide3(parms, node1->next1, node2, node3, 0);
	AutoInterSubDivide3(parms, node1->next2, node2, node3, 0);
    }

   
    return;
    
}

#undef TREE
