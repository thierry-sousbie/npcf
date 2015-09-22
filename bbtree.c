#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "types.h"
#include "bbtree.h"
#include "tools.h"

static int *glob_index;
static float *glob_distance;

void free_BBtree(BBTree **tree)
{
    long i;

    free((*tree)->x0);
    free((*tree)->delta);
    free((*tree)->pos);
    if ((*tree)->weight!=NULL) free((*tree)->weight);
    for (i=(*tree)->npart;i<(*tree)->nnodes;i++)
	if ((*tree)->node[i].com!=NULL) free((*tree)->node[i].com);
    free((*tree)->node);
    free(*tree);

    *tree=NULL;
}

int Compareint_a(const void *ap, const void *bp)
{
    int a=*((int*)ap);
    int b=*((int*)bp);
    float da = glob_distance[a];
    float db = glob_distance[b];
            
    if (da<db) return -1;
    else if (da>db) return 1;
    else return 0;
}

// return the index of the far particle
int SetNodeStat(BBTree *tree, int start_index, int npart, BBTreeNode *node)
{
    long i,j;
    int k;
    float *pos=tree->pos;
    float *weight=tree->weight; 
    float *val=tree->val; 
    double d;
    float *com;
    double dcom[tree->ndims];
    double w=0;
   
    node->start_index=start_index;
    node->npart=npart;
    
    node->size=0;
    node->weight=0;

    if ((node->com==NULL)&&(npart>1)) node->com=calloc(tree->ndims,sizeof(FLOAT));
    for (k=0;k<tree->ndims;k++) dcom[k]=0;

    if (val==NULL)
    {
	if (weight!=NULL)
	{
	    for (i=(long)start_index;i<(long)(npart+start_index);i++)
	    {
		for (k=0;k<tree->ndims;k++)
		    dcom[k] += pos[(long)i*tree->ndims+k];
		w+=weight[i];
	    }
	    if (npart>1)
		for (k=0;k<tree->ndims;k++) node->com[k]=dcom[k]/npart;
	    else
		for (k=0;k<tree->ndims;k++) node->com=&pos[(long)start_index*tree->ndims];

	    node->weight=w;

#ifdef TREE_SHUFFLE
	    if (tree->nshuffle) 
	    {
		node->wei_shuffle=malloc(sizeof(FLOAT)*tree->nshuffle);
		for (j=0;j<tree->nshuffle;j++)
		{
		    for (i=(long)start_index,w=0;i<(long)(npart+start_index);i++)
		    {
			w+=tree->wei_shuffle[j][i];
		    }
		    node->wei_shuffle[j]=w;
		}
	    }
	    else node->wei_shuffle=NULL;
#endif

	}
	else
	{
	    for (i=(long)start_index*tree->ndims;i<(long)(npart+start_index)*tree->ndims;i+=tree->ndims)
	    {
		for (k=0;k<tree->ndims;k++)
		    dcom[k] += pos[i+k];
	    }
	    if (npart>1)
		for (k=0;k<tree->ndims;k++) node->com[k]=dcom[k]/npart;
	    else
		node->com=&pos[(long)start_index*tree->ndims];

	    node->weight=npart;
	}
    }

    com=node->com;
 
 
    for (i=(long)start_index*tree->ndims;i<(long)(npart+start_index)*tree->ndims;i+=tree->ndims)
    {
	for (k=0,d=0;k<tree->ndims;k++)
	    d+=(pos[i+k]-com[k])*(pos[i+k]-com[k]);
	if (d>node->size) {node->size=d;j=i;}
    }
    node->size=sqrt(node->size);
    return j/tree->ndims;

}

int splitnode(BBTree *tree, BBTreeNode *node, int far_index)//,int level)
{
    double *v;
    float *pos;
    int i,j,k;
    double d;
    int n;
    float tmpf;
    
    v=malloc(tree->ndims*sizeof(double));

    if (node->npart==1) return 1;

    pos = &tree->pos[(long)tree->ndims*far_index];
    for (i=0;i<tree->ndims;i++)
	v[i]=(double)pos[i]-(double)node->com[i];
    for (i=0,d=0;i<tree->ndims;i++)
	d+=v[i]*v[i];
    d=sqrt(d);
    for (i=0;i<tree->ndims;i++)
	v[i]/=d;
    

    for (i=0;i<node->npart;i++)
    {
	glob_index[i]=i;
    }
    
    pos=&tree->pos[(long)node->start_index*tree->ndims];
    for (i=0;i<node->npart;i++,pos+=tree->ndims)
    {
	for (j=0,d=0;j<tree->ndims;j++)
	    d+=v[j]*(pos[j]-node->com[j]);
	glob_distance[i]=d;
    }
    free(v);
    pos=&tree->pos[(long)node->start_index*tree->ndims];
    qsort(glob_index,node->npart,sizeof(int),Compareint_a);

    for (i=0;i<node->npart/2;i++)
    {
	float tmpv[tree->ndims];

	j=glob_index[i];
	if (j<i)
	{
	    do
	    {
		j=glob_index[j];
	    } while(j<i);
	}
	if (i!=j)
	{
	    memcpy(tmpv,&pos[i*tree->ndims],tree->ndims*sizeof(float));
	    memcpy(&pos[i*tree->ndims],&pos[j*tree->ndims],tree->ndims*sizeof(float));
	    memcpy(&pos[j*tree->ndims],tmpv,tree->ndims*sizeof(float));
	    if (tree->weight!=NULL)
	    {
		tmpf=tree->weight[(long)node->start_index+i];
		tree->weight[(long)node->start_index+i]=tree->weight[(long)node->start_index+j];
		tree->weight[(long)node->start_index+j]=tmpf;
#ifdef TREE_SHUFFLE
		for (k=0;k<tree->nshuffle;k++)
		{
		    tmpf=tree->wei_shuffle[k][(long)node->start_index+i];
		    tree->wei_shuffle[k][(long)node->start_index+i]=tree->wei_shuffle[k][(long)node->start_index+j];
		    tree->wei_shuffle[k][(long)node->start_index+j]=tmpf;
		}
#endif	
	    }
	    if (tree->val!=NULL)
	    {
		tmpv[0]=tree->val[(long)node->start_index+i];
		tree->val[(long)node->start_index+i]=tree->val[(long)node->start_index+j];
		tree->val[(long)node->start_index+j]=tmpv[0];
	    }
	}
    }
    
    n=(int) node->npart/2;
    if (n==1)
	node->next1=&tree->node[node->start_index];
    else
	node->next1=&tree->node[tree->nnodes++];
   
    if (node->npart - n ==1)
	node->next2=&tree->node[node->start_index+n];
    else
	node->next2=&tree->node[tree->nnodes++];
    /*
    if (level==6)
    {
	static int tmp=0;
	FILE *f;
	char name[255];
	sprintf(name,"left.%d.%d.dat",level,tmp);
	f=fopen(name,"w");
	for (i=0;i<n;i++){
	    fprintf(f,"%.10e %.10e %.10e\n",
		    tree->pos[(node->start_index+i)*tree->ndims+0],
		    tree->pos[(node->start_index+i)*tree->ndims+1],
		    tree->pos[(node->start_index+i)*tree->ndims+2]);
	}
	fclose(f);
	sprintf(name,"right.%d.%d.dat",level,tmp);
	f=fopen(name,"w");
	for (i=0;i<node->npart-n;i++){
	    fprintf(f,"%.10e %.10e %.10e\n",
		    tree->pos[(node->start_index+n+i)*tree->ndims+0],
		    tree->pos[(node->start_index+n+i)*tree->ndims+1],
		    tree->pos[(node->start_index+n+i)*tree->ndims+2]);
	}
	fclose(f);
	tmp++;
    } 
    */
    //node1->prev=node2->prev=node_index;

    i = SetNodeStat(tree, node->start_index, n, node->next1);
    j = SetNodeStat(tree, node->start_index+n, node->npart - n, node->next2);
    
    n =splitnode(tree, node->next1, i);//,level+1);
    n+=splitnode(tree, node->next2, j);//,level+1);
    
    return n;
}

BBTree *ComputeBBTree(float *pos,int npart,int ndims,float *weight, float *val,
		      double *x0,double *delta,int Opt_periodic, int inplace, int shuffle)
{
    BBTree *tree;
    BBTreeNode *node;
    int far_index;
    clock_t time;
    int nsplit=0;
    double mem;
    double memadd;
    
    
    mem = ((2.*npart-1)*sizeof(BBTreeNode) + 
	   ndims*sizeof(FLOAT)*(npart-1))/1024./1024.;
    memadd = ((long)npart*sizeof(int) +
	      (long)npart*sizeof(FLOAT))/1024./1024.;
#ifdef TREE_SHUFFLE
    if (weight!=NULL)
    {
	memadd += ((long)sizeof(float)*shuffle*npart) /1024. /1024.;
	mem += ((long)shuffle*sizeof(FLOAT)*(2.*npart-1)) / 1024. /1024.;
    }
#endif
    tree=malloc(sizeof(BBTree));

    printf ("Computing %dD BBtree for %d particles (%.1f+%.1f Mo)... ",ndims,npart,mem,memadd);fflush(0);
    time=clock();
    
    tree->ndims=ndims;
    tree->x0=malloc(sizeof(double)*tree->ndims);
    memcpy(tree->x0,x0,sizeof(double)*tree->ndims);
    tree->delta=malloc(sizeof(double)*tree->ndims);
    memcpy(tree->delta,delta,sizeof(double)*tree->ndims);

    tree->npart=npart;

    if (inplace)
	tree->pos=pos;
    else
    {
	tree->pos=malloc((size_t)npart*sizeof(float)*tree->ndims);
	memcpy(tree->pos,pos,(size_t)npart*sizeof(float)*tree->ndims);
    }

    if ((weight==NULL)||(inplace)) 
	tree->weight=weight;
    else
    {
	tree->weight=malloc((size_t)npart*sizeof(float));
	memcpy(tree->weight,weight,(size_t)npart*sizeof(float));
    }

    if ((val==NULL) ||(inplace)) 
	tree->val=val;
    else
    {
	tree->val=malloc((size_t)npart*sizeof(float));
	memcpy(tree->val,val,(size_t)npart*sizeof(float));
    }

    tree->nnodes=npart+1;
    tree->node=calloc((long)2*npart-1,sizeof(BBTreeNode));
    
    node = &(tree->node[npart]);
       
    glob_index = malloc((long)tree->npart*sizeof(int));
    glob_distance = malloc((long)tree->npart*sizeof(FLOAT));

    far_index = SetNodeStat(tree,0,tree->npart,node);
    //node->prev=0;

#ifdef TREE_SHUFFLE
    if (tree->weight!=NULL)
    {
	tree->nshuffle=shuffle;
	tree->wei_shuffle = GeneratePermutation(tree->npart,tree->nshuffle,tree->weight,0);
    }
    else tree->nshuffle=0;
#endif

    nsplit = splitnode(tree,node,far_index);//,0);

#ifdef TREE_SHUFFLE
    if (tree->weight!=NULL)
    {
	int i;

	for (i=0;i<tree->nshuffle;i++)
	    free(tree->wei_shuffle[i]);
	free(tree->wei_shuffle);
	tree->wei_shuffle=NULL;
    }
#endif

    free(glob_index);
    free(glob_distance);

    printf ("done. (%.2f s)\n",(double)(clock()-time)/CLOCKS_PER_SEC);
    //printf ("%d levels\n",nsplit);
    return tree;
}
