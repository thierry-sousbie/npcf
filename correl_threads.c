#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "correl_threads.h"
#include "bbtree.h"
#include "2points.h"
#include "3points.h"
#include "npoints.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#define TREE(n) (parms->tree[ n ])

extern int glob_NThreads; //declared in npcf.c

void *pthread_subDivide(void *ptr)
{
    threaded_subDivide_parms *thparms = (threaded_subDivide_parms *) ptr;
    subDivide_parms *parms = &(thparms->parms);
    BBTreeNode **node = thparms->node;
/*
#ifdef USE_OPENMP
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    printf("\nThread %d/%d starting...\n",tid,nthreads-1);
#endif
*/
  

    if (parms->npoints==2)
	SubDivide2(parms,node);
    else if (parms->npoints==3)
	SubDivide3(parms,node);
/*
    else
	SubDivideN(parms,index);
*/
    thparms->ret = 1;
 
   return NULL;
}

void CreateThreadsParams_recursive(threaded_subDivide_parms **th_parms, threaded_subDivide_parms *model,int *nparms, 
			   BBTreeNode **node, int level, int maxlevel)
{
    subDivide_parms *parms = &model->parms;
    int npoints =parms->npoints;
    BBTreeNode *newnode[npoints]; 
    int i,j,tmp;

    if (level==maxlevel)
    {
	(*th_parms) = realloc ((*th_parms), ((*nparms)+ 1)*sizeof(threaded_subDivide_parms));
	(*th_parms)[(*nparms)].node = calloc(npoints,sizeof(BBTreeNode *));

	alloc_copy_parms(&((*th_parms)[(*nparms)].parms),&(model->parms));
	memcpy((*th_parms)[(*nparms)].node,node,npoints*sizeof(BBTreeNode *));

	(*nparms)++;
	
	return;
    }

    for (i=0,j=0,tmp=1;i<parms->npoints;i++)
    {
	if ((i!=parms->npoints-1)&&
	    (node[i]==node[i+1]))
	    j++;
	else if (j>0)
	{
	    int k,l;
	    
	    for (k=0,tmp=0;k<=j+1;k++)
	    {
		if ((k>1)&&
		    (node[i]->next1->npart<2)) 
		    continue;
		
		if ((j>k)&&
		    (node[i]->next2->npart<2))
		    continue;

		memcpy(newnode,node,parms->npoints*sizeof(BBTreeNode *));
		for (l=i-j;l<i-j+k;l++)
		    newnode[l]=node[i]->next1;
		for (l=i-j+k;l<i+1;l++)
		    newnode[l]=node[i]->next2;

		CreateThreadsParams_recursive(th_parms,model,nparms,newnode,level+1,maxlevel);
	    }
	    j=0;
	}
    }
    if (tmp==0) return ;

    memcpy(newnode,node,parms->npoints*sizeof(BBTreeNode *));

    if (node[0]->next1->npart<2)
    {
	printf ("ERROR: threadlevel too high for this dataset\n");
	printf ("Lower it using option -threadlevel <N>\n");
	exit(0);
    }

    newnode[0]=node[0]->next1;
    CreateThreadsParams_recursive(th_parms,model,nparms,newnode,level+1,maxlevel);
    
    if (node[0]->next2->npart<2)
    {
	printf ("ERROR: threadlevel too high for this dataset\n");
	printf ("Lower it using option -threadlevel <N>\n");
	exit(0);
    }

    newnode[0]=node[0]->next2;
    CreateThreadsParams_recursive(th_parms,model,nparms,newnode,level+1,maxlevel);
    

    return;
}

int CreateThreadsParams(threaded_subDivide_parms **th_parms, threaded_subDivide_parms *model)
{
    subDivide_parms *parms =&model->parms;
    int npoints = parms->npoints;
    BBTreeNode *node[npoints];
    int i,nparms;
    
    nparms=0;
    for (i=0;i<npoints;i++)
	node[i]=&TREE(i)->node[TREE(i)->npart];

    CreateThreadsParams_recursive(th_parms,model,&nparms,node,0,glob_NThreads);
    
    return nparms;
}

    
int threaded_subDivide(subDivide_parms *parms)
{
    int i,j,k;
    int ret;

    pthread_t *thread;
    int *thread_ret;
    threaded_subDivide_parms *thread_parms=NULL;
    threaded_subDivide_parms thread_model;
    int nparms=0;
    
   
    alloc_copy_parms(&(thread_model.parms),parms);
    nparms = CreateThreadsParams (&thread_parms,&thread_model);
			 
    printf ("  Started %d threads ...\n",nparms);fflush(0);
    
    thread=calloc(nparms,sizeof(pthread_t));
    thread_ret=calloc(nparms,sizeof(int));

#ifndef USE_OPENMP  
    for (i=0;i<nparms;i++)
    {
	thread_ret[i] = pthread_create( &thread[i], NULL, 
					pthread_subDivide, 
					(void*) (&thread_parms[i]));
    }

    for (i=0;i<nparms;i++)
    {
	pthread_join(thread[i],NULL);
	
    }
#endif
#ifdef USE_OPENMP
#pragma omp parallel for	
	for (i=0;i<nparms;i++)
	{
	    pthread_subDivide((void*) (&thread_parms[i]));
	}
#endif

    for (i=0;i<nparms;i++)
	for (j=0;j<parms->nbins_tot;j++)
	{
	    parms->val[j]+=thread_parms[i].parms.val[j];
	    parms->wei[j]+=thread_parms[i].parms.wei[j];
	    if (parms->nperm)
	    {
		for (k=0;k<parms->nperm;k++)
		{
		    parms->wei_shuffle[k][j]+=thread_parms[i].parms.wei_shuffle[k][j];
		}
	    }
	}

    for (i=0;i<nparms;i++)
    {
	free_parms(&(thread_parms[i].parms));
	free(thread_parms[i].node);
    }

    free(thread_parms);
    free_parms(&(thread_model.parms));
    free(thread);
    free(thread_ret);

    return ret;
}

#undef TREE
