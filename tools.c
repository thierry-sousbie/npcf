#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gadget_io.h"
#include "mystring.h"
#include "NDfield.h"
#include "bbtree.h"
#include "tools.h"

int countNPairs(int npoints, int ndims)
{
    int i=0;
    int j=0;

    for (i=1;i<npoints;i++)
    	j+=(i<ndims)?i:ndims;
    
    return j;
}

double Variance(double **val, int nreal, int index)
{
    int i;
    double result=0;
    double avg=0;

    if (val==NULL) return 0.;
    
    for (i=0;i<nreal;i++)
    	avg += val[i][index];
    avg/=nreal;

    for (i=0;i<nreal;i++)
	result += ((val[i][index]-avg)*(val[i][index]-avg));
    
    result =sqrt(result/nreal);

    return result;
}

void decimateDist(float **pos, INT *npart, int ndims, double fraction)
{
    long i,j;
    INT n=0;
    double v;

    for (i=0,j=0;i<(*npart)*ndims;i+=ndims)
    {
	v=(double)rand()/RAND_MAX;
	if (v<fraction)
	{
	    if (i!=j) memcpy(&(*pos)[j],&(*pos)[i],ndims*sizeof(float));
	    j+=ndims;
	    n++;
	}
    }

    *npart=n;
    *pos=realloc(*pos,n*sizeof(float)*ndims);
} 

int loadField(char *fname,float **pos, INT *npart, int *ndims, double **x0,double **delta, double decimate)
{
    snapshot_data *snap;
    NDfield *field;
    long i;

    if (!IsNDfield(fname))
    {
	snap=calloc(1,sizeof(snapshot_data));
	ReadGadget(fname,snap,FLAG_POS);
	
	//WriteGadget(fname,snap,NULL,0);exit(0);
	*pos=snap->Pos;
	*npart=snap->header.npartTotal[1];
	*ndims=3;
	*x0=realloc(*x0,(*ndims)*sizeof(double));
	*delta=realloc(*delta,(*ndims)*sizeof(double));
	for (i=0;i<(*ndims);i++) (*x0)[i]=0;
	for (i=0;i<(*ndims);i++) (*delta)[i]=snap->header.BoxSize;
	
	free(snap);
    }
    else
    {
	field=Load_NDfield(fname);
	if (field->ndims!=field->dims[0])
	    fprintf (stderr,"ERROR: NDfield must contain a [ndims,npart,1,1,...] ND array with (N-2) last dims of size 1.\n");
	Convert_NDfield(field,ND_FLOAT);
	*ndims=field->dims[0];
	*npart=field->dims[1];
	*pos=(float*)field->val;
	
	*x0=realloc(*x0,(*ndims)*sizeof(double));
	*delta=realloc(*delta,(*ndims)*sizeof(double));
	memcpy(*x0,field->x0,sizeof(double)*(*ndims));
	memcpy(*delta,field->delta,sizeof(double)*(*ndims));
	
	free(field);
    }

    if (decimate<1)
	decimateDist(pos, npart, *ndims, decimate);

    return 0;
}

void printPeriodicity(int Opt_periodic, int ndims)
{
    int i;

    if (Opt_periodic==0) printf ("Using NON-PERIODIC boundary conditions.\n");
    else if (Opt_periodic==0xffffffff) printf ("Using PERIODIC boundary conditions.\n");
    else
    {
	printf ("Using PERIODIC boundary conditions along axis ");
	for (i=0;i<ndims;i++)
	{
	    if (Opt_periodic&(1<<i))
	    {
		if ((i>0)&&((Opt_periodic&(((1<<(i))-1)))!=0))
		{
		    if ((Opt_periodic&(~((1<<(i+1))-1)))==0) printf (" and ");
		    else printf (", ");
		}
		
		printf ("%d",(int)i+1);
	    }
	}
	printf (".\n");
    }
}

double SphereVolume(int ndims)
{
    double val;
    double fac;
    int i;

    if (ndims&1)
    {
	for (i=3,fac=1;i<=ndims;i+=2) fac*=i;
	val = 2.*pow(2.*PI,(ndims-1)/2)/fac;
    }
    else
    {
	for (i=2,fac=1;i<=ndims;i+=2) fac*=i;
	val = pow(2.*PI,ndims/2)/fac;
    }
    
    return val;
}

INT **GeneratePermutationIndex(INT N,int Nperm)
{
    INT **perm;
    INT i,j;
    INT a,b,c;
    double ratio;
   
    if ((Nperm==0)||(N==0)) return NULL;

    printf ("Generating %d permutations of %d indices ... ",Nperm,N);fflush(0);
   
    ratio = (double)(N-1)/ (double)RAND_MAX;

    perm = calloc(Nperm,sizeof(INT*));

    for (i=0;i<Nperm;i++)
    {
	perm[i]=malloc(N*sizeof(INT));
	for (j=0;j<N;j++) perm[i][j]=j;
	for (j=0;j<N;j++)
	{
	    a=(INT)((double)rand() * ratio);
	    b=(INT)((double)rand() * ratio);
	    c=perm[i][a];
	    perm[i][a]=perm[i][b];
	    perm[i][b]=c;
	}
    }

    printf ("done.\n");
    return perm;
}


float **GeneratePermutation(INT N,int Nperm, float *val,int verbose)
{
    float **perm;
    INT i,j;
    INT a,b,c;
    double ratio;
   
    if ((Nperm==0)||(N==0)||(val==NULL)) return NULL;

    if (verbose) printf ("Generating %d permutations of %d indices ... ",Nperm,N);fflush(0);
   
    ratio = (double)(N-1)/ (double)RAND_MAX;

    perm = calloc(Nperm,sizeof(float*));

    for (i=0;i<Nperm;i++)
    {
	perm[i]=malloc(N*sizeof(float));
	for (j=0;j<N;j++) perm[i][j]=val[j];
	for (j=0;j<N;j++)
	{
	    a=(INT)((double)rand() * ratio);
	    b=(INT)((double)rand() * ratio);
	    c=perm[i][a];
	    perm[i][a]=perm[i][b];
	    perm[i][b]=c;
	}
    }

    if (verbose) printf ("done.\n");
    return perm;
}

BBTree *GenerateRandomRealTree(BBTree *tree, int sizeFactor, NDfield *mask)
{
    BBTree *rtree;
    float *pos;
    long i,j;
    
    pos=malloc((long)tree->npart*sizeFactor*tree->ndims*sizeof(float));
    for (i=0;i<(long)tree->npart*sizeFactor*tree->ndims;i+=tree->ndims)
	for (j=0;j<tree->ndims;j++)
	    pos[i+j]=(double)(rand())*tree->delta[j]/(double)RAND_MAX+tree->x0[j];

    rtree = ComputeBBTree(pos,tree->npart*sizeFactor,tree->ndims,NULL,NULL,tree->x0,tree->delta,0,1,0);
    
    return rtree;
}

void *rotateArray(void *array, int datasize, int nval, int delta)
{
    char tmp[datasize*nval];
    int d=delta;
    

    if (d>=nval) d-=nval;
    if (d<0) d+=nval;

    memcpy(tmp,array,datasize*nval);
    memcpy(array,&tmp[datasize*d],(nval-d)*datasize);
    memcpy(array+(nval-d)*datasize,tmp,d*datasize);
    
    return array;

}

int *rotateArrayOld(int *array, int n, int delta)
{
    int i,j;
    int tmp;

    
    for (i=0,j=delta;i<n-1;i++,j++)
    {
	
	if (j>=n) j-=n;
	if (j<0) j+=n;
	tmp = array[i];
	array[i]=array[j];
	array[j]=tmp;
    }
    
    return array;
}

int coord2index(int *coords, int *size, int ndims)
{
    int i;
    int result=0;
    int fact=1;

    for (i=0;i<ndims;i++)
    {
	result += coords[i]*fact;
	fact*=size[i];
    }

    return result;
}

int index2coords(int index, int *coord,int *dims, int ndims)
{
    INT val;INT old_val;
    int i;
    
    val=1;
    for (i=0;i<ndims;i++)
    {     
	old_val=val;
	val*=dims[i];
	coord[i]=(index%val)/old_val;
    }
    
    return;
}


double nPointsNorm(int npairs, double density,double volume, int ndims, int index,double **bin, int *nbins)
{
    int coord[npairs];
    int dims[npairs];
    double norm=1;
    int i;
    
    for (i=0;i<npairs;i++) dims[i]=nbins[i];
    index2coords(index,coord,dims,npairs);

    if (npairs==0) 
	return density*volume;
   
    if (npairs==1)
    	return density*volume*SphereVolume(ndims)*(pow(bin[0][coord[0]+1],ndims)-pow(bin[0][coord[0]],ndims));
    
    if (npairs==3)
    {
	
	double h;
	double dist[npairs];
	double size[npairs];
	int sort[npairs];
	int j,tmp;
	int a,b,c;

	return density;

	for (i=0;i<npairs;i++)
	{
	    sort[i]=i;
	    dist[i]=0.5*(bin[i][coord[i]]+bin[i][coord[i]+1]);
	    size[i]=bin[i][coord[i]+1]-bin[i][coord[i]];
	}

	for (i=1;i<npairs;i++)
	    for (j=0;j<i;j++)
		if (dist[sort[i]]<dist[sort[j]])
		{
		    tmp=sort[i];
		    sort[i]=sort[j];
		    sort[j]=tmp;
		}
	a=sort[0];b=sort[1];c=sort[2];

	h=0.5*(dist[a]+dist[b]+dist[c]);
	
	h=h*(h-dist[a])*(h-dist[b])*(h-dist[c]); //area
	if (h<0) return -1.;
	h=(2.*sqrt(h))/dist[a];// height from d0

	norm =density*volume;
	norm*=SphereVolume(ndims)*(pow(bin[2][coord[a]+1],ndims)-pow(bin[2][coord[a]],ndims)); 
	norm*=2*PI*h * size[b]*size[c]* (dist[b]*dist[c])/(h*dist[a]);
	//printf ("dist[a],h = (%g,%g)\n",dist[a],h);

/*
	{
	    double r,u,v;
	    double dr,du,dv;
	    
	    r=dist[a];
	    u=dist[b]/dist[a];
	    v=(dist[c]-dist[b])/dist[a];

	    dr=size[a];
	    du=size[b]/dist[a]-u*size[a]/(dist[a]*dist[a]);
	    dv=(size[c]-size[b])/dist[a]-(dist[c]-dist[b])/(dist[a]*dist[a])*size[a];

	    if (dist[a]!=dist[b])
	    {
		if (dist[b]!=dist[c])
		    r = density*volume*8*PI*PI*r*r*r*r*r*u*(u+v)*dr*du*dv;
		else
		    r = density*volume*8*PI*PI*r*r*r*r*r*u*(u+v)*dr*du;
	    }
	    else
	    {
		if (dist[b]!=dist[c])
		    r = density*volume*8*PI*PI*r*r*r*r*r*u*(u+v)*dr*dv;
		else
		    r = density*volume*8*PI*PI*r*r*r*r*r*u*(u+v)*dr;
	    }
		    
	    printf ("%g %g %g : %g = %g \n",dist[a],dist[b],dist[c],norm,r);
	    norm = 8*PI*PI*r*r*r*r*r*u*(u+v)*dr*du*dv;
	}
*/

	return norm;
    }
    
    return norm;
}

double *collectData(double *arr, int *size, int ndims)
{
    int i,j,oldj,k,factor;
    int c[ndims+1];
    int sort[ndims];
    int tmp;

    return arr;

    if (ndims<=1) return arr;

    for (i=0;i<=ndims;i++) c[i]=0;
        
    do {
	for (i=0;i<ndims;i++) sort[i]=i;

	for (i=1;i<ndims;i++)
	    for (j=0;j<i;j++)
		if (c[sort[i]]<c[sort[j]])
		{
		    tmp=sort[i];
		    sort[i]=sort[j];
		    sort[j]=tmp;
		}
	for (i=0,j=0,factor=1;i<ndims;i++,factor*=size[i])
	    j+=c[sort[i]]*factor;

	for (i=0,k=0,factor=1;i<ndims;i++,factor*=size[i])
	    k+=c[i]*factor;
	
	//printf ("%d (%d %d %d) -> %d (%d %d %d)\n",k,c[0],c[1],c[2],j,c[sort[0]],c[sort[1]],c[sort[2]]);
		
	if (j!=k) 
	{
	    arr[j]+=arr[k];
	    arr[k]=arr[j];
	}
	
	j=0;
	do {
	    oldj=j;
	    c[j]++;
	    if (c[j]==size[j%ndims])
		c[j++]=0;
	} while ((oldj!=j)&&(j<=ndims));

	
    } while (c[ndims]==0);
    
    return arr;
}

char *OutputResults(int npoints,int ndims, int *nbins,double **bin_val, 
		    double *val, double *val_w, double **val_ws,double **LS,
		    int shuffle,int index,char *out)
{
    int npairs = countNPairs(npoints,ndims);
    int coord[npairs];
    int dims[npairs];
    double norm=1;
    int i,j,k;
    
    sprintf(out,"");

    if (index<0) 
    {
	sprintf(out,"# ");
	sprintf(out,"%s d",out);
	if (npairs>1) 
	    for (i=1;i<npairs;i++)
		sprintf(out,"%s a%d",out,i);

	sprintf(out,"%s C%d",out,npoints);
	if (val_w!=NULL) sprintf(out,"%s C%d_w",out,npoints);
	if (val_ws!=NULL) sprintf(out,"%s Delta_w",out);

	return out;
    }
	
    
    for (i=0;i<npairs;i++) dims[i]=nbins[i];
    index2coords(index,coord,dims,npairs);
       
    for (k=0;k<npairs;k++)
    	sprintf(out,"%s %14.10g",out,bin_val[k][coord[k]]);
    
    //for (i=0;i<npairs;i++) 
//	sprintf(out,"%s %14.10g",out,bin_val[coord[i]]);

    if (LS==NULL) 
	sprintf(out,"%s %14.10g",out,val[index]-1);
    else
	sprintf(out,"%s %14.10g",out,val[index]);

    if (val_w!=NULL) sprintf(out,"%s %14.10g",out,val_w[index]-1);
    if ((val_ws!=NULL)&&(shuffle>1)) 
	sprintf(out,"%s %14.10g",out,Variance(val_ws,shuffle,index));

    if (LS!=NULL) 
	for (i=0;i<npoints+1;i++)
	    sprintf (out,"%s %14.10g",out,LS[i][index]);
    return out;
}
