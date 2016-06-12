#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mystring.h"
#include "NDfield.h"
#include "bbtree.h"
#include "correl.h"
#include "tools.h"

#ifdef USE_THREADS
#include "correl_threads.h"
int glob_NThreads=4;
#endif


void Usage(char *fname)
{
  int i;
  
  fprintf (stderr,"\nUsage:\n  %s <dist1 filename> [-dist2 <dist2 fname>] [-npoints <np=2>]\n",CutName(fname));
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-decimate <val=1>] [-periodic] [-periodicity <val>] [-theta <T=0.5>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-weight <fname>] [-weight2 <fname>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-LS [<r_factor=10>]] [-shuffle <Nshuffle=100>] [-outdir <dir>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-lmin <lmin=auto>] [-lmax <lmax=auto>] [-amin <amin=0>] [-amax <amax=90>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-nbins <N=20>] [-nangbins <N=20>] [-linearbins]\n");

#ifdef USE_THREADS
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-threadlevel <level=3>]\n");
  fprintf (stderr,"\n");
  fprintf (stderr,"   Threadlevel is NOT the number of cores, but the thread pool size! Tune it for maximum efficiency.\n");
  fprintf (stderr,"\n");
#else
  fprintf (stderr,"\n");
  fprintf (stderr,"   Threads are not enabled. Set option USE_THREADS in Makefile to use them.\n");
  fprintf (stderr,"\n");
#endif

  fprintf (stderr,"\n");
  fprintf (stderr,"           * '-dist2 <fname>' to compute intercorrelation between dist1 and dist2.\n");
  fprintf (stderr,"             If dist2 is omitted, the autocorellation function of dist1 is computed.\n");
  fprintf (stderr,"             If option '-LS' is also set, dist2 is considered as the random sample.\n");
  fprintf (stderr,"\n");
  fprintf (stderr,"           * '-LS' to use Landy-Szalay estimator (set dist2 to force a random sample).\n");
    fprintf (stderr,"\n");
  fprintf (stderr,"           * '-periodicity 1011' will enable periodic boundary conditions\n");
  fprintf (stderr,"           for dimensions N-3, N-1 and N of a ND field.\n");
  fprintf (stderr,"           By default, boundary conditions are non periodic (use -periodic for \n");
  fprintf (stderr,"           fully periodic boundaries).\n");
  fprintf (stderr,"\n");

  exit(0);
}

int main(int argc, char **argv)
{
    char fullname[255];
    char fname[255];
    char outdir[255];
    char str[255];
    int ndims;
    INT npart;
    float *pos;
    double *x0=NULL;
    double *delta=NULL;
    long i,j;
    
    int Opt_fname=0;
    int Opt_periodic=0;
    int Opt_weight=0;
    int Opt_weight2=0;
    int Opt_value=0;
    //int Opt_nonscalarWeight=0;
    int Opt_outdir=0;
    int Opt_LS=0;
    int Opt_dist2=0;
    int Opt_volume=0;
    int Opt_linearbins=0;
    int Opt_decimate=0;
    int Opt_shuffle=100;
    int Opt_npoints=2;
    

    double lmin=-1;
    double lmax=-1;
    double amin=0;
    double amax=90;
    int nbins=20;
    int nangbins=20;
    double theta=0.5;
    char weight_fname[255];
    char weight2_fname[255];
    char value_fname[255];
    char fullname2[255];
    char fname2[255];
    float *weight=NULL;
    float *weight2=NULL;
    float *value=NULL;
    NDfield *field;
    double volume=0;
    int random_factor=10;
    double decimate=1;

    BBTree *tree[20];
    double *h=NULL;
    double *hw=NULL;
    double **hws=NULL;
    double **x=NULL;
    double **hLS=NULL;

    int npairs=0;
    int nbins_all[20];
    int nbins_tot;
    double min[20];
    double max[20];

    FILE *f;

    for (i=1;i<argc;)
    {
      	if (argv[i][0]!='-')
	{
	    Opt_fname=1;
	    strcpy(fullname,argv[i]);
	    strcpy(fname,CutName(fullname));
	    i++;
	}
	else if (!strcmp(argv[i],"-dist2"))
	{
	    Opt_dist2=1;
	    i++;
	    if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-dist2' needs an argument.\n");
	      Usage(argv[0]);
	    }
	    strcpy(fullname2,argv[i]);
	    strcpy(fname2,CutName(fullname));
	    i++;
	}
	else if (!strcmp(argv[i],"-periodicity"))
	{
	  i++;
	  if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-periodicity' needs an argument.\n");
	      Usage(argv[0]);
	    }
	  Opt_periodic=0;
	  for (j=0;j<strlen(argv[i]);j++)
	    {
	      if (argv[i][j]-'0')
		Opt_periodic |= (1<<j);
	    }
	  i++;
	}
	else if (!strcmp(argv[i],"-periodic"))
	{
	    Opt_periodic=0xffffffff;
	    i++;
	}
	else if (!strcmp(argv[i],"-linearbins"))
	{
	    Opt_linearbins=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-LS"))
	{
	    Opt_LS=1;
	    i++;
	    if ((argc!=i)&&(argv[i][0]!='-'))
	    {
		random_factor=atof(argv[i]);
		i++;
	    }
	}
	else if (!strcmp(argv[i],"-weight"))
	{
	    Opt_weight=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-weight'\n");
		Usage(argv[0]);
	    }
	    strcpy(weight_fname,argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-weight2"))
	{
	    Opt_weight2=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-weight'\n");
		Usage(argv[0]);
	    }
	    strcpy(weight2_fname,argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-value"))
	{
	    Opt_value=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-value'\n");
		Usage(argv[0]);
	    }
	    strcpy(value_fname,argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-outdir"))
	{
	    Opt_outdir=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-outdir'\n");
		Usage(argv[0]);
	    }
	    strcpy(outdir,argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-nbins"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-nbins'\n");
		Usage(argv[0]);
	    }
	    nbins=atoi(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-nangbins"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-nangbins'\n");
		Usage(argv[0]);
	    }
	    nangbins=atoi(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-npoints"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-npoints'\n");
		Usage(argv[0]);
	    }
	    Opt_npoints=atoi(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-shuffle"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-shuffle'\n");
		Usage(argv[0]);
	    }
	    Opt_shuffle=atoi(argv[i]);
	    i++;
	}
#ifdef USE_THREADS
	else if (!strcmp(argv[i],"-threadlevel"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-threadlevel'\n");
		Usage(argv[0]);
	    }
	    glob_NThreads=atoi(argv[i]);
	    i++;
	}
#endif
	else if (!strcmp(argv[i],"-decimate"))
	{
	    Opt_decimate=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-decimate'\n");
		Usage(argv[0]);
	    }
	    decimate=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-theta"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-theta'\n");
		Usage(argv[0]);
	    }
	    theta=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-lmin"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-lmin'\n");
		Usage(argv[0]);
	    }
	    lmin=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-lmax"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-lmax'\n");
		Usage(argv[0]);
	    }
	    lmax=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-amin"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-amin'\n");
		Usage(argv[0]);
	    }
	    amin=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-amax"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-amax'\n");
		Usage(argv[0]);
	    }
	    amax=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-volume"))
	{
	    Opt_volume=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-volume'\n");
		Usage(argv[0]);
	    }
	    volume=atof(argv[i]);
	    i++;
	}
	else 
	{
	    printf ("\nWhat is %s ???\n",argv[i]);
	    Usage(argv[0]);
	}
	
    }
    
    if (!Opt_fname) 
    {
	printf ("\nI need a file name !!!\n");
	Usage(argv[0]);
    }
/*
    if (Opt_dist2&&Opt_LS)
    {
	printf ("\nOptions '-dist2' and '-LS' are not compatible.\n");
	Usage(argv[0]);
    }
*/
    if (!Opt_outdir)
    {
	strcpy(outdir,"./");
    }
    else
    {
	if (outdir[strlen(outdir)-1]!='/')
	    sprintf(outdir,"%s/",outdir);
    }
    
    loadField(fullname,&pos,&npart,&ndims,&x0,&delta,decimate);
    
    if (Opt_weight)
    {
	field=Load_NDfield(weight_fname);

	if (field->ndims!=1)
	    fprintf (stderr,"ERROR: weight must be scalar, not %d dimensional.\n",field->ndims);
	if (field->dims[0]!=npart)
	    fprintf (stderr,"ERROR: I need as many weights val as particles (%d part. != %d wei.)\n",npart,field->dims[0]);

	Convert_NDfield(field,ND_FLOAT);
	weight = (float*) field->val;
	free(field);
    }

    if (Opt_decimate)
	sprintf(fname,"%s.d%.2f",fname,decimate);

    tree[0] = ComputeBBTree(pos,npart,ndims,weight,value,x0,delta,Opt_periodic,1,Opt_shuffle);
    
    if (Opt_dist2)
    {
	loadField(fullname2,&pos,&npart,&ndims,&x0,&delta,decimate);
	if (Opt_weight2)
	{
	    field=Load_NDfield(weight2_fname);
	    
	    if (field->ndims!=1)
		fprintf (stderr,"ERROR: weight2 must be scalar, not %d dimensional.\n",field->ndims);
	    if (field->dims[0]!=npart)
		fprintf (stderr,"ERROR: I need as many weights val as particles (%d part. != %d wei.)\n",npart,field->dims[0]);
	    
	    Convert_NDfield(field,ND_FLOAT);
	    weight2 = (float*) field->val;
	    free(field);
	}
    	tree[1] = ComputeBBTree(pos,npart,ndims,weight2,value,x0,delta,Opt_periodic,1,Opt_shuffle);

    }
    else
	tree[1]=tree[0];

    printPeriodicity(Opt_periodic,ndims);
    
        
    if (volume<=0)
    {
	volume=1.;
	for (i=0;i<tree[0]->ndims;i++) volume*=tree[0]->delta[i];
    }
    
    sprintf (str,"%s%s",outdir,fname);

    npairs = countNPairs(Opt_npoints,(*tree)->ndims);
    nbins_all[0]=nbins;
    *min=lmin;
    *max=lmax;
    nbins_tot=1;
    for (i=1;i<npairs;i++) 
    {
	nbins_all[i]=nangbins;
	min[i] = amin;
	max[i] = amax;
    }

    for (i=0;i<npairs;i++) nbins_tot*=nbins_all[i];

    if (Opt_LS)
    {
	if (!Opt_dist2)
	    tree[1] = GenerateRandomRealTree(tree[0],random_factor,NULL);
	hLS = NPointCorLS(Opt_npoints,tree[0], tree[1], theta, min, max, nbins_all,
			  Opt_linearbins, &h,&hw,&hws, &x,Opt_shuffle, Opt_periodic);
	sprintf (str,"%s.LS",str);
    }
    else
    {
	for (i=2;i<Opt_npoints;i++) tree[i]=tree[1];
	NPointCor(Opt_npoints,tree, theta, min, max, nbins_all, Opt_linearbins, 
		  &h,&hw,&hws, &x,Opt_shuffle, Opt_periodic);
	if (!Opt_dist2)
	    sprintf(str,"%s.auto",str);
	else
	    sprintf(str,"%s_%s.inter",str,fname2);
    }
    
    sprintf(str,"%s.B%2.2dT%2.2d.%dPCF",str,nbins,(int)(100*theta),Opt_npoints);
    
    printf("Saving results to file %s\n",str);
    f=fopen(str,"w");
    for (i=-1;i<nbins_tot;i++)
	fprintf(f,"%s\n",OutputResults(Opt_npoints,(*tree)->ndims,nbins_all,x,h,hw,hws,hLS,Opt_shuffle,i,str));
    fclose(f);
/*

    if (Opt_LS)
    {
	if (!Opt_dist2)
	    tree[1] = GenerateRandomRealTree(tree[0],random_factor,NULL);

	hLS = NPointCorLS(Opt_npoints,tree[0], tree[1], theta, lmin, lmax, nbins, Opt_linearbins, &h,&hw,&hws, &x,Opt_shuffle, Opt_periodic);
	//hLS = TwoPointCorLS(tree[0], tree[1], theta, lmin, lmax, nbins, Opt_linearbins, &h,&hw,&hws, &x,Opt_shuffle, Opt_periodic);

	sprintf(str,"%s%s.LS.B%2.2dT%2.2d.%dPC",outdir,fname,nbins,(int)(100*theta),Opt_npoints); 

	printf("Saving results to file %s\n",str);
	f=fopen(str,"w");
	fprintf (f,"# x Xi Xi_w D_Xi_w (DD) (DR) (RR) (DD_w) (DR_w) (RR_w)\n");
	for (i=0;i<nbins;i++)
	    fprintf(f,"%10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",
		    x[i],h[i]-1,hw[i]-1,Variance(hws,Opt_shuffle,i),hLS[0][i],hLS[1][i],hLS[2][i],hLS[3][i],hLS[4][i],hLS[5][i]);
	fclose(f);
    }
    else
    {
	//TwoPointCor(tree[0], tree[1], theta, lmin, lmax, nbins, Opt_linearbins, &h,&hw,&hws, &x,Opt_shuffle, Opt_periodic);
	for (i=2;i<Opt_npoints;i++) tree[i]=tree[1];
	
		
	NPointCor(Opt_npoints,tree, theta, lmin, lmax, nbins, Opt_linearbins, &h,&hw,&hws, &x,Opt_shuffle, Opt_periodic);
	
	if (!Opt_dist2)
	    sprintf(str,"%s%s.auto.%dPC",outdir,fname,Opt_npoints);
	else
	    sprintf(str,"%s%s_%s.B%2.2dT%2.2d.%dPC",outdir,fname,fname2,nbins,(int)(theta*100),Opt_npoints);
	
	printf("Saving results to file %s\n",str);
	f=fopen(str,"w");
	fprintf (f,"# x Xi Xi_w D_Xi_w\n");
	for (i=0;i<nbins;i++)
	    fprintf(f,"%10.10e %10.10e %10.10e %10.10e\n",x[i],h[i]-1.,hw[i]-1.,Variance(hws,Opt_shuffle,i));
	fclose(f);
    }
*/

    return 0;
}
