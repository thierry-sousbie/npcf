#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "endian.h"
#include "mystring.h"
#include "NDfield.h"

int sizeof_NDfield(int type)
{
  switch (type)
    {
    case ND_UCHAR: return sizeof(unsigned char); break;
    case ND_CHAR: return sizeof(char); break;
    case ND_USHORT: return sizeof(unsigned short); break;
    case ND_SHORT: return sizeof(short); break;
    case ND_UINT: return sizeof(unsigned int); break;
    case ND_INT: return sizeof(int); break;
    case ND_ULONG: return sizeof(unsigned long); break;
    case ND_LONG: return sizeof(long); break;
    case ND_FLOAT: return sizeof(float); break;
    case ND_DOUBLE: return sizeof(double); break;
    }

  return 0;
}

int Get_NDtype(int size,int is_integer, int not_signed)
{
  if (size==sizeof(double))
    {
      if (is_integer)
	return ND_LONG;
      else
	return ND_DOUBLE;
    }
  else
    {
      if (is_integer)
	{
	  if (not_signed)
	    {
	      if (size==sizeof(char))
		return ND_UCHAR;
	      else if (size==sizeof(short))
		return ND_USHORT;
	      else if (size==sizeof(int))
		return ND_UINT;
	    }
	  else
	    {
	      if (size==sizeof(char))
		return ND_CHAR;
	      else if (size==sizeof(short))
		return ND_SHORT;
	      else if (size==sizeof(int))
		return ND_INT;
	    }
	}
      else
	return ND_FLOAT;
    }

  return -1;
}

NDfield *Create_NDfield(int *dims,int ndims,int fdims_index,int type,double *x0, double *delta, void *data,char *comment)
{
  NDfield *field;
  int i;

  if (dims==NULL)
    {
      fprintf (stderr,"Error in Create_NDfield: dimensions needed.");
      return NULL;
    }
  if (ndims==0)
    {
      fprintf (stderr,"Error in Create_NDfield: numberof of dimensions must be >0.");
      return NULL;
    }
 
  field=calloc(1,sizeof(NDfield));
  if (fdims_index) memcpy(field->dims,dims,(ndims-fdims_index)*sizeof(int));
  else memcpy(field->dims,dims,ndims*sizeof(int));
  memcpy(field->dims,dims,ndims*sizeof(int));
  field->ndims = ndims;

  field->nval=1;
  if (fdims_index) field->nval=(long)field->dims[0]*(long)field->dims[1];
  else
    for (i=0;i<ndims;i++)
      field->nval *= field->dims[i];

  field->fdims_index = fdims_index;
  // if (field->fdims_index) field->n_dims=2;
  // else field->n_dims=field->ndims;
  /*
  if (fdims_index>=0)
    field->fdims_index = fdims_index;
  else
    field->fdims_index = ndims-1;
  */
  
  field->datatype=type;
  if ((type==ND_LONG)&&(sizeof(long)==4)) field->datatype=ND_INT;
  if ((type==ND_ULONG)&&(sizeof(long)==4)) field->datatype=ND_UINT;

  field->datasize = sizeof_NDfield(field->datatype);

  if (x0!=NULL)
    memcpy(field->x0,x0,sizeof(double)*ndims);
  if (delta!=NULL)
    memcpy(field->delta,delta,sizeof(double)*ndims);

  if (data==NULL)
    field->val = calloc(field->nval,field->datasize);
  else
    field->val=data;

  if (comment!=NULL)
    strcpy(field->comment,comment);

  return field;
}
/*
NDfield *Create_NDfield(int *dims,int ndims,int fdims_index,int type,double *x0, double *delta, void *data, char *comment)
{
  NDfield *field;
  int i;

  if (dims==NULL)
    {
      fprintf (stderr,"Error in Create_NDfield: dimensions needed.");
      return NULL;
    }
  if (ndims==0)
    {
      fprintf (stderr,"Error in Create_NDfield: numberof of dimensions must be >0.");
      return NULL;
    }
 
  field=calloc(1,sizeof(NDfield));
  memcpy(field->dims,dims,ndims*sizeof(int));
  field->ndims = ndims;

  field->nval=1;
  for (i=0;i<ndims;i++)
    field->nval *= field->dims[i];

  if (fdims_index>=0)
    field->fdims_index = fdims_index;
  else
    field->fdims_index = ndims-1;

  
  field->datatype=type;
  if ((type==ND_LONG)&&(sizeof(long)==4)) field->datatype=ND_INT;
  if ((type==ND_ULONG)&&(sizeof(long)==4)) field->datatype=ND_UINT;

  field->datasize = sizeof_NDfield(field->datatype);

  if (x0!=NULL)
    memcpy(field->x0,x0,sizeof(double)*ndims);
  if (delta!=NULL)
    memcpy(field->delta,delta,sizeof(double)*ndims);

  if (data==NULL)
    field->val = calloc(field->nval,field->datasize);
  else
    field->val=data;

  if (comment!=NULL)
    strcpy(field->comment,comment);

  return field;
}
*/
int Init_NDfield(NDfield *field,int *dims,int ndims,int fdims_index,int type, double *x0,double *delta, void *data,char *comment)
{
  int i;

  memset(field,0,sizeof(NDfield));

  if (dims==NULL)
    {
      fprintf (stderr,"Error in Create_NDfield: dimensions needed.");
      return -1;
    }
  if (ndims==0)
    {
      fprintf (stderr,"Error in Create_NDfield: numberof of dimensions must be >0.");
      return -1;
    }
 
  memcpy(field->dims,dims,ndims*sizeof(int));
  field->ndims = ndims;

  field->nval=1;
  for (i=0;i<ndims;i++)
    field->nval *= field->dims[i];

  if (fdims_index>=0)
    field->fdims_index = fdims_index;
  else
    field->fdims_index = ndims-1;

  
  field->datatype=type;
  if ((type==ND_LONG)&&(sizeof(long)==4)) type=ND_INT;
  if ((type==ND_ULONG)&&(sizeof(long)==4)) type=ND_UINT;

  field->datasize = sizeof_NDfield(field->datatype);

  if (x0!=NULL)
    memcpy(field->x0,x0,sizeof(double)*ndims);
  if (delta!=NULL)
    memcpy(field->delta,delta,sizeof(double)*ndims);

  if (data==NULL)
    field->val = calloc(field->nval,field->datasize);
  else
    field->val=data;

  if (comment!=NULL)
    strcpy(field->comment,comment);

  return 0;
}

int Free_NDfield(NDfield **field)
{
  if (*field == NULL) return -1;

  free((*field)->val);
  free(*field);
  *field=NULL;

  return 0;
}

int Save_NDfield(NDfield *field,char *filename)
{
  int i;
  char tag[16];
  char dummy[160];
  FILE *f;

  printf ("Saving %dD field to file %s ...",field->ndims,filename);fflush(0);

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  f=fopen(filename,"w");

  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),16,f);
  fwrite(&i,sizeof(int),1,f);

  i=sizeof(int)*(NDFIELD_MAX_DIMS+3) + sizeof(double)*(2*NDFIELD_MAX_DIMS) + (160+80)*sizeof(char);
  fwrite(&i,sizeof(int),1,f);
  fwrite(field->comment,sizeof(char),80,f);
  fwrite(&field->ndims,sizeof(int),1,f);
  fwrite(field->dims,sizeof(int),field->ndims,f);
  if (field->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(int),NDFIELD_MAX_DIMS-field->ndims,f);
  fwrite(&field->fdims_index,sizeof(int),1,f);
  fwrite(&field->datatype,sizeof(int),1,f);
  fwrite(field->x0,sizeof(double),field->ndims,f);
  if (field->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-field->ndims,f);
  fwrite(field->delta,sizeof(double),field->ndims,f);
  if (field->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-field->ndims,f);
  fwrite(field->dummy,sizeof(char),160,f);
  fwrite(&i,sizeof(int),1,f);
  
  i=field->nval*sizeof_NDfield(field->datatype);
  fwrite(&i,sizeof(int),1,f);
  fwrite(field->val,sizeof_NDfield(field->datatype),field->nval,f);
  fwrite(&i,sizeof(int),1,f);

  fclose(f);
  printf (" done.\n");
  return 0;
}

int LoadDensity_CIC(char *fname,OLD_density_grid *density)
{
    FILE *f;
    unsigned int i;
    char test[30];
    int swap = 0;
    long offs=0;
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }

    fread(&i,sizeof(int),1,f);
    fread(test,sizeof(char)*30,1,f);
    fread(&i,sizeof(int),1,f);

    if (i!=30*sizeof(char)) swap=1;
 
    if (strcmp(test,"Density grid file")) 
	{
	    printf("-->");
	    return -1;
	}
    else printf ("Loading density grid from %s ... ",fname);fflush(0);

    if (swap) printf ("(Swapping)");fflush(0);
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->Nx),sizeof(int),1,f,swap);
    fread_sw(&(density->Ny),sizeof(int),1,f,swap);
    fread_sw(&(density->Nz),sizeof(int),1,f,swap);
    fread_sw(&(density->NNodes),sizeof(int),1,f,swap);

    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);

    if (density->NNodes != (INT)density->Nx*(INT)density->Ny*(INT)density->Nz)
	density->NNodes = (INT)density->Nx*(INT)density->Ny*(INT)density->Nz;

    if (!(density->grid=(FLOAT *)malloc((size_t)density->NNodes*sizeof(FLOAT))))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }
    
    fread(&i,sizeof(int),1,f);
    for (i=0;i<density->Nx;i++)
      {
	fread_sw(&(density->grid[offs]),sizeof(FLOAT),density->NNodes/density->Nx,f,swap);
	offs += density->NNodes/density->Nx;
      }
    fread(&i,sizeof(int),1,f);
    
    i=sizeof(float);
    fread(&i,sizeof(int),1,f);
    fread_sw(&density->redshift,sizeof(float),1,f,swap);
    fread_sw(&density->smoothing,sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);
    
    density->HasGradiant=0;
    printf ("done.\n");
    return 0;
}



inline int C2I(INT *coord, INT *index,int *dims, int ndims,int periodic)
{
  INT val;
  INT dec;
  int i;
  int out=0;

  *index=0;
  dec=1;
  for (i=0;i<ndims;i++)
    {   
      val=coord[i];
      if (periodic&(1<<i))
	{
	  if (val==dims[i]) val=0;
	  if (val>=dims[i]) val = val%(dims[i]);
	  //else if (val<0) val = (100*(dims[i]) + val)%(dims[i]);
	  else if (val<0) val = dims[i]-((-val)%dims[i]);
	}
      else
	{
	    if (val>=dims[i]) {val = dims[i]-1;out=1;}
	    if (val<0) {val=0;out=1;}
	}
      *index += val*dec;
      dec*=dims[i];
    }
  
  return out;
}

int Save_NDfieldPartial(char *filename, NDfield *header, NDfield *field, int periodic)
{
  int i,k;
  long n;
  long j;
  char tag[16];
  char dummy[160];
  FILE *f;
  long pos;
  long long tmp=0;
  INT index,old_index;
  double coord[NDFIELD_MAX_DIMS];
  INT coordi[NDFIELD_MAX_DIMS];
  INT x0[NDFIELD_MAX_DIMS];
  int elSize;
  int newfile=0;

  printf ("Saving partial %dD field to file %s ...",field->ndims,filename);fflush(0);

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if ((f=fopen(filename,"r"))==NULL) 
    newfile=1;
  else fclose(f);

  if (!newfile) 
    f=fopen(filename,"r+");
  else
    f=fopen(filename,"w");

  fseek(f,0,SEEK_SET);

  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),16,f);
  fwrite(&i,sizeof(int),1,f);

  i=sizeof(int)*(NDFIELD_MAX_DIMS+3) + sizeof(double)*(2*NDFIELD_MAX_DIMS) + (160+80)*sizeof(char);
  fwrite(&i,sizeof(int),1,f);
  fwrite(header->comment,sizeof(char),80,f);
  fwrite(&header->ndims,sizeof(int),1,f);
  fwrite(header->dims,sizeof(int),header->ndims,f);
  if (header->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(int),NDFIELD_MAX_DIMS-header->ndims,f);
  fwrite(&header->fdims_index,sizeof(int),1,f);
  fwrite(&header->datatype,sizeof(int),1,f);
  fwrite(header->x0,sizeof(double),header->ndims,f);
  if (header->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-header->ndims,f);
  fwrite(header->delta,sizeof(double),header->ndims,f);
  if (header->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-header->ndims,f);
  fwrite(header->dummy,sizeof(char),160,f);
  fwrite(&i,sizeof(int),1,f);

  elSize = sizeof_NDfield(field->datatype);
  for (i=0,n=1;i<header->ndims;i++) n*=(long)header->dims[i];

  i=header->nval*sizeof_NDfield(header->datatype);
  fwrite(&i,sizeof(int),1,f);
  
  if (newfile)
    {
      printf(" (new file)");fflush(0);
      pos = ftell(f);
      for (j=0;j<n;j++) 
	fwrite(&tmp,elSize,1,f);
      i=header->nval*sizeof_NDfield(header->datatype);
      fwrite(&i,sizeof(int),1,f);
      fseek(f,pos,SEEK_SET);
    }

  //do the partial write HERE ...
  for (i=0;i<header->ndims;i++) coord[i]=field->x0[i];
  for (i=0;i<header->ndims;i++) x0[i]=coordi[i] = header->dims[i]*(coord[i]-header->x0[i])/header->delta[i];
  for (i=0,n=1;i<field->ndims;i++) n*=field->dims[i];

  j=0;
  old_index=0;
  
  for (j=0;j<n*elSize;j+=elSize)
    {
	k=C2I(coordi,&index,header->dims,header->ndims,periodic);
	
	if (index-old_index) fseek(f,elSize*(index-old_index),SEEK_CUR);
	if (!k) fwrite(field->val + j,elSize,1,f);
	
	i=0;coordi[i]++;
	if (*coordi>=x0[i]+field->dims[i])
	    do {coordi[i]=x0[i];coordi[++i]++;} while(coordi[i]>=x0[i]+header->dims[i]);
	old_index=index+1;
    } 
 
  fclose(f);

  printf (" done.\n");

  return 1;
}

NDfield *Load_NDfieldHeader(char *filename)
{
  int i;
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  NDfield *field;

  char comment[80];
  int dims[NDFIELD_MAX_DIMS];
  int ndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double x0[NDFIELD_MAX_DIMS];
  double delta[NDFIELD_MAX_DIMS];
  int isgrafic=0;
 

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      float tmpf[8];
      fclose(f);
      f=fopen(filename,"r");
      fread_sw(&i,sizeof(int),1,f,swap);
      if (i!=44) {
	swap=1-swap;
	fclose(f);f=fopen(filename,"r");
	fread_sw(&i,sizeof(int),1,f,swap);
      }

      if (i==44)
	{
	  fread_sw(dims,sizeof(int),3,f,swap);
	  fread_sw(tmpf,sizeof(float),8,f,swap);
	  fread_sw(&i,sizeof(int),1,f,swap);
	  strcpy(comment,"GRAFIC file");
	  ndims=3;fdims_index=0;
	  datatype=ND_FLOAT;
	  for (i=0;i<3;i++) delta[i] = tmpf[0]*dims[i]*tmpf[7]/100.;
	  for (i=0;i<3;i++) x0[i]=0;
	  isgrafic=1;
	}
      else
	{
	  fclose(f);
	  fprintf (stderr,"File %s has an unknown format.\n",filename);
	  return NULL;
	}
    }
  else
    {
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(comment,sizeof(char),80,f,swap);
      fread_sw(&ndims,sizeof(int),1,f,swap);
      fread_sw(dims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(&fdims_index,sizeof(int),1,f,swap);
      fread_sw(&datatype,sizeof(int),1,f,swap);
      fread_sw(x0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(delta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(dummy,sizeof(char),160,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
    }

  field = Create_NDfield(dims,ndims,fdims_index,datatype,x0,delta,(void*)(&i),comment);
  field->val=NULL;
  fclose(f);

  memcpy(field->dummy,dummy,sizeof(char)*160);

  return field;
}

NDfield *Load_NDfieldChunk(char *filename, double *x0, double *delta, int periodic)
{
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  char comment[80];
  int ldims[NDFIELD_MAX_DIMS];
  int lndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double lx0[NDFIELD_MAX_DIMS];
  double ldelta[NDFIELD_MAX_DIMS];

  NDfield *field;
  int i;
  INT index;
  INT oldindex;
  int imin[NDFIELD_MAX_DIMS];
  int imax[NDFIELD_MAX_DIMS];
  int deltai[NDFIELD_MAX_DIMS];
  double deltax[NDFIELD_MAX_DIMS];
  double xmin[NDFIELD_MAX_DIMS];
  double xmax[NDFIELD_MAX_DIMS];
  INT coord[NDFIELD_MAX_DIMS];
  double tmp;
  INT tmpi;
  INT n;
  int elSize;
  int isgrafic=0;
  
  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }
 
  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      fclose(f);
      if ((field=Load_NDfieldHeader(filename))==NULL) return NULL;
      strcpy(comment,field->comment);
      lndims=field->ndims;
      memcpy(ldims,field->dims,sizeof(int)*NDFIELD_MAX_DIMS);
      memcpy(lx0,field->x0,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(ldelta,field->delta,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(dummy,field->dummy,sizeof(char)*160);
      fdims_index=field->fdims_index;
      datatype=field->datatype;
      f=fopen(filename,"r");
      fseek(f,44+3*sizeof(int),SEEK_SET);  
      isgrafic=1;
    }
  else
    {
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(comment,sizeof(char),80,f,swap);
      fread_sw(&lndims,sizeof(int),1,f,swap);
      fread_sw(ldims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(&fdims_index,sizeof(int),1,f,swap);
      fread_sw(&datatype,sizeof(int),1,f,swap);
      fread_sw(lx0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(ldelta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(dummy,sizeof(char),160,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
    }
  /*
  for (i=0;i<lndims;i++)
    if (ldelta[i] <= delta[i])
      {
	fprintf (stderr,"ERROR in Load_NDfieldChunk, chunk larger than array.\n");
	return NULL;
      }
  */
  for (i=0;i<lndims;i++)
    {
      tmpi = (int)(1.E-8 + (x0[i]-lx0[i])/(ldelta[i])*ldims[i]) - 1;
      //printf("%d: [%e %e %e-> %d]",i,(x0[i]-lx0[i]),(ldelta[i]),((x0[i]-lx0[i])/(ldelta[i])*ldims[i]),tmpi);
      
      imin[i] = tmpi;
      xmin[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);

      tmpi = (int)(1.E-8 + (x0[i]+delta[i]-lx0[i])/(ldelta[i])*ldims[i]) + 1;
      //printf(" [%e %e %e-> %d]\n",(x0[i]+delta[i]-lx0[i]),(ldelta[i]),((x0[i]+delta[i]-lx0[i])/(ldelta[i])*ldims[i]),tmpi);      
      
      imax[i] = tmpi;
      xmax[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);
         
      if (periodic&(1<<i))
	{
	  tmp = (xmin[i]-lx0[i])/ ldelta[i];
	  tmpi = - ((int)tmp);
	  if (tmpi)
	    {
	      imax[i] += tmpi * ldims[i];
	      imin[i] += tmpi * ldims[i];
	      xmax[i] += tmpi * ldelta[i];
	      xmin[i] += tmpi * ldelta[i];
	    }
	}
      else
	{
	  if (imax[i] > ldims[i]) {imax[i]=ldims[i];xmax[i]=lx0[i]+ldelta[i];}
	  if (imin[i] < 0) {imin[i]=0;xmin[i]=lx0[i];}
	}
      
      deltai[i] = imax[i]-imin[i];
      deltax[i] = xmax[i]-xmin[i];
    }


  
  tmpi=1;
  for (i=0;i<lndims;i++) coord[i]=imin[i];
  for (i=0;i<lndims;i++) tmpi*=deltai[i];

  /*
  for (i=0;i<lndims;i++) printf ("[%f,%f]",x0[i],x0[i]+delta[i]);
  printf("\n");
  for (i=0;i<lndims;i++) printf ("[%f,%f]",xmin[i],xmin[i]+deltax[i]);
  printf("\n");
  for (i=0;i<lndims;i++) printf ("[%d,%d]",imin[i],imax[i]);
  printf("\n");
  */

  field = Create_NDfield(deltai,lndims,fdims_index,datatype,xmin,deltax,NULL,comment);
  memcpy(field->dummy,dummy,160*sizeof(char));
  //for (i=0;i<2;i++) printf ("final %f %f\n",xmin[i],deltax[i]);
  //for (i=0;i<2;i++) printf ("cut %f %f\n",x0[i],delta[i]);
  //for (i=0;i<2;i++) printf ("init %f %f\n",lx0[i],ldelta[i]);
  

  oldindex=0;
  elSize = sizeof_NDfield(field->datatype);
  printf ("Loading %dD chunk from file %s ...",lndims,filename);fflush(0);
  for(n=0;n<tmpi*elSize;n+=elSize)
    {
      C2I(coord,&index,ldims,lndims,periodic);
      
      if (index-oldindex) fseek(f,elSize*(index-oldindex),SEEK_CUR);

      if (isgrafic)
	{
	  int a,b,c;
	  c=ldims[1]*ldims[0];
	  a=(int)(index/c);
	  b=(int)(oldindex/c);
	  
	  if (a!=b) 
	    fseek(f,sizeof(int)*2*(a-b),SEEK_CUR);
	  
	}

      fread_sw(field->val + n,elSize,1,f,swap);
      i=0;coord[i]++;
      if (*coord>=*imax)
	do {coord[i]=imin[i];coord[++i]++;} while(coord[i]>=imax[i]);
      oldindex=index+1;
    }
 
  fclose(f);
  printf(" done.\n");
  return field;
}

NDfield *Load_NDfieldChunkHeader(char *filename, double *x0, double *delta, int periodic)
{
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  char comment[80];
  int ldims[NDFIELD_MAX_DIMS];
  int lndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double lx0[NDFIELD_MAX_DIMS];
  double ldelta[NDFIELD_MAX_DIMS];

  NDfield *field;
  int i;
  int imin[NDFIELD_MAX_DIMS];
  int imax[NDFIELD_MAX_DIMS];
  int deltai[NDFIELD_MAX_DIMS];
  double deltax[NDFIELD_MAX_DIMS];
  double xmin[NDFIELD_MAX_DIMS];
  double xmax[NDFIELD_MAX_DIMS];
  double tmp;
  INT tmpi;
  int isgrafic=0;
    
  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }
 
  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      fclose(f);
      if ((field=Load_NDfieldHeader(filename))==NULL) return NULL;
      strcpy(comment,field->comment);
      lndims=field->ndims;
      memcpy(ldims,field->dims,sizeof(int)*NDFIELD_MAX_DIMS);
      memcpy(lx0,field->x0,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(ldelta,field->delta,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(dummy,field->dummy,sizeof(char)*160);
      fdims_index=field->fdims_index;
      datatype=field->datatype;
      f=fopen(filename,"r");
      fseek(f,44+3*sizeof(int),SEEK_SET);  
      isgrafic=1;
    }
  else
    {
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(comment,sizeof(char),80,f,swap);
      fread_sw(&lndims,sizeof(int),1,f,swap);
      fread_sw(ldims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(&fdims_index,sizeof(int),1,f,swap);
      fread_sw(&datatype,sizeof(int),1,f,swap);
      
      fread_sw(lx0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(ldelta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(dummy,sizeof(char),160,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
    }

  for (i=0;i<lndims;i++)
    {
      tmpi = (int)(1.E-8 + (x0[i]-lx0[i])/(ldelta[i])*ldims[i]) - 1;
      imin[i] = tmpi;
      xmin[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);
      
      tmpi = (int)(1.E-8 + (x0[i]+delta[i]-lx0[i])/(ldelta[i])*ldims[i]) + 1;
      imax[i] = tmpi;
      xmax[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);
      
    
      if (periodic&(1<<i))
	{
	  tmp = (xmin[i]-lx0[i])/ ldelta[i];
	  tmpi = - ((int)tmp);
	  if (tmpi)
	    {
	      imax[i] += tmpi * ldims[i];
	      imin[i] += tmpi * ldims[i];
	      xmax[i] += tmpi * ldelta[i];
	      xmin[i] += tmpi * ldelta[i];
	    }
	}
      else
	{
	  if (imax[i] > ldims[i]) {imax[i]=ldims[i];xmax[i]=lx0[i]+ldelta[i];}
	  if (imin[i] < 0) {imin[i]=0;xmin[i]=lx0[i];}
	}
      
      deltai[i] = imax[i]-imin[i];
      deltax[i] = xmax[i]-xmin[i];
    }

  field = Create_NDfield(deltai,lndims,fdims_index,datatype,xmin,deltax,(void *)(&tmpi),comment);
  memcpy(field->dummy,dummy,160*sizeof(char));
  field->val=NULL;
  
  fclose(f);
  
  return field;
}

int IsNDfield_ASCII(char *filename)
{
  //if (fname==std::string("")) return true;
  FILE *f = fopen(filename,"r");
  if (f==NULL) return 0;
  char *line=NULL;
  int n;
  Mygetline(&line,&n,f);      
  fclose(f);
      
  if (strstr(line,NDFIELD_ASCII_TAG)!=NULL) 
    {
      free(line);
      return 1;
    }
  free(line);
  return 0;
}

int IsNDfield(char *filename)
{
  int i;
  char tag[16];
  FILE *f;
  int swap=0;
  
  memset(tag,0,16*sizeof(char));
   
  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return 0;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  if (i!=16) swap=1-swap;
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  
  fclose(f); 
  tag[15]='\0';

  if (!strcmp(tag,NDFIELD_TAG)) return 1;
  if (IsNDfield_ASCII(filename)) return 2;

  return 0;
}

NDfield *Load_NDfield_BIN(char *filename)
{
  int i;
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  NDfield *field;

  char comment[80];
  int dims[NDFIELD_MAX_DIMS];
  int ndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double x0[NDFIELD_MAX_DIMS];
  double delta[NDFIELD_MAX_DIMS];
  OLD_density_grid density;

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;
  
  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      fclose(f);      
      if (LoadDensity_CIC(filename,&density))
	{
	  fprintf (stderr,"File %s has an unknown format.\n",filename);
	  return NULL;
	}
      if (sizeof(FLOAT)==8) datatype=ND_DOUBLE;
      else datatype=ND_FLOAT;
      strcpy(comment,"CIC format");
      ndims=3;dims[0]=density.Nx;dims[1]=density.Ny;dims[2]=density.Nz;
      fdims_index=0;
      x0[0]=density.x0;x0[1]=density.y0;x0[2]=density.z0;
      delta[0]=density.dx*density.Nx;delta[1]=density.dy*density.Ny;delta[2]=density.dz*density.Nz;
      field = Create_NDfield(dims,ndims,fdims_index,datatype,x0,delta,density.grid,comment);
      //printf ("d=%e %e\n",x0[0],delta[1]);
      return field;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(comment,sizeof(char),80,f,swap);
  fread_sw(&ndims,sizeof(int),1,f,swap);
  fread_sw(dims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
  fread_sw(&fdims_index,sizeof(int),1,f,swap);
  fread_sw(&datatype,sizeof(int),1,f,swap);
  fread_sw(x0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
  fread_sw(delta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
  fread_sw(dummy,sizeof(char),160,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  
  printf ("Loading %dD field from file %s ...",ndims,filename);fflush(0);
  if (swap) printf ("(swapping)\n");
  field = Create_NDfield(dims,ndims,fdims_index,datatype,x0,delta,NULL,comment);
  memcpy(field->dummy,dummy,160*sizeof(char));

  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(field->val,sizeof_NDfield(field->datatype),field->nval,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);

  fclose(f);
  printf (" done.\n");
  return field;
}

int Convert_NDfield(NDfield *field,int type)
{
  void *new_tab;
  long i;

  unsigned char *ptruc;
  char *ptrc;
  unsigned short *ptrus;
  short *ptrs;
  unsigned int*ptrui;
  int *ptri;
  unsigned long *ptrul;
  long *ptrl;
  float *ptrf;
  double *ptrd;

  unsigned char *newuc;
  char *newc;
  unsigned short *newus;
  short *news;
  unsigned int*newui;
  int *newi;
  unsigned long *newul;
  long *newl;
  float *newf;
  double *newd;
  

  if (type==field->datatype)
    return 0;

  if (sizeof_NDfield(type) > sizeof_NDfield(field->datatype))
    new_tab = calloc(field->nval,sizeof_NDfield(type));
  else
    new_tab = field->val;

  switch (field->datatype)
    {
    case ND_UCHAR: ptruc = (unsigned char*)field->val;break;
    case ND_CHAR: ptrc = (char*)field->val;break;
    case ND_USHORT: ptrus = (unsigned short*)field->val;break;
    case ND_SHORT: ptrs = (short*)field->val;break;
    case ND_UINT: ptrui = (unsigned int*)field->val;break;
    case ND_INT: ptri = (int*)field->val;break;
    case ND_ULONG:ptrul = (unsigned long*)field->val;break;
    case ND_LONG: ptrl = (long*)field->val;break;
    case ND_FLOAT: ptrf = (float*)field->val;break;
    case ND_DOUBLE: ptrd = (double*)field->val;break;
    }

  switch (type)
    {
    case ND_UCHAR: newuc = (unsigned char*)new_tab;break;
    case ND_CHAR: newc = (char*)new_tab;break;
    case ND_USHORT: newus = (unsigned short*)new_tab;break;
    case ND_SHORT: news = (short*)new_tab;break;
    case ND_UINT: newui = (unsigned int*)new_tab;break;
    case ND_INT: newi = (int*)new_tab;break;
    case ND_ULONG:newul = (unsigned long*)new_tab;break;
    case ND_LONG: newl = (long*)new_tab;break;
    case ND_FLOAT: newf = (float *)new_tab;break;
    case ND_DOUBLE: newd = (double*)new_tab;break;
    }

  switch (field->datatype)
    {
    case ND_UCHAR:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptruc[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptruc[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptruc[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptruc[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptruc[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptruc[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptruc[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptruc[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptruc[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptruc[i]; break; 
	  }
      break;
    case ND_CHAR:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrc[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrc[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrc[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrc[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrc[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrc[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrc[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrc[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrc[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrc[i]; break; 
	  }
      break;
    case ND_USHORT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrus[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrus[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrus[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrus[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrus[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrus[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrus[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrus[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrus[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrus[i]; break; 
	  }
      break;
    case ND_SHORT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrs[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrs[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrs[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrs[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrs[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrs[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrs[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrs[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrs[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrs[i]; break; 
	  }
      break;
    case ND_UINT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrui[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrui[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrui[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrui[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrui[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrui[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrui[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrui[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrui[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrui[i]; break; 
	  }
      break;
    case ND_INT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptri[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptri[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptri[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptri[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptri[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptri[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptri[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptri[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptri[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptri[i]; break; 
	  }
      break;
    case ND_ULONG:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrul[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrul[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrul[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrul[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrul[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrul[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrul[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrul[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrul[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrul[i]; break; 
	  }
      break;
    case ND_LONG:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrl[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrl[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrl[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrl[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrl[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrl[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrl[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrl[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrl[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrl[i]; break; 
	  }
      break;
    case ND_FLOAT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrf[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrf[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrf[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrf[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrf[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrf[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrf[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrf[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrf[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrf[i]; break; 
	  }
	break;
    case ND_DOUBLE:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) newuc[i]=(unsigned char)ptrd[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) newc[i]=(char)ptrd[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) newus[i]=(unsigned short)ptrd[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) news[i]=(short)ptrd[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) newui[i]=(unsigned int)ptrd[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) newi[i]=(int)ptrd[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) newul[i]=(unsigned long)ptrd[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) newl[i]=(long)ptrd[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) newf[i]=(float)ptrd[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) newd[i]=(double)ptrd[i]; break; 
	  }
      break;
    }
  
  

  if (sizeof_NDfield(type) > sizeof_NDfield(field->datatype))
    {
      free(field->val);
      field->val = new_tab;
    }
  else if (sizeof_NDfield(type) < sizeof_NDfield(field->datatype))
    field->val = realloc(field->val,field->nval*sizeof_NDfield(type));

  field->datatype=type;

  return 0;
}

// delta_bbox and x0_bbox are the coordinates of the GLOBAL bounding box
// This can be NULL if periodic=0
// result stored in x0_cut and x0_delta
int NDIntersection(double *x0_a,double *delta_a,double *x0_b,double *delta_b, double *x0_bbox,double *delta_bbox, double *x0_cut,double *delta_cut, int ndims, int periodic)
{
  int i,j;
  double xmin[ndims];
  double xmax[ndims];
  double x[ndims];
  int tmpa[ndims];
  int tmpb[ndims];

  for (j=0;j<ndims;j++) { xmin[j]=1.E40; xmax[j]=-1.E40;}

  //printf ("a:[%f,%f][%f,%f]\n",x0_a[0],x0_a[1],x0_a[0]+delta_a[0],x0_a[1]+delta_a[1]);
  //printf ("b:[%f,%f][%f,%f]\n",x0_b[0],x0_b[1],x0_b[0]+delta_b[0],x0_b[1]+delta_b[1]);
  
  /*
  if (periodic)
    for (j=0;j<ndims;j++)
      {
	if (periodic&(1<<j))
	  {
	    tmpa[j]=0;tmpb[j]=0;
	    if (x0_a[j]<x0_bbox[j])
	      tmpa[j] = (int)(((x0_bbox[j]+delta_bbox[j])-x0_a[j])/delta_bbox[j]);
	    else if (x0_a[j]>x0_bbox[j]+delta_bbox[j])
	      tmpa[j] = -(int)((x0_a[j]-x0_bbox[j])/delta_bbox[j]);

	    if (x0_b[j]<x0_bbox[j])
	      tmpb[j] = (int)(((x0_bbox[j]+delta_bbox[j])-x0_b[j])/delta_bbox[j]);
	    else if (x0_b[j]>x0_bbox[j]+delta_bbox[j])
	      tmpb[j] = -(int)((x0_b[j]-x0_bbox[j])/delta_bbox[j]);

	    if (tmpa[j]) x0_a[j]+=tmpa[j]*delta_bbox[j];
	    if (tmpb[j]) x0_b[j]+=tmpb[j]*delta_bbox[j];
	  }
      }
  */
  if (periodic)
    for (j=0;j<ndims;j++)
      {
	if (periodic&(1<<j))
	  {
	    tmpa[j]=0;tmpb[j]=0;
	    if (x0_a[j]<x0_bbox[j])
	      tmpa[j] = (int)(((x0_bbox[j]+delta_bbox[j])-x0_a[j])/delta_bbox[j]);
	    else if (x0_a[j]>x0_bbox[j]+delta_bbox[j])
	      tmpa[j] = -(int)((x0_a[j]-x0_bbox[j])/delta_bbox[j]);

	    if (tmpa[j]) x0_a[j]+=tmpa[j]*delta_bbox[j];
	    //printf ("%f %f, %f %d\n",x0_a[j],x0_b[j],(x0_a[j]-x0_b[j]),(int)(x0_a[j]-x0_b[j]));

	    if (x0_b[j]+delta_b[j]<x0_a[j]) 
	      tmpb[j] = (int)((x0_a[j]-x0_b[j]-delta_bbox[j])/delta_bbox[j]) +1;
	      
	    else if (x0_a[j]+delta_a[j]<x0_b[j])
	      tmpb[j] = (int)((x0_a[j]-x0_b[j]+delta_bbox[j])/delta_bbox[j]) -1;
	      
	  
	     if (tmpb[j]) x0_b[j]+=tmpb[j]*delta_bbox[j];
	  }
      }
  //printf ("a2:[%f,%f][%f,%f]\n",x0_a[0],x0_a[1],x0_a[0]+delta_a[0],x0_a[1]+delta_a[1]);
  //printf ("b2:[%f,%f][%f,%f]\n",x0_b[0],x0_b[1],x0_b[0]+delta_b[0],x0_b[1]+delta_b[1]);
  
  for (i=0;i<(1<<ndims);i++)
    for (j=0;j<ndims;j++)
      {
	x[j]=x0_a[j];
	if (i&(1<<j)) x[j] += delta_a[j]; 
	
	if ((x[j]>=x0_b[j])&&(x[j]<=x0_b[j]+delta_b[j]))
	  {
	    if (x[j]<xmin[j]) xmin[j]=x[j];
	    if (x[j]>xmax[j]) xmax[j]=x[j];
	  }
	
	x[j]=x0_b[j];
	if (i&(1<<j)) x[j] += delta_b[j]; 
	
	if ((x[j]>=x0_a[j])&&(x[j]<=x0_a[j]+delta_a[j]))
	  {
	    if (x[j]<xmin[j]) xmin[j]=x[j];
	    if (x[j]>xmax[j]) xmax[j]=x[j];
	  }
      }

  if (periodic)
    for (j=0;j<ndims;j++)
      if (periodic&(1<<j))
	{
	  if (tmpa[j]) x0_a[j]-=tmpa[j]*delta_bbox[j];
	  if (tmpb[j]) x0_b[j]-=tmpb[j]*delta_bbox[j];
	  //printf ("%d %d\n",tmpa[j],tmpb[j]);
	}
 
  for (j=0;j<ndims;j++) if (xmin[j]>xmax[j]) return 0;
  for (j=0;j<ndims;j++) {x0_cut[j]=xmin[j];delta_cut[j]=xmax[j]-xmin[j];}
  //for (j=0;j<ndims;j++) if ((tmpa[j])||(tmpb[j])) printf ("PERIODICITY USEDDDD !!!!!\n");
  return 1;
}


NDfield *Load_NDfield_ASCII(char *filename)
{
  int verbose=1;
  long i;      
  FILE *f = fopen(filename,"r");
  if (f==NULL) return NULL;
  char *line=NULL;
  int n;
  long dummy;
  char comment[80];

  strcpy(comment,"");
      
  Mygetline(&line,&n,f);  
  if (strstr(line,NDFIELD_ASCII_TAG)==NULL) 
    {
      free(line);
      return NULL;
    }   
      
  int coord=0;
  char *str[1000];
      
  if (strstr(line,"COORDS")!=NULL)
    coord=1;
      
  Mygetline(&line,&n,f);
    
  if ((strstr(line,"[")==NULL)||(strstr(line,"]")==NULL))
    {
      fprintf(stderr,"ERROR trying to load file '%s'.\n",filename);
      fprintf(stderr,"   Not a valid ASCII NDfield format.\n");
      fprintf(stderr,"   --> could not read dimensions.\n");
      free(line);
      return NULL;
    }

  if (strstr(line,"BBOX"))
    {
      fprintf(stderr,"ERROR trying to load file '%s'.\n",filename);
      fprintf(stderr,"   bounding box must be defined AFTER dimensions.\n");	  
      free(line);
      exit(-1);
    }
  if (line[0]=='#')
    {
      fprintf(stderr,"ERROR trying to load file '%s'.\n",filename);
      fprintf(stderr,"   comments must be defined AFTER dimensions.\n");	  
      free(line);
      exit(-1);
    }

  char *c=strstr(line,"[");
  int dims[NDFIELD_MAX_DIMS];
  int dimsCount=0;
  while (*c!=']')
    {
      if (!isdigit(*c)) c++;
      else
	{
	  int v;
	  sscanf(c,"%d",&v);
	  dims[dimsCount++]=v;
	  while(isdigit(*c)) c++;
	}
    };

  int ndims=coord?dims[0]:dimsCount;
  long nval=1;
  for (i=0;i<dimsCount;i++) nval*=dims[i];
  double x0[NDFIELD_MAX_DIMS];
  double delta[NDFIELD_MAX_DIMS];
  int nX0=0;
  int nDelta=0;
  //std::vector<double> x0;
  //std::vector<double> delta;
  dummy=Mygetline(&line,&n,f);
  if (line[0]=='#')
    {
      if (dummy>80) line[79]='\0';
      strcpy(comment,line);
      dummy=Mygetline (&line,&n,f);
    }

  if ((c=strstr(line,"BBOX"))!=NULL)
    {
      int ntok=str2tok2(c,"[ ],;|\t",7,0,0,str);
      if (ntok != ndims*2+1)
	{
	  fprintf(stderr,"ERROR trying to load file '%s'.\n",filename);
	  fprintf(stderr,"  Bad definition of the bounding box, should have ndims=%d\n",ndims);
	  exit(-1);
	}
      for (i=1;i<ntok;i++)
	{
	  if (nX0<ndims) x0[nX0++]=strtod(str[i],NULL);
	  //if ((long)x0.size()<ndims) x0.push_back(strtod(str[i],NULL));
	  else delta[nDelta++]=strtod(str[i],NULL);
	}
      Mygetline(&line,&n,f);
    }
      
  char dimstr[256];
  sprintf(dimstr,"[%d",dims[0]);
  for (i=1;i<(long)dimsCount;i++) sprintf(dimstr,"%s,%d",dimstr,dims[i]);
  strcat(dimstr,"]");
  if (verbose>=1) printf("Reading %s %s from NDfied ASCII file '%s' ... ",dimstr,coord?"coords":"grid",filename);
  fflush(0);
      
  double *tab=(double*)malloc(sizeof(double)*nval);
  long nread=0;
  do {
	
    //Mygetline(&line,&n,f);
    //printf ("line : %s\n",line);
    if (feof(f)) break;
    //printf ("line : %s\n",line);
    int ntok=str2tok2(line," ,;|[]\t",7,'#',0,str);
    if ((nread+ntok) > nval) 
      {
	ntok=nval-nread;
	fprintf(stderr,"WARNING: Too many values to read,stopping.\n");
      }
    for (i=0;i<ntok;i++) tab[nread++]=strtod(str[i],NULL);
    Mygetline(&line,&n,f);
  } while (nread<nval);
  free(line);

  if (nread!=nval)
    {
      fprintf(stderr,"ERROR: too few values to read, aborting. (%ld/%ld)\n",nread,nval);
      free(tab);
      return NULL;
    }

  if (!nX0)
    {
      if (!coord)
	{
	  int i;
	  for (i=0;i<ndims;++i)
	    {
	      x0[i]=ndims;
	      delta[i]=dims[i];
	    }
	  /* x0.assign(ndims,0); */
	  /* delta.assign(dims.begin(),dims.end()); */
	}
      else
	{
	  long j;
	  int i;
	  for (i=0;i<ndims;++i)
	    {
	      x0[i]=1.E30;
	      delta[i]=-1.E30;
	    }
	  /* x0.assign(ndims,std::numeric_limits<double>::max()); */
	  /* delta.assign(ndims,-std::numeric_limits<double>::max()); */
	  for (i=0;i<nval;i+=dims[0])
	    for (j=0;j<dims[0];j++)
	      {
		if (tab[i+j]<x0[j]) x0[j]=tab[i+j];
		if (tab[i+j]>delta[j]) delta[j]=tab[i+j];
	      }
	      
	  for (j=0;j<dims[0];j++) delta[j]-=x0[j];
	      
	}
    }
  NDfield *field=Create_NDfield(&dims[0],ndims,coord?1:0,ND_DOUBLE,&x0[0],&delta[0],(void*)(tab),comment);

  if (verbose>=1) printf("done. \n");
  return field;
}

NDfield *Load_NDfield(char *filename)
{
  int type=IsNDfield(filename);
  if (type==1) return Load_NDfield_BIN(filename);
  if (type==2) return Load_NDfield_ASCII(filename);

  fprintf(stderr,"ERROR:File %s is not an NDfield file !\n",filename);
  return 0;
}


/*
  NDfield *ComputeNodeField(NDfield *field,int periodic)
  {
  //int *dims=field->dims;
  //int ndims=field->ndims;

  INT coord[field->ndims];
  INT dm[field->ndims];
  INT dp[field->ndims];
  INT dx[field->ndims+1];
  INT delta;
  int i,j,k;
  //INT tmp[3];
  
  FLOAT *val=field->val;
  FLOAT *new_val;

  NDfield *new_field;

  new_field=Create_NDfield(field->dims,field->ndims,field->fdims_index,field->datatype,field->x0, field->delta, NULL, "Node value");
  new_val=(FLOAT *)new_field->val;

  memset(coord,0,sizeof(INT)*field->ndims);
  memset(new_val,0,field->nval*sizeof(FLOAT));
  coord[0]=-1;

  //dims=field->dims;
  //ndims=field->ndims;

  for (i=0;i<field->nval;i++)
  {
  for (j=0,coord[0]++;coord[j]>=field->dims[j];coord[j]=0,coord[j+1]++,j++) {}
      
  dx[0]=1;
      
  for (j=0;j<field->ndims;j++)
  {
  dx[j+1]=dx[j]*field->dims[j];
	  
  if (coord[j]+1>=field->dims[j]) dp[j]=dx[j]-dx[j+1];
  else dp[j]=dx[j];
	  
  if (coord[j]-1<0) dm[j]=dx[j+1]-dx[j];
  else dm[j]=-dx[j];
  }
      
  for (j=0;j<(1<<field->ndims);j++)
  {
  delta=i;
  for (k=0;k<field->ndims;k++)
  {
  if (j&(1<<k)) delta += dm[k];
  }
  new_val[i]+=val[delta];
  }
  new_val[i]/=(1<<field->ndims);
  }

  return new_field;
  }
*/

