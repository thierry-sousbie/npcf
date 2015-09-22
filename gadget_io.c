#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <xlocale.h>
#include "gadget_io.h"
#include "mystring.h"

//checks if file is a gadget file
int IsGadgetFile(const char *fname)
{
    FILE *fd;
    int i;
	char *line=NULL;
   
    //check if file exists
    if(!(fd=fopen(fname,"r")))
	return 1;

    fread(&i,sizeof(int),1,fd);

    if ((i!=256)&&(swapI(i)!=256))
    {
	fclose(fd);
	fd=fopen(fname,"r");
	Mygetline(&line,&i,fd);
	fclose(fd);
	if ((strstr(line,"x")!=NULL)&&(strstr(line,"y")!=NULL)&&(strstr(line,"z")!=NULL))
	{
	    free(line);
	    return 2;
	}

	free(line);
	return 0;
	}
	

    fclose(fd);
    
    return 1;
}

int ReadASCII2Gadget(const char *fname,snapshot_data *P)
{
int npart;
FILE *fd;
int i,j;
char *line=NULL;
int linesize=0;
char *token[50];
int nfields=0;
char fieldname[50][25];
int indexof[50];

int idid=-1,idx=-1,idy=-1,idz=-1,idvx=-1,idvy=-1,idvz=-1,idmass=-1;

printf("reading %s (ASCII) ...",fname);fflush(0);
memset(P,0,sizeof(snapshot_data));
if(!(fd=fopen(fname,"r")))
	return 0;

Mygetline(&line,&linesize,fd);
nfields=str2tok(line,' ',token);
for (i=0;i<nfields;i++) strcpy(fieldname[i],token[i]);

P->nSupFields=0;
P->Pos=NULL;P->Vel=NULL;
P->Id=NULL;

for (i=0;i<nfields;i++)
{
    indexof[i]=-1;

    if (!strncmp(fieldname[i],"x",1)) idx=i;
    else if (!strncmp(fieldname[i],"y",1)) idy=i;
    else if (!strncmp(fieldname[i],"z",1)) idz=i;
    else if (!strncmp(fieldname[i],"vx",2)) idvx=i;
    else if (!strncmp(fieldname[i],"vy",2)) idvy=i;
    else if (!strncmp(fieldname[i],"vz",2)) idvz=i;
    else if (!strncmp(fieldname[i],"id",2)) idid=i;
    //else if (!strncmp(fieldname[i],"mass",4)) idmass=i;
    else {
	indexof[i]=P->nSupFields;
	P->nSupFields++;
    }
}
P->supField = calloc (P->nSupFields , sizeof(float*));

npart=0;
if (idmass!=-1) P->Mass = malloc((npart+10000)*sizeof(float));
if (idid!=-1) P->Id = malloc((npart+10000)*sizeof(int));
if ((idvx!=-1)||(idvy!=-1)||(idvz!=-1)) P->Vel = malloc(3*(npart+10000)*sizeof(float));
P->Pos = malloc(3*(npart+10000)*sizeof(float));
for (j=0;j<P->nSupFields;j++)
    P->supField[j] = malloc((npart+10000)*sizeof(float));

Mygetline(&line,&linesize,fd);
do
{
    if (nfields != str2tok(line,' ',token))
    {
	printf ("Line %d is incomplete, SKIPPED.\n",npart+2);
	continue;
    }
    for (i=0;i<nfields;i++)
    {	
	if (indexof[i]!=-1) P->supField[indexof[i]][npart] = atof(token[i]);
	else if (i==idx) P->Pos[3*npart]= atof(token[i]);
	else if (i==idy) P->Pos[3*npart+1]= atof(token[i]);
	else if (i==idz) P->Pos[3*npart+2]= atof(token[i]);
	else if (i==idvx) P->Vel[3*npart]= atof(token[i]);
	else if (i==idvy) P->Vel[3*npart+1]= atof(token[i]);
	else if (i==idvz) P->Vel[3*npart+2]= atof(token[i]);
	else if (i==idmass) P->Mass[npart]= atof(token[i]);
	else if (i==idid) P->Id[npart]= atoi(token[i]);
    }
    npart++;
    if (npart%10000 == 0)
    {
	if (idmass!=-1) P->Mass = realloc(P->Mass,(npart+10000)*sizeof(float));
	if (idid!=-1) P->Id = realloc(P->Id,(npart+10000)*sizeof(int));
	if ((idvx!=-1)||(idvy!=-1)||(idvz!=-1)) P->Vel = realloc(P->Vel,3*(npart+10000)*sizeof(float));
	P->Pos = realloc(P->Pos,3*(npart+10000)*sizeof(float));
	for (j=0;j<P->nSupFields;j++)
	    P->supField[j] = realloc(P->supField[j],(npart+10000)*sizeof(float));
	printf("\rreading %s (ASCII) ... (%d lines)",fname,npart);fflush(0);
    }
    
    Mygetline(&line,&linesize,fd);
    
} while (!feof(fd));

if (idmass!=-1) P->Mass = realloc(P->Mass,(npart)*sizeof(float));
if (idid!=-1) P->Id = realloc(P->Id,(npart)*sizeof(int));
if ((idvx!=-1)||(idvy!=-1)||(idvz!=-1)) P->Vel = realloc(P->Vel,3*(npart)*sizeof(float));
P->Pos = realloc(P->Pos,3*(npart)*sizeof(float));
for (j=0;j<P->nSupFields;j++)
    P->supField[j] = realloc(P->supField[j],(npart)*sizeof(float));

P->header.npart[1]=P->header.npartTotal[1]=npart;
if (idmass == -1) P->header.mass[1]=1;
P->header.BoxSize=0;

free(line);
fclose(fd);

printf("\rreading %s (ASCII) ... done. (%d lines)\n",fname,npart);fflush(0);

return 1;
}

/* POINTERS in P MUST BE SET TO NULL BEFORE CALLING */
/* OR TO PREVIOUSLY ALLOCATED DATA  */
/* SHOULD CGANGE REALLOC-> MALLOC OR CALLOC */
/*
 * reads gadget format (multiple files) 
 * if filename is "" stdin is read
 */
int ReadGadget(const char *fname, snapshot_data *P, int flags_p)
{
    FILE *fd;
    char buf[200];
    long i,k,offs,ntot_withmasses,nr;
    long n,pc,pc_new,pc_sph;
    int Ngas;
    long NumPart=0;
    int files;
    int tmp;
    int flags = flags_p;

    if (!IsGadgetFile(fname))
      {
	fprintf(stderr,"%s is not a gadget file, will try SIMPLE format \n",fname);
	return 0;
      }
	  
    if (IsGadgetFile(fname)==2)
    {
	ReadASCII2Gadget(fname,P);
	return flags;
    }
    //guess how many files there are
    sprintf(buf,"%s",fname);
    if(!(fd=fopen(buf,"r")))
    {
	sprintf(buf,"%s.%d",fname,0);
	files=0;
	while ((fd=fopen(buf,"r")))
	{
	    fclose(fd);
	    sprintf(buf,"%s.%d",fname,++files);
	}
	if(files<2)
	{
	    printf("can't open file '%s' nor `%s`\n",fname,buf);
	    exit(0);
	}
    }
    else 
    {
	files=1;
	fclose(fd);
    }

    for(i=0, pc=0; i<files; i++, pc=pc_new)
    {
	if(files>1)
	    sprintf(buf,"%s.%d",fname,(int)i);
	else
	    sprintf(buf,"%s",fname);
	
	if (strlen(buf)==0) fd=stdin;
	else
	    if(!(fd=fopen(buf,"r")))
	    {
		printf("can't open file `%s`\n",buf);
		exit(0);
	    }
	
	printf("reading %s ...",buf); fflush(stdout);
	
	//SKIP(fd);
	//Guess if endianness changes or not 
	{
	    int skpdummy;
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    if (skpdummy==256) 
	    {
		flags &= (~FLAG_SWAPENDIAN);
	    }
		else
	    {
		flags |= FLAG_SWAPENDIAN;
		printf ("(Endianness differs)");fflush(0);
	    }
	}
		
		P->nSupFields = 0;
		P->supField = NULL;

	fread(&P->header, sizeof(P->header), 1, fd);
	//eventually swap endianness
	if (flags&FLAG_SWAPENDIAN)
	{
	    Dswap4BArr(P->header.npart,6);
	    Dswap4BArr(P->header.npartTotal,6);
	    Dswap8BArr(P->header.mass,6);

	    Dswap4B(&P->header.flag_sfr);Dswap4B(&P->header.flag_feedback);Dswap4B(&P->header.flag_cooling);
	    Dswap4B(&P->header.num_files);
	    Dswap8B(&P->header.time);Dswap8B(&P->header.redshift);Dswap8B(&P->header.BoxSize);
	    Dswap8B(&P->header.Omega0);Dswap8B(&P->header.OmegaLambda);Dswap8B(&P->header.HubbleParam);
	}
	SKIP(fd);
	
	if(files==1)
	{
	    for(k=0, NumPart=0, ntot_withmasses=0; k<6; k++)
		NumPart+= (long)P->header.npart[k];
	    Ngas= P->header.npart[0];
	}
	else
	{
	    for(k=0, NumPart=0, ntot_withmasses=0; k<6; k++)
		NumPart+= (long)P->header.npartTotal[k];
	    Ngas= P->header.npartTotal[0];
	}
	
      for(k=0, ntot_withmasses=0; k<6; k++)
	{
	  if(P->header.mass[k]==0)
	      ntot_withmasses+= (long)P->header.npart[k];
	}

      P->N = NumPart;
      
      if(i==0)
      {
	  if (flags&FLAG_POS)
	      if(!(P->Pos=(float *)realloc(P->Pos,(long)3*(long)NumPart*(long)sizeof(float))))
	      {
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
	      }
	  
	  if (flags&FLAG_VEL)
	      if(!(P->Vel=(float *)realloc(P->Vel,(long)3*(long)NumPart*(long)sizeof(float))))
	      {
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
	      }

	  if (flags&FLAG_ID)
	      if(!(P->Id=(int *)realloc(P->Id,(long)NumPart*(long)sizeof(int))))
	      {
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
	      }

	  if (flags&FLAG_TYPE)
	      if(!(P->Type=(char *)realloc(P->Type,NumPart*sizeof(char))))
	      {
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
	      }
	  

	  if (P->header.npartTotal[0])
	  {
	      if (flags&FLAG_GAS)
	      {
			  if(!(P->Rho=(float *)realloc(P->Rho,P->header.npartTotal[0]*sizeof(float))))
			  {
				  fprintf(stderr,"failed to allocate memory.\n");
				  exit(0);
			  }
			  
			  if(!(P->U=(float *)realloc(P->U,P->header.npartTotal[0]*sizeof(float))))
			  {
				  fprintf(stderr,"failed to allocate memory.\n");
				  exit(0);
			  }
			  
			  if(!(P->Temp=(float *)realloc(P->Temp,P->header.npartTotal[0]*sizeof(float))))
			  {
				  fprintf(stderr,"failed to allocate memory.\n");
				  exit(0);
			  }
			  
			  if(!(P->Ne=(float *)realloc(P->Ne,P->header.npartTotal[0]*sizeof(float))))
			  {
				  fprintf(stderr,"failed to allocate memory.\n");
				  exit(0);
			  }
	      }
	  }
	  else
	  {
	      P->Rho=NULL;
	      P->U=NULL;
	      P->Temp=NULL;
	      P->Ne=NULL;
	  }

	  if ((ntot_withmasses)&&(flags&FLAG_MASS))
	      if(!(P->Mass=(float *)realloc(P->Mass,NumPart*sizeof(float))))
	      {
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
	      }

      }
	  
      SKIP(fd);
      if ((flags&FLAG_POS)&&(feof(fd))) 
      {
	  flags &= ~FLAG_POS;
	  free(P->Pos);P->Pos=NULL;
      }
      for(k=0,tmp=0;k<6;k++)
      {
	   
	  if (flags&FLAG_POS)
	  {
	      //printf ("Hello %d\n",(size_t)P->header.npart[k]);
	    //printf("read %d \n",fread(&P->Pos[(long)3*pc], 1+0*3*sizeof(float), (size_t)P->header.npart[k], fd));
	    //printf("read %d \n",fread(&P->Pos[(long)3*pc], 1+0*3*sizeof(float), (size_t)P->header.npart[k], fd));
	    nr = 16;
	    if (P->header.npart[k]!=0)
	    {
	      if (P->header.npart[k]<150000000)
		fread(&P->Pos[(long)3*pc+tmp],3*sizeof(float),P->header.npart[k],fd);
	      else
		for (offs=0;offs<nr;offs++)
		  {
		    //fseek(fd,(long)P->header.npart[k]*off,SEEK_CUR);
		    fread(&P->Pos[(long)3*pc +tmp+ 3*(long)offs*(long)P->header.npart[k]/nr],3*sizeof(float),P->header.npart[k]/nr,fd);
		  }
	    }

	    if (flags&FLAG_SWAPENDIAN)
	      Dswap4BArr(&P->Pos[(long)3*pc+tmp],(long)3*(long)P->header.npart[k]);

	    tmp+=(long)3*(long)P->header.npart[k];
	  }
	  else
	      fseek(fd,(long)3*(long)P->header.npart[k]*sizeof(float),SEEK_CUR);
      }

      SKIP(fd);

      SKIP(fd);
      if ((flags&FLAG_VEL)&&(feof(fd)))
      {
	  flags &= ~FLAG_VEL;
	  free(P->Vel);P->Vel=NULL;
      }
      for(k=0,tmp=0;k<6;k++)
      {
	  if (flags&FLAG_VEL)
	  {
	    nr = 16;
	    if (P->header.npart[k]<150000000)
	      fread(&P->Vel[(long)3*pc+tmp],3*sizeof(float),P->header.npart[k],fd);
	    else
	      for (offs=0;offs<nr;offs++)
		{
		  //fseek(fd,(long)P->header.npart[k]*off,SEEK_CUR);
		  fread(&P->Vel[(long)3*pc +tmp+ 3*(long)offs*(long)P->header.npart[k]/nr],3*sizeof(float),P->header.npart[k]/nr,fd);
		}
	    
	      if (flags&FLAG_SWAPENDIAN)
		  Dswap4BArr(&P->Vel[(long)3*pc+tmp],(long)3*(long)P->header.npart[k]);
	      tmp+=(long)3*(long)P->header.npart[k];
	  }
	  else
	      fseek(fd,(long)3*(long)P->header.npart[k]*sizeof(float),SEEK_CUR);
      }
      SKIP(fd);

      SKIP(fd);
      if ((flags&FLAG_ID)&&(feof(fd)))
      {
	  flags &= ~FLAG_ID;
	  free(P->Id);P->Id=NULL;
      }
      for(k=0,tmp=0;k<6;k++)
      {
	  if (flags&FLAG_ID)
	  {
	      //printf ("pc=%d\n",pc);
	    nr=16;
	    if (P->header.npart[k]<400000000)
	      fread(&P->Id[(long)pc+tmp],sizeof(int),P->header.npart[k],fd);
	    else
	      for (offs=0;offs<nr;offs++)
		{
		  //fseek(fd,(long)P->header.npart[k]*off,SEEK_CUR);
		  fread(&P->Id[(long)pc +tmp+ (long)offs*(long)P->header.npart[k]/nr],sizeof(int),P->header.npart[k]/nr,fd);
		}
	    
	      if (flags&FLAG_SWAPENDIAN)
		  Dswap4BArr(&P->Id[pc+tmp],(long)P->header.npart[k]);
	      tmp+=P->header.npart[k];
	  }
	  else
	      fseek(fd,(long)P->header.npart[k]*sizeof(int),SEEK_CUR);
      }
      SKIP(fd);

      if(ntot_withmasses>0)
	  SKIP(fd);

		for(k=0, pc_new=pc; k<6; k++)
		{
			for(n=0;n<(long)P->header.npart[k];n++)
			{
				if (flags&FLAG_TYPE)
					P->Type[pc_new]=(char)k;
				
				if (flags&FLAG_MASS)
				{
					if (P->header.mass[k]==0)
					{
						fread(&P->Mass[pc_new], sizeof(float), 1, fd);
						if (flags&FLAG_SWAPENDIAN)
							Dswap4B(&P->Mass[pc_new]);
					}
					else if(ntot_withmasses>0)
						P->Mass[pc_new]= P->header.mass[k];
				}
				else if (P->header.mass[k]==0)
					fseek(fd,sizeof(float),SEEK_CUR);
				
				pc_new++;
			}
		}
		
		if(ntot_withmasses>0)
	  SKIP(fd);

      if ((P->header.npart[0]>0)&&(flags&FLAG_GAS))
	  {
		  SKIP(fd);
		  fread(&P->U[pc], sizeof(float), P->header.npart[0], fd);
		  if (flags&FLAG_SWAPENDIAN)
			  Dswap4BArr(&P->U[pc],P->header.npart[0]);
		  SKIP(fd);
		  
		  SKIP(fd);
		  fread(&P->Rho[pc], sizeof(float), P->header.npart[0], fd);
		  if (flags&FLAG_SWAPENDIAN)
			  Dswap4BArr(&P->Rho[pc],P->header.npart[0]);
		  SKIP(fd);
		  
		  SKIP(fd);
		  fread(&P->Temp[pc], sizeof(float), P->header.npart[0], fd);
		  if (flags&FLAG_SWAPENDIAN)
			  Dswap4BArr(&P->Temp[pc],P->header.npart[0]);
		  SKIP(fd);
		  
		  if(P->header.flag_cooling)
		  {
			  SKIP(fd);
			  fread(&P->Ne[pc], sizeof(float), P->header.npart[0], fd);
			  if (flags&FLAG_SWAPENDIAN)
				  Dswap4BArr(&P->Ne[pc],P->header.npart[0]);
			  SKIP(fd);
		  }
		  else
			  for(n=0, pc_sph=pc; n<P->header.npart[0];n++)
			  {
				  P->Ne[pc_sph]= 1.0;
				  pc_sph++;
			  }
	  }
		
		P->nSupFields = 0;
		SKIP(fd);
		while ((!feof(fd))&&(flags&FLAG_SUP)) {
			
			if (i==0)
			{
				P->supField = realloc (P->supField,sizeof(float*)*(P->nSupFields+1));
				P->supField[P->nSupFields] = malloc(sizeof(float)*P->N);
			}
			
			for(k=0,tmp=0;k<6;k++)
			{

					//printf ("pc=%d\n",pc);
					nr=16;
					if (P->header.npart[k]<400000000)
						fread(&(P->supField[P->nSupFields][(long)pc+tmp]),sizeof(float),P->header.npart[k],fd);
					else
						for (offs=0;offs<nr;offs++)
						{
							//fseek(fd,(long)P->header.npart[k]*off,SEEK_CUR);
							fread(&(P->supField[P->nSupFields][(long)pc +tmp+ (long)offs*(long)P->header.npart[k]/nr]),sizeof(float),P->header.npart[k]/nr,fd);
						}
					
					if (flags&FLAG_SWAPENDIAN)
						Dswap4BArr(&(P->supField[P->nSupFields][pc+tmp]),(long)P->header.npart[k]);
					tmp+=P->header.npart[k];

			}			
			SKIP(fd);SKIP(fd);
			P->nSupFields++;
			
		}
		
		fclose(fd);
		printf(" done.\n"); fflush(stdout);
    }
	
    return flags;
    

}

void freegadgetstruct(snapshot_data *snap)
{
    free(snap->Pos);
    free(snap->Vel);
    free(snap->Id);
    free(snap);
}

int WriteGadget(const char *fname, snapshot_data *gadget, int *selection, int selection_size)
{
    int i,index,n,ncur;
    int npart[6];
    int npartTotal[6];
    FILE *f;
    int dummy;
    int ntotal;

    memcpy(npart,gadget->header.npart,6*sizeof(int));
    memcpy(npartTotal,gadget->header.npartTotal,6*sizeof(int));
   
    if (selection_size)
    {
	index=i=ncur=0;
	n=gadget->header.npartTotal[index];
	
	do
	{
	    while (selection[i]>=3*n) 
	    {
		gadget->header.npartTotal[index] = ncur;
		n+=npartTotal[++index];
		ncur=0;
	    }
	    i++;
	    ncur++;
	} while (i<selection_size);
	
	for (;index<6;index++) {gadget->header.npartTotal[index] = ncur;ncur=0;}
	for (i=0;i<6;i++) gadget->header.npart[i] = gadget->header.npartTotal[i];
    }

    if (!selection_size)
    {
	for (i=0,ntotal=0;i<6;i++) ntotal+= gadget->header.npartTotal[i];
	for (i=0;i<6;i++) gadget->header.npart[i]= gadget->header.npartTotal[i];
    }
    else
	ntotal = selection_size;

    f=fopen(fname,"w");
    dummy=256;
    fwrite(&dummy,sizeof(int),1,f);
    fwrite(&(gadget->header),sizeof(char),256,f);
    fwrite(&dummy,sizeof(int),1,f);

    dummy = 3*ntotal*sizeof(float);
    if (gadget->Pos!=NULL)
    {
	fwrite(&dummy,sizeof(int),1,f);
	if (selection_size)
	{
	    for (i=0;i<selection_size;i++)
		fwrite(&(gadget->Pos[selection[i]]),sizeof(float),3,f);
	}
	else
	    fwrite(gadget->Pos,sizeof(float),3*ntotal,f);
	fwrite(&dummy,sizeof(int),1,f);
    }

    if (gadget->Vel!=NULL)
    {
	fwrite(&dummy,sizeof(int),1,f);
	if (selection_size)
	{
	    for (i=0;i<selection_size;i++)
		fwrite(&(gadget->Vel[selection[i]]),sizeof(float),3,f);
	}
	else
	    fwrite(gadget->Vel,sizeof(float),3*ntotal,f);
	fwrite(&dummy,sizeof(int),1,f);
    }

    dummy = ntotal*sizeof(int);
    if (gadget->Id!=NULL)
    {
	fwrite(&dummy,sizeof(int),1,f);
	if (selection_size)
	{
	    for (i=0;i<selection_size;i++)
		fwrite(&(gadget->Id[selection[i]/3]),sizeof(int),1,f);
	}
	else
	    fwrite(gadget->Id,sizeof(int),ntotal,f);
	fwrite(&dummy,sizeof(int),1,f);
    }

    if (gadget->header.npartTotal[0])
    {
	dummy = gadget->header.npartTotal[0]*sizeof(float);
	fwrite(&dummy,sizeof(int),1,f);
	
	if (selection_size)
	    for (i=0;i<gadget->header.npartTotal[0];i++)
		fwrite(&(gadget->U[selection[i]/3]), sizeof(float), 1, f);
	else
	    fwrite(gadget->U, sizeof(float), npartTotal[0], f);

	if (selection_size)
	    for (i=0;i<gadget->header.npartTotal[0];i++)
		fwrite(&(gadget->Rho[selection[i]/3]), sizeof(float), 1, f);
	else
	    fwrite(gadget->Rho, sizeof(float), npartTotal[0], f);

	if (selection_size)
	    for (i=0;i<gadget->header.npartTotal[0];i++)
		fwrite(&(gadget->Temp[selection[i]/3]), sizeof(float), 1, f);
	else
	    fwrite(gadget->Temp, sizeof(float), npartTotal[0], f);

	if (gadget->header.flag_cooling)
	{
	    if (selection_size)
		for (i=0;i<gadget->header.npartTotal[0];i++)
		    fwrite(&(gadget->Ne[selection[i]/3]), sizeof(float), 1, f);
	    else
		fwrite(gadget->Ne, sizeof(float), npartTotal[0], f);
	}

	fwrite(&dummy,sizeof(int),1,f);
    }

    fclose(f);

    memcpy(gadget->header.npart,npart,6*sizeof(int));
    memcpy(gadget->header.npartTotal,npartTotal,6*sizeof(int));

    return 0;
}

