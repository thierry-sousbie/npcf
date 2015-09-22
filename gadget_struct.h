#ifndef __Struct_IO__
#define __Struct_IO__


#ifdef __cplusplus
extern "C" {
#endif


typedef struct snapshot_header_str
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];  /* fills to 256 Bytes */
} snapshot_header;


typedef struct snapshot_data_str 
{
    snapshot_header header;
    float  *Pos;
    float  *Vel;
    float  *Mass;
    int    *Id;
    char   *Type;
    
    float  *Rho, *U, *Temp, *Ne;
    int N;
    int flags;
	
    int nSupFields;
    float **supField;
    char **supFieldName;
} snapshot_data;


#ifdef __cplusplus
}
#endif


#endif
