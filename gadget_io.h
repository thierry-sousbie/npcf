#ifndef __GADGET_IO
#define __GADGET_IO

#include <stdlib.h>
#include "gadget_struct.h"
#include "endian.h"

#ifdef __cplusplus
extern "C" {
#endif
    
#define FLAG_POS  ((int)1<<0)
#define FLAG_VEL  ((int)1<<1)
#define FLAG_MASS ((int)1<<2)
#define FLAG_ID   ((int)1<<3)
#define FLAG_TYPE ((int)1<<4)
#define FLAG_GAS  ((int)1<<5)
#define FLAG_SUP  ((int)1<<6)
#define FLAG_ALL  (((int)1<<7)-1)
#define FLAG_DM   ((int)FLAG_POS|FLAG_VEL|FLAG_MASS|FLAG_ID)
//if the file is not of the same endian type as the computer
//this is done automatically
#define FLAG_SWAPENDIAN ((int)1<<31)
    
#define SKIP(filename) {int skpdummy;fread(&skpdummy,sizeof(skpdummy),1,filename);}

    int IsGadgetFile(const char *);
    int ReadGadget(const char*,snapshot_data*,int);
    void freegadgetstruct(snapshot_data *snap);
    int WriteGadget(const char *fname, snapshot_data *gadget, int *selection, int selection_size);
    
#ifdef __cplusplus
}
#endif


#endif
