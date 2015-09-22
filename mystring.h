#ifndef __MY_STRING_H__
#define __MY_STRING_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include <stdlib.h>

  ssize_t Mygetline (char **lineptr, int *n, FILE *stream);
  int str2tok(char *str, char sep, char** tok);
  int str2tok2(char *str,const char *sep, int nsep, char comment, char ignore, char** tok);
  char *CutName(char *Name);

#ifdef __cplusplus
}
#endif


#endif
