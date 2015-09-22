#include <stdio.h>
#include "mystring.h"

//#ifndef _GNU_SOURCE
/* Default value for line length.  */
static const int line_size = 1024;

ssize_t
Mygetdelim (char **lineptr, int *n, int delim, FILE *stream)
{
  int indx = 0;
  int c;

  /* Sanity checks.  */
  if (lineptr == NULL || n == NULL || stream == NULL)
    return -1;

  /* Allocate the line the first time.  */
  if (*lineptr == NULL)
    {
      *lineptr = malloc (line_size);
      if (*lineptr == NULL)
        return -1;
      *n = line_size;
    }

  /* Clear the line.  */
  memset (*lineptr, '\0', *n);

  while ((c = getc (stream)) != EOF)
    {
      /* Check if more memory is needed.  */
      if (indx >= *n)
        {
          *lineptr = realloc (*lineptr, *n + line_size);
          if (*lineptr == NULL)
            {
              return -1;
            }
          /* Clear the rest of the line.  */
          memset(*lineptr + *n, '\0', line_size);
          *n += line_size;
        }

      /* Push the result in the line.  */
      (*lineptr)[indx++] = c;

      /* Bail out.  */
      if (c == delim)
        {
          break;
        }
    }
  return (c == EOF) ? -1 : indx;
}
//#endif


ssize_t Mygetline (char **lineptr, int *n, FILE *stream)
{
  //#ifndef _GNU_SOURCE
  return Mygetdelim (lineptr, n, '\n', stream);
  //#else
  //return getline(lineptr,  (size_t *)n,  stream);
  //#endif

}


int str2tok(char *str, char sep, char** tok)
{
    int n=0;
    int i;
    int len = strlen(str);
    int last=0;
    
    for (i=0;i<len;i++)
    {
	if (str[i]==sep) 
	    str[i]='\0';
    }

    for (i=0;i<len;i++)
    {
	if (str[i]!='\0')
	{
	    if (last == 0) 
		tok[n++] = &(str[i]);
	    last=1;
	}
	else
	    last=0;
    }
    return n;
}

int str2tok2(char *str,const char *sep, int nsep, char comment, char ignore, char** tok)
{
    int n=0;
    int i,j;
    int len = strlen(str);
    int last=0;
    
    for (i=0;i<len;i++)
    {
	for (j=0;j<nsep;j++)
	    if (str[i]==sep[j])
	    {
		str[i]='\0';
		continue;
	    }
	    else if (comment&(str[i]==comment))
	    {
		len=i;
		continue;
	    }
    }

    for (i=0;i<len;i++)
    {
	if ((str[i]!='\0')&&(str[i]!='\n')&&(str[i]!=ignore))
	{
	    if (last == 0) 
		tok[n++] = &(str[i]);
	    last=1;
	}
	else
	    last=0;
    }
  
    return n;
}

char *CutName(char *Name)
{
    int i,j; 

    for (i=0,j=-1;i<strlen(Name);i++) 
	if (Name[i]=='/') j=i;

    if (j!=-1)
	return &Name[j+1];
    else
	return Name;
}
