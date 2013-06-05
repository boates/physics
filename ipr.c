#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#define TRUE 1
#define FALSE 0

#define MaxChars 256

char *getline (char *s, int n, FILE *iop);
char *getline (char *s, int n, FILE *iop)
{
  register int c;
  register char *cs;
  
  cs = s;
  while (--n > 0 && (c = getc(iop)) != EOF && c != '\n' && c != '\r')
    *cs++ = c;
  
  *cs = '\0';
  return (c == EOF && cs == s) ? NULL : s;
}

int main(int argc, char *argv[])
{
  /* Read PARCHG and evaluate IPR

  IPR(psi_n)=N * Sum_1^N |psi_n|**4 / [Sum_1^N |psi_n|**2]**2
        N=nx*ny*nz

	The IPR is large for highly localized states and small for delocalized states
	-> a delta function would give N
	-> a constant function would give 1

  */
  FILE *inFile1 = NULL;
  char *comment;
  char *tok;
  double d1,d2,d3,d4,d5;
  double sum1,sum2;
  int i1, i2, i3;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  int numAtoms1;
  int numTypes1;
  long numPoints;
  int extraPoints = 0;
  long i;

  inFile1 = fopen (argv[1], "r");
  if (inFile1 == NULL)
    {
      printf ("Error opening file %s.\n", argv[1]);
      return 1;
    }
  comment = (char *) malloc (MaxChars);
  if (getline (comment, MaxChars, inFile1) == NULL)
    {  printf ("error reading comment line\n"); fclose(inFile1);
    free(comment); return 3;
    }
  comment[strlen(comment)] = '\0';/* strip '\r' from end of string */
  fscanf (inFile1, "%lg", &d1);
  //for (i = 0; i < 3; i++)
  //{
  //fscanf (inFile1, "%lg %lg %lg", &d1, &d2, &d3);
  //}
  fscanf (inFile1, "%lg %lg %lg", &a11, &a12, &a13);
  fscanf (inFile1, "%lg %lg %lg", &a21, &a22, &a23);
  fscanf (inFile1, "%lg %lg %lg", &a31, &a32, &a33);

  tok = (char *) malloc (MaxChars);
  getline (comment, MaxChars, inFile1); // first get end of last line
  getline (comment, MaxChars, inFile1);
  tok = strtok(comment, " \t");
  sscanf (tok, "%d", &numAtoms1);
  while ((tok = strtok(NULL, " \t")) != NULL)
    {
      sscanf (tok, "%d", &i1);
      numAtoms1 += i1;
    }
  free (tok);
  //comment = (char *) malloc (MaxChars);
  getline (comment, MaxChars, inFile1);
  comment[strlen(comment)] = '\0';/* strip '\r' from end of string */
  free (comment);
  
  for (i = 0; i < numAtoms1; i++)
    {
      fscanf (inFile1, "%lg %lg %lg", &d1, &d2, &d3);
    }
  fscanf (inFile1, "%d %d %d", &i1, &i2, &i3);

  numPoints = i1*i2*i3;
  extraPoints = numPoints%5;
  numPoints /= 5;
  sum1 = 0;
  sum2 = 0;
  for (i = 0; i < numPoints; i++)
    {
      fscanf (inFile1, "%lg %lg %lg %lg %lg", &d1,&d2,&d3,&d4,&d5);
      sum1 += pow(d1,2)+pow(d2,2)+pow(d3,2)+pow(d4,2)+pow(d5,2);
      sum2 += d1+d2+d3+d4+d5;
    }
  for (i = 0; i < extraPoints; i++)
    {
      fscanf (inFile1, "%lg", &d1);
      sum1 += pow(d1,2);
      sum2 += d1;
    }
  printf ("%15.8E\n",i1*i2*i3*sum1/pow(sum2,2));
  //printf ("%15.8E\n",sum2/i1/i2/i3);
  
  fclose (inFile1);
  return 0;
}
