/*
** memoryhandling.c - C code of the package CDRVine  
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/

#include "vine.h"

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates a double matrix.
//  Input: Dimension of the matrix to be created//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
double **create_matrix(int rows, int columns)
{
  double **a;
  int i=0;
  a = (double**) Calloc(rows, double*);
  for(i=0;i<rows;i++) a[i] = (double*) Calloc(columns,double);
  return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that a double matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_matrix(double **a, int rows)
{
  int i=0;
  for(i=0;i<rows;i++) Free(a[i]);
  Free(a);
}

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates an int matrix.
//  Input: Dimension of the matrix to be created.
//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
int **create_intmatrix(int rows, int columns)
{
  int **a;
  int i=0;
  a = (int**) Calloc(rows,int*);
  for(i=0;i<rows;i++) a[i] = (int*) Calloc(columns,int);
  return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that an int matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_intmatrix(int **a, int rows)
{
  int i=0;
  for(i=0;i<rows;i++) Free(a[i]);
  Free(a);
}