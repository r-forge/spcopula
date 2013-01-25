/*
** tools.c - C code of the package CDRVine  
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

//////////////////////////////////////////////////////////////////////////////
// Print error text and return to R
//////////////////////////////////////////////////////////////////////////////
void printError(char *text, char filename[200])
{
  Rprintf(text);
  Rprintf(": %s ", filename);
  Rprintf(" !!!\n");
}

//////////////////////////////////////////////////////////////////////////////
// Compute gamma division in a more stable way than gamma(x1)/gamma(x2)
// Input:
// x1 - Divisor
// x2 - Denominator
//-------------------
// a1 - modulus of a
// a2 - integer part of a
// b1 - modulus of b
// b2 - integer part of b
//-------------------
// We are after gamma(x1)/gamma(x2) and this computation will hopefully make it more numerically 
// stable. gamma(x1)/gamma(x2) will sometimes be INF/INF.
//////////////////////////////////////////////////////////////////////////////
double StableGammaDivision(double x1, double x2)
{
  int i;
  double a1, a2, b1, b2, sum=1.0;
  a1 = fmod(MAX(x1,x2),1.0);
  a2 = MAX(x1,x2)-a1;
  b1 = fmod(MIN(x1,x2),1.0);
  b2 = MIN(x1,x2)-b1;
  if(a1==0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=b2 ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
  }
  else if(a1>0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=(int)b2 ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1);
  }
  else if(a1==0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum /= gammafn(b1);
  }
  else if(a1>0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1)/gammafn(b1);
  }
  if(x2 > x1) sum = 1.0/sum;
  return sum;
}