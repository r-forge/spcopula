/*
** cdvine.c - C code of the package CDRVine  
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


#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4


//////////////////////////////////////////////////////////////
// Function to simulate from a pair-copula construction (vine)
// Input:
// n         sample size
// d         dimension (>= 2)
// type      vine type (1=Canonical vine, 2=D-vine)
// family    copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1)
// par       parameter values (at least d*(d-1)/2 parameters)
////////////////////////////////////////////////////////////////

void pcc(int* n, int* d, int* family, int* type, double* par, double* nu, double* out)
{
  int i, j, in=1, k, **fam;
  double *w, **v, t, **theta, **x, **ny;

  GetRNGstate();
  //Allocate memory:
  w = Calloc((*d+1),double);

  v = create_matrix(*d+1,2*(*d)-1);
  theta = create_matrix(*d,*d);
  x = create_matrix(*n+1,*d+1);
  ny = create_matrix(*d,*d);
  fam = create_intmatrix(*d,*d);
  //Initialize dependency parameters
  k = 0;
  for(i=1;i<=*d-1;i++)
  {
    for(j=1;j<=*d-i;j++)
    {
      fam[i][j] = family[k];
      ny[i][j] = nu[k];
      theta[i][j] = par[k];
      k ++;
    }
  }
  //Simulate:
  if(*type==1) //Canonical vine
  {
    for(j=1;j<=*n;j++)
    {
      for(i=1;i<=*d;i++) w[i] = runif(0,1);
      x[j][1] = w[1];
      for(i=2;i<=*d;i++)
      {
        t = w[i];
        for(k=i-1;k>=1;k--)
        {
		  Hinv1(&fam[k][i-k],&in, &t,&w[k],&theta[k][i-k],&ny[k][i-k],&t);
        }
        x[j][i] = t;
      }
    }
  }
  else if(*type==2) //D-vine
  {
    for(j=1;j<=*n;j++)
    {
      for(i=1;i<=*d;i++) { w[i] = runif(0,1);}
      v[1][1] = w[1];
      v[2][1] = w[2];
      Hinv1(&fam[1][1],&in,&w[2],&v[1][1],&theta[1][1],&ny[1][1],&v[2][1]);
      Hfunc2(&fam[1][1],&in, &v[1][1],&v[2][1],&theta[1][1],&ny[1][1],&v[2][2]);
      for(i=3;i<=*d;i++)
      {
        v[i][1] = w[i];
	
        for(k=i-1;k>=2;k--) { 
	  Hinv1(&fam[k][i-k],&in, &v[i][1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][1]);
	}
        Hinv1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][1]);
        // Compute conditional cdf's needed in next step:
        if(i<*d)
        {
          Hfunc2(&fam[1][i-1],&in, &v[i-1][1],&v[i][1],&theta[1][i-1],&ny[1][i-1],&v[i][2]);
          Hfunc1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][3]);
          if(i>3)
          {
            for(k=2;k<=(i-2);k++)
            {
              Hfunc2(&fam[k][i-k],&in, &v[i-1][2*k-2],&v[i][2*k-1],&theta[k][i-k],&ny[k][i-k],&v[i][2*k]);
              Hfunc1(&fam[k][i-k],&in, &v[i][2*k-1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][2*k+1]);
            }
          }
          Hfunc2(&fam[i-1][1],&in, &v[i-1][2*i-4],&v[i][2*i-3],&theta[i-1][1],&ny[i-1][1],&v[i][2*i-2]);
        }
      }
      for(i=1;i<=*d;i++) x[j][i] = v[i][1];
    }
  }
  //Write to output vector:
  k = 0;
  for(i=1;i<=*d;i++)
  {
    for(j=1;j<=*n;j++)
    {
      out[k] = x[j][i];
      k ++;
    }
  }
  PutRNGstate();
  //Free memory:
  Free(w); free_matrix(v,*d+1); free_matrix(theta,*d); free_matrix(ny,*d); free_intmatrix(fam,*d); free_matrix(x,*n+1);
}