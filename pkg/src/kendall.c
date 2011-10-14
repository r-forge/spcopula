/*
Vn must be sorted
v the argument
n number of observations
z the value of the function at v
*/

void emp_kendallDist(double *Vn, double *v, int *n, double *z)
{
  int i;
  double s=0.0;
    
  for (i=0;i<*n;i++)
    {
      if ((double) Vn[i] <= (double) v[0]) s++;
    }
  *z = s / (double) *n;
}


/*void empCop(double *un, double *vn, int *n, int *d, double *p)
{
  int i;
  double s=0.0;
   
  for (i=0;i<*n;i++)
    {
      if ((double) un[i] <= (double) v[0]) s++;
    }
  *p = 
}

 res <- NULL
    for(i in 1:nrow(u)) {
      bool <- t(t(data) <= u[i,])
      for (i in 2:ncol(data)) bool[,1] <- bool[,1] * bool[,i]
      res <- c(res,sum(bool[,1]))
*/