#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "hrandom++.h"

void init_hrand()         /* 乱数初期化 */
{
  time_t hrhr;
  long t;
  t = (long)time(&hrhr);
  srand48(t);

  rand_param_pn = lrand48();
}

void hrand_seed(long t)  /* 乱数初期化 */
{
  srand48(t);

  rand_param_pn = lrand48();
}

long hrandom()     /* random()と同じ */
{
  return lrand48();
}

double hrand_01m() /* [0,1) の間で乱数発生 */
{
  return drand48();
}

double hrand_01()  /* [0,1] の間で乱数発生 */
{
  return ((double)lrand48()/2147483647.0);
}

double hrand_11()  /* [-1,1) の間で乱数発生 */
{
  return (drand48()*2-1);
  //return ((double)mrand48()/2147483648.0);
}

int hrand_zp()            /* 0 と 1 の乱数発生 */
{
  if (rand_param_pn & HRANDIB31) {
    rand_param_pn = ((rand_param_pn ^ HRANDMASK) << 1) | HRANDIB1;
    return 1;
  } else {
    rand_param_pn <<= 1;
    return 0;
  }
  //return ((int)lrand48()&01);
}

int hrand_pm()     /* -1 と 1 の乱数発生 */
{
  return (hrand_zp()*2-1);
  //return ((int)lrand48()&01)*2-1);
}


double hrand(double Min, double Max)
{
  double Absval, Randval;
  
  Absval = fabs(Max-Min);
  Randval = (double)lrand48()/2147483647.0 * Absval;
  
  if (Max < Min)
    Randval += Max;
  else
    Randval += Min;
  
  return (Randval);
}

/* ガウス分布によるランダム変数発生 */
double hrand_gauss(double sigma)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if (iset==0){
    do{
      v1=hrand_11();
      v2=hrand_11();
      rsq=v1*v1+v2*v2;
    }while (rsq>=1.0 || rsq==0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    
    gset=v1*fac*sigma;
    iset=1;
    return (v2*fac*sigma);
  }
  else{
    iset=0;
    return (gset);
  }
}


/*
main()
{
  int i,b;
  double a,sum=0.0;

  init_hrand();
  for(i=0;i<1000;i++)
    {
      a=hrand_11();
      sum+=a;
      printf("%d\t%f\n",i,a);
    }
  printf("# ave %f\n",sum/1000.0);
}
*/

