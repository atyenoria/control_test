//
// hrandom++.h
//

#ifndef __HRANDOM_H__
#define __HRANDOM_H__

#define HRANDIB1 1
#define HRANDIB3 4
#define HRANDIB31 1073741824

#define HRANDMASK (HRANDIB3)

static unsigned long rand_param_pn=1699223980;


void init_hrand();         /* 乱数初期化 */
void hrand_seed(long t);         /* 乱数初期化 */

long hrandom();     /* random()と同じ */

double hrand_01m(); /* [0,1) の間で乱数発生 */

double hrand_01();  /* [0,1] の間で乱数発生 */

double hrand_11();  /* [-1,1) の間で乱数発生 */

int hrand_zp();            /* 0 と 1 の乱数発生 */

int hrand_pm();     /* -1 と 1 の乱数発生 */

double hrand(double Min, double Max);

/* ガウス分布によるランダム変数発生 */
double hrand_gauss(double sigma);

#endif // __HRANDOM_H__
