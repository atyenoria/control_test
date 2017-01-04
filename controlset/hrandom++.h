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


void init_hrand();         /* �������� */
void hrand_seed(long t);         /* �������� */

long hrandom();     /* random()��Ʊ�� */

double hrand_01m(); /* [0,1) �δ֤����ȯ�� */

double hrand_01();  /* [0,1] �δ֤����ȯ�� */

double hrand_11();  /* [-1,1) �δ֤����ȯ�� */

int hrand_zp();            /* 0 �� 1 �����ȯ�� */

int hrand_pm();     /* -1 �� 1 �����ȯ�� */

double hrand(double Min, double Max);

/* ������ʬ�ۤˤ��������ѿ�ȯ�� */
double hrand_gauss(double sigma);

#endif // __HRANDOM_H__
