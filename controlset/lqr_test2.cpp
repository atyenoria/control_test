
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <limits.h>
#include <typeinfo>
#include <fstream>
#include <ctime>
#include <sstream>
#include <cstring>
#include "hrandom++.h"
#include "DTCS.h"
#include "PendulumR.h"
#include "LQ.h"


using namespace std;

// template function for calculating static variety
template<typename TYPE,std::size_t SIZE>std::size_t array_length(const TYPE(&)[SIZE]){return SIZE;}



int main(int argc, char **argv) {




	/***************defualt setting***************/
  int i, k, l,f;
  const int seed = 1;
  const int Loop = 200;
  const double Hz = 100;
  const double Ts = 1.0/Hz;
  const int Time_Range = 100;
  stringstream ss;




  const double Arm_Angle_Period = 10;
  const double Pendulum_Angle_Range = M_PI/2;





  dMatrix<4,4> Ac; dVector<4> Bc;
  pendulum(Ac, Bc);
  dMatrix<4,4> Q = (double[4][4]){{10, 0,  0,  0},
				  {0,  2,  0,  0},
				  {0,  0, 0.5, 0},
				  {0,  0,  0,  0}};
  double R = 0.05;
  const dIdentityMatrix<4,4> I;
  dMatrix<4,4> A = Ac; dVector<4> B = Bc;
  dscr(A, B, Ts);
  /************ end of default setting ***************/








  const int Sample_Range = (int)(Time_Range*Hz);
  dVector<4> xref[Sample_Range+2];

  for (k=0; k<=Sample_Range; k++) {
    xref[k] = 0.0;
    if ( fmod(k*Ts,Arm_Angle_Period) < Arm_Angle_Period/2) {
      xref[k][ARM_ANGLE] = M_PI/2;
    } else {
      xref[k][ARM_ANGLE] = 0;
    }
  }



  dVector<4> xest[Sample_Range+2];
  double u[Sample_Range+2];
  double u_hat[Sample_Range+2];
  dVector<4> x[Sample_Range+2];
  dVector<4> x_hat[Sample_Range+2];
  dVector<4> w;
  dVector<4> xideal[Sample_Range+2];
  xideal[0] = 0.0;
  const dVector<4> K = dlqr(A,B,Q,R);
  for (k=0; k<=Sample_Range; k++) {
    xideal[k+1] = A*xideal[k] + B*(K*(xref[k] - xideal[k]));
  }

  cout<<"# Pendulum:"<<PENDULUM_MODEL<<endl;
  printf("Q=diag:");printf("%2.lf",diag(Q));
  cout<<"# R =diag:"<<R<<endl;
  cout<<"# sampling rate:"<<Hz<<endl;
  cout<<"# simulation time:"<<Time_Range<<endl;
  cout<<"# Loop cout:"<<Loop<<"\n\n\n\n\n\n\n"<<endl;




//designate disturbance
double sigma_w[]={1e-2,1e-4,1e-6};
int z=array_length(sigma_w);




//designate packet error rate
  double p[]={
		  0.0000001,			//0 indicator
		  0.0001,		//1
		  0.0005,		//2
		  0.2,		//3
		  0.3,		//4
		  0.4,			//5
		  0.5,			//6
		  0.6,			//7
		  0.7,			//8
  	  	  0.8};			//9
int q=array_length(p);
double error[q][4];

/* for saving
 0:packet error rate
 1:fall down rate
 2:arm angel RMSE
 3:pendulum angel RMSE
*/








/**********starting point**********/


for(int number_disturbance=0;number_disturbance<z;number_disturbance++){



  for(int f=0;f<q;f++){
  int loop=0;
  int count_falldown=0;
  dVector<4> count_mse = 0.0;
  const dVector<4> K_tcp = dlqr_tcp(A,B,Q,R,p[f]);
  const dVector<4> K_udp = dlqr_udp(A,B,Q,R,p[f],p[f]);


  do{
    init_hrand();

    u[0] = 0.0;
    x[0] = xest[0] = 0.0;

    int noerror_U=1;
    int noerror_X=1;
    int falldown=0;
    int falldown_time=Sample_Range+1;

    k=0;
    do {
    			if (k>0) {

								/*
								// TCP-case
								if (noerror_X) {
								  xest[k] = A*x_hat[k-1] + noerror_U*B*u[k-1];
								} else {
								  xest[k] = A*xest[k-1] + noerror_U*B*u[k-1];
								}
								*/

								// UDP-case
							   if (noerror_X) {
								  xest[k] = x_hat[k-1];
								} else {
								  xest[k] = 0.0;
								}
      	  	  	  	  	  }


/********************start of forward channel***************************/

      //u[k] = K*(xref[k] - xest[k]);
      //u[k] = K_tcp*(xref[k] - xest[k]);
      u[k] = K*(xref[k] - xest[k]);


      noerror_U = (hrand_01() >= p[f]);

      if (noerror_U) {
	u_hat[k] = u[k];
      } else {
	u_hat[k] = 0;
      }
 /********************end of forward channel***************************/




      for (i=0; i<4; i++) w[i] = hrand_gauss(sigma_w[number_disturbance]);

      x[k+1] = A*x[k] + B*u_hat[k] + w;

      noerror_X = (hrand_01() >= p[f]);

      if (noerror_X) {
	x_hat[k] = x[k];
      } else {
	x_hat[k] = 0.0;
      }

      falldown =
    		(fabs(x[k+1][PENDULUM_ANGLE]) >= Pendulum_Angle_Range);
      if (falldown) {
    	  break;
      }


      k++;
    } while (k <= Sample_Range+1);

    if (falldown) {
      count_falldown++;
    } else {
      for (k=1; k<=Sample_Range; k++) {
	count_mse += square(x[k] - xideal[k]);
      }

    }

/***********display details regarding specific packet error rate and the number of loop********/
/*
    printf("# t\t Pen_ang\t Arm_ang\t Pen_vel\t Arm_vel\t Voltage\t Arm_ref\n");
    for (k=0; k<falldown_time; k++) {
      printf("%.3lf\t % e\t % e\t % e\t % e\t % e\t % e\n",
    		  k/Hz,
	     x[k][PENDULUM_ANGLE],
	     x[k][ARM_ANGLE],
	     x[k][PENDULUM_VELOCITY],
	     x[k][ARM_VELOCITY],
	     u_hat[k],
	     xref[k][ARM_ANGLE]);
	      }
*/

    loop++;
  } while (loop<Loop);


count_mse/= (loop-count_falldown)*Sample_Range;
cout<<"#packet error rate"<<"\t"<<p[f]<<"\n"<<endl;
cout<<"#falldown rate\t"<<count_falldown<<"/"<<loop<<"%"<<endl;
cout<<"#RMSE\t"<<"pendulum_angle\t"<<"arm_angel\t"<<"pendulum_velocity\t"<<"arm_velocity\t"<<"\n"<<endl;
cout<<"#rad\t"<<sqrt(count_mse[PENDULUM_ANGLE])<<"\t"<<sqrt(count_mse[ARM_ANGLE])<<
		"\t"<<sqrt(count_mse[PENDULUM_VELOCITY])<<"\t"<<sqrt(count_mse[ARM_VELOCITY])<<endl;
cout<<"#deg\t"<<180*sqrt(count_mse[PENDULUM_ANGLE])/M_PI<<"\t"<<180*sqrt(count_mse[ARM_ANGLE])/M_PI
		<<"\t"<<sqrt(count_mse[PENDULUM_VELOCITY])*180/M_PI<<"\t"<<sqrt(count_mse[ARM_VELOCITY])*180/M_PI<<"\n\n\n\n\n\n\n"<<endl;




/**************write data*********/
error[f][0]	=	p[f];
error[f][1]	=	double(count_falldown/loop);
error[f][2]	=	180*sqrt(count_mse[ARM_ANGLE])/M_PI;
error[f][3]	=	180*sqrt(count_mse[PENDULUM_ANGLE])/M_PI;
     /*    //designate packet error rate
           double p[]={
         		  0,			//0 indicator
         		  0.0001,		//1
         		  0.0005,		//2
         		  0.001,		//3
         		  0.005,		//4
         		  0.01,			//5
         		  0.05,			//6
         		  0.1,			//7
         		  0.001,			//8
           	  	  0.0001};			//9
         int q=array_length(p);
         double error[q][4];

          for saving
          0:packet error rate
          1:fall down rate
          2:arm angel RMSE
          3:pendulum angel RMSE
         */
/*********end of write data**********/


  }

/**** up to here changing error rate ****/


ss.clear();
ss.str("");

ss<<number_disturbance<<"datasdfaa"<<sigma_w[number_disturbance];
	std::ofstream ofs(ss.str().c_str());
	ofs<<"#packet error rate\t"<<"#fall down rate\t"<<"#arm angel EMSE\t"<<"#pendulum angel RMSE\t"<<"#all degree"<<endl;
for(i=0;i<9;i++){
for(k=0;k<4;k++){
	ofs<<error[i][k]<<",";
}
ofs<<endl;
}
ofs.close();


}
  /*****up to here changing disturbance****/




}

