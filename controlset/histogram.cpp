/*
 * histo *
 *  Created on: 2014/04/22
 *      Author: root
 */



#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <math.h>
#include <limits.h>
#include <typeinfo>
#include <fstream>
#include <ctime>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include "hrandom++.h"
#include "DTCS.h"
#include "PendulumR.h"
#include "LQ.h"
#define CONTROLE_MODE 1 // 0 for omission angel, 1 for omission velocity,2 for normal control
#define ERROR_MODE 1//0 for packet error rate ,1 for bit error rate
#define RMSE_MODE 1//0 for difference of the optimal value,		1 for difference of the reference



using namespace std;

// template function for calculating static variety
template<typename TYPE,std::size_t SIZE>std::size_t array_length(const TYPE(&)[SIZE]){return SIZE;}



int main(int argc, char **argv) {




	/***************defualt setting***************/
  int i, k, l,f;
  const int seed = 1;
  const int Loop = 3000;
  const double Hz = 100;
  const double Ts = 1.0/Hz;
  const int Time_Range = 60;
  stringstream ss;




  const double Arm_Angle_Period = 20;
  const double Pendulum_Angle_Range = M_PI/3;





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
        xref[k][ARM_ANGLE] = 0;
    } else {
    	 xref[k][ARM_ANGLE] = M_PI/4;
    }
  }



  dVector<4> xest[Sample_Range+2];
  double u[Sample_Range+2];
  double u_hat[Sample_Range+2];
  dVector<4> x[Sample_Range+2];
  dVector<4> x_hat[Sample_Range+2];
  dVector<4> x_controller[Sample_Range+2];
  int x_controller_index;
  dVector<4> w;
  dVector<4> xideal[Sample_Range+2];
  xideal[0] = 0.0;
  const dVector<4> K = dlqr(A,B,Q,R);
  for (k=0; k<=Sample_Range; k++) {
    xideal[k+1] = A*xideal[k] + B*(K*(xref[k] - xideal[k]));
  }

  cout<<"# Pendulum:"<<PENDULUM_MODEL<<endl;
  printf("#Q=diag:");printf("%2.lf",diag(Q));
  cout<<"# R =diag:"<<R<<endl;
  cout<<"# sampling rate:"<<Hz<<endl;
  cout<<"# simulation time:"<<Time_Range<<endl;
  cout<<"# ARM_ANGLE_PERIOD:"<<Arm_Angle_Period<<endl;
  cout<<"# Loop cout:"<<Loop<<endl;
  cout<<"# ERROR_MODE= "<<ERROR_MODE<<"//0 for packet error rate ,1 for bit error rate "<<endl;
  cout<<"# CONTROLE_MODE= "<<CONTROLE_MODE<<"/ 0 for omission angel, 1 for omission velocity,2 for normal control "<<endl;
  cout<<"# RMSE_MODE= "<<ERROR_MODE<<"//0 for difference of the optimal value,1 for difference of the reference\n\n\n"<<endl;

//designate disturbance
double sigma_w[]={1e-2,1e-3,1e-4};
int z=array_length(sigma_w);




//designate packet error rate or bit error rate

	  double p[]={
			  0.0000001,			//0 indicator
			  0.000001,		//1
			  0.00001,		//2
			  0.0001,		//3
			  0.001,		//4
			  0.01,			//5
			  0.05,			//6
			  0.1,			//7
			  0.2,			//8
	  	  	  0.3};			//9
   int p_packet_length=array_length(p);
//ERROR_MODE_0










	double p_bit[]={//bit error rate
			0,//0
			1.9073486328125E-6,//1
			1.1444091796875E-5,//2
			1.9073486328125E-5,//3
			2.47955322265625E-5,//4
			7.43865966796875E-5,//5
			0.0001220703125,//6
			0.000171661376953125,
			0.0002460479736328125,
			0.0003681182861328125,
			0.0007419586181640625//7
	};
	int p_bit_length=array_length(p_bit);
	double pbit_per_forward[p_bit_length];
	double pbit_per_backward_estimate[p_bit_length];
	double pbit_per_backward_not_estimate[p_bit_length];
//calculate forward channel's packet error rate

	for(k=0;k<p_bit_length;k++){
		pbit_per_forward[k]=1-pow((1-p_bit[k]),48);
		pbit_per_backward_estimate[k]=1-pow((1-p_bit[k]),72);
		pbit_per_backward_not_estimate[k]=1-pow((1-p_bit[k]),104);


	}


	cout<<"forward error rate\n";
	for(k=0;k<p_bit_length;k++){
		cout<<pbit_per_forward[k]<<"\n";
	}

	cout<<"backward error rate(estimate)\n";
	for(k=0;k<p_bit_length;k++){
			cout<<pbit_per_backward_estimate[k]<<"\n";
		}

	cout<<"backward error rate(not estimate)\n";
	for(k=0;k<p_bit_length;k++){
				cout<<pbit_per_backward_not_estimate[k]<<"\n";
			}







//ERROR_MODE_1




/* for saving
 0:packet error rate
 1:fall down rate
 2:arm angel RMSE
 3:pendulum angel RMSE
*/

	int count_error_index;

	switch(ERROR_MODE){

	case 0:
		count_error_index=p_packet_length;
		break;
	case 1:
		count_error_index=p_bit_length;
		break;

	}

	double error[count_error_index][5];






/**********starting point**********/


for(int number_disturbance=0;number_disturbance<z;number_disturbance++){




  for(int f=0;f<count_error_index;f++){
  int loop=0;
  int count_falldown=0;
  int count_reduction_number=0;
  dVector<4> count_mse = 0.0;
  const dVector<4> K_tcp = dlqr_tcp(A,B,Q,R,p[f]);
  const dVector<4> K_udp = dlqr_udp(A,B,Q,R,p[f],p[f]);
  int count_estimation=0;


  do{
    init_hrand();

    u[0] = 0.0;
    x[0] = xest[0] =x_controller[0]=0.0;

    int noerror_U=1;
    int noerror_X=1;
    int estimate=0;
    int falldown=0;
    int falldown_time=Sample_Range+1;
    int x_controller_index=0;

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


/********************start of forward channel***************************/


					// UDP-case
				   if (noerror_X) {
					  xest[k] = x_hat[k-1];
					  u[k] = K*(xref[k] - xest[k]);
					} else {
					  xest[k] = 0.0;
					  u[k] = 0.0;
					}
	  	  	  	  }


      //u[k] = K*(xref[k] - xest[k]);
      //u[k] = K_tcp*(xref[k] - xest[k]);

switch(ERROR_MODE){
case 0:
	noerror_U = (hrand_01() >= p[f]);
	break;
case 1:
	noerror_U=(hrand_01() >=pbit_per_forward[f]);
	break;

}



      if (noerror_U) {
	u_hat[k] = u[k];
      } else {
	u_hat[k] = 0;
      }
 /********************end of forward channel***************************/




      for (i=0; i<4; i++) w[i] = hrand_gauss(sigma_w[number_disturbance]);

      x[k+1] = A*x[k] + B*u_hat[k] + w;



      if(estimate>=5){estimate=0;}


      switch(ERROR_MODE){
      case 0:
    	   noerror_X = (hrand_01() >= p[f]);
      	break;
      case 1:
    	  if(estimate>=1){
      	noerror_U=(hrand_01() >=pbit_per_backward_estimate[f]);}
    	  else{
        noerror_U=(hrand_01() >= pbit_per_backward_not_estimate[f]);
      }
    	  break;
      }





      if (noerror_X) {
    	  switch(CONTROLE_MODE){

    	  case 0://omission angle

    		  	  if(estimate>=1){
    		  	  x_controller_index=k;
    		  	  x_controller[k][ARM_ANGLE]=x_controller[k-1][ARM_ANGLE]+Ts*x[k][ARM_VELOCITY];
    		  	  x_controller[k][ARM_VELOCITY]=x[k][ARM_VELOCITY];
    		  	  x_controller[k][PENDULUM_ANGLE]=x_controller[k-1][PENDULUM_ANGLE]+Ts*x[k][PENDULUM_VELOCITY];
    		  	  x_controller[k][PENDULUM_VELOCITY]=x[k][PENDULUM_VELOCITY];
    		  	  x_hat[k]=x_controller[k];
    		  	  count_reduction_number++;
    		  	  estimate++;
    		  	  }else{
    		  		  x_controller_index=k;
    		  		  x_controller[k]=x[k];
    		  		  x_hat[k]=x_controller[k];

    		  		  estimate++;
    		  	  }
    		  	  	  	  	  	  break;
    	  case 1:

    	  	  if(estimate>=1){
        		  	  x_controller_index=k;
        		  	  x_controller[k][ARM_ANGLE]=x[k][ARM_ANGLE];
        		  	  x_controller[k][ARM_VELOCITY]=(x[k][ARM_ANGLE]-x[k-1][ARM_ANGLE])/0.01;
        		  	  x_controller[k][PENDULUM_ANGLE]=x[k][PENDULUM_ANGLE];
        		  	  x_controller[k][PENDULUM_VELOCITY]=(x[k][PENDULUM_ANGLE]-x[k-1][PENDULUM_ANGLE])/0.01;
        		  	  x_hat[k]=x_controller[k];
        		  	  count_reduction_number++;
        		  	  estimate++;
        		  	  }else{
        		  		  x_controller_index=k;
        		  		  x_controller[k]=x[k];
        		  		  x_hat[k]=x_controller[k];
        		  		  estimate++;

        		  	  }


    		  	  	  	  	  	  break;


    	  case 2:
    		  	  	  	  x_controller_index=k;
    		          	  x_controller[k]=x[k];
    		          	  x_hat[k]=x_controller[k];


    		          	  break;


    	  }
      } else {
	x_hat[k] = 0.0;
	estimate=0;
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

    	switch(RMSE_MODE){

    	case 0:
					  for (k=1; k<=Sample_Range; k++) {
					count_mse += square(x[k] - xideal[k]);
					  }
					  	  break;
    	case 1:
					  for (k=1; k<=Sample_Range; k++) {
								count_mse += square(x[k] - xref[k]);
								  }
					      break;



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

  double reductionrate=((double)count_reduction_number*(double)16)/((double)loop*(double)Sample_Range*(double)144);
count_mse/= ((double)loop-(double)count_falldown)*(double)Sample_Range;

switch(ERROR_MODE){
case 0:
cout<<"#packet error rate"<<"\t"<<p[f]<<endl;
                                    break;
case 1:
cout<<"#packet error rate(bit error mode)"<<"\t"<<pbit_per_forward[f]<<endl;
break;



}


cout<<"#disturbance"<<sigma_w[number_disturbance]<<"\n"<<endl;
cout<<"#falldown rate\t"<<count_falldown<<"/"<<loop<<"%"<<endl;
cout<<"#RMSE\t"<<"pendulum_angle\t"<<"arm_angel\t"<<"pendulum_velocity\t"<<"arm_velocity\t"<<"\n"<<endl;
cout<<"#rad\t"<<sqrt(count_mse[PENDULUM_ANGLE])<<"\t"<<sqrt(count_mse[ARM_ANGLE])<<
		"\t"<<sqrt(count_mse[PENDULUM_VELOCITY])<<"\t"<<sqrt(count_mse[ARM_VELOCITY])<<endl;
cout<<"#deg\t"<<180*sqrt(count_mse[PENDULUM_ANGLE])/M_PI<<"\t"<<180*sqrt(count_mse[ARM_ANGLE])/M_PI
		<<"\t"<<sqrt(count_mse[PENDULUM_VELOCITY])*180/M_PI<<"\t"<<sqrt(count_mse[ARM_VELOCITY])*180/M_PI<<"\n\n\n"<<endl;
cout<<"#redction rate\t"<<setprecision(10)<<count_reduction_number<<"\t"<<reductionrate<<"%\n\n\n\n"<<endl;






/**************write data*********/
switch(ERROR_MODE){
case 0:
	error[f][0]	=	p[f];
	break;
case 1:
	error[f][0]	=pbit_per_forward[f];

break;



}

error[f][1]	=	double(count_falldown/loop);
error[f][2]	=	180*sqrt(count_mse[ARM_ANGLE])/M_PI;
error[f][3]	=	180*sqrt(count_mse[PENDULUM_ANGLE])/M_PI;
error[f][4] =   reductionrate;
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

ss<<"omission_reference_biterror_"<<sigma_w[number_disturbance];
	std::ofstream ofs(ss.str().c_str());
	ofs<<"#packet error rate\t"<<"#fall down rate\t"<<"#arm angel EMSE\t"<<"#pendulum angel RMSE\t"<<"#all degree"<<endl;
for(i=0;i<count_error_index;i++){
for(k=0;k<4;k++){
	ofs<<error[i][k]<<",";
}
ofs<<endl;
}
ofs.close();


cout<<"\n\n\n\n********************************************************************************************\n\n\n\n\n\n"<<endl;

}
  /*****up to here changing disturbance****/




}

