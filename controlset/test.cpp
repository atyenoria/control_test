/*
 * calculate.cpp
 *
 *  Created on: May 21, 2014
 *      Author: atyenoria
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
using namespace std;

//パラメーター一覧
 const double PI=3.14159265;

 const double mp=0.016; //振り子の質量[kg]
 const double lp=0.213;  //振り子の長さ[m]
 const double Lcp=0.085; //振り子の重心までの距離[m]
 const double Jp=1;	//振り子の重心周りの慣性モーメント[kg*m^2]
 const double Dp=1;   //振り子の回転軸の粘性摩擦[kg*m^2/s]


 const double Jr=0.0048; //アームの回転軸周りの慣性モーメント[kg*m^2]
 const double Lr=0.206;   //アームの長さ[m]

 const double g=9.797;   //重力加速度[m/s^2]
 const double R=8.3;    //モータトルク直流抵抗[Ω]
 const double Km=0.023; //モータトルク定数[Nm/A,Vs/rad]
 // トルク:τ(t) = amVm(t) − bmθ ̇(t)
 //am=Kg*Km/R 	bm=Kg^2*Km^2/R
 const double Kg=7.5;   //ギア比 120:16
 const double Dr=1;   //アームの回転軸の粘性摩擦[kg*m^2/s]

//アームの角度　θ(t),
//振り子の角度　φ(t)
//アームのの角速度θ'(t)
//振り子の角速度φ'(t)


int main(int argc, char **argv) {

   //JpとDpを計算する


	const double am=Kg*Km/R;
	cout<<"am:"<<am<<endl;
	const double bm=Kg*Kg*Km*Km/R;
	cout<<"bm:"<<bm<<endl;
	const double Jt=Jp*mp*Lr*Lr+Jr*Jp+Jr*mp*Lcp*Lcp;




	//A行列の係数を決定する
	const double a=mp*mp*Lcp*Lcp*Lr*g/Jt;
	const double b=-(Jp+mp*Lcp*Lcp)*(Dr+bm)/Jt;
	const double c=-mp*Lcp*Lr*Dp/Jt;
	const double d=mp*Lcp*g*(Jr+mp*Lr*Lr)/Jt;
	const double e=mp*Lcp*Lr*(Dr+bm)/Jt;
	const double f=-(Jr+mp*Lr*Lr)*Dp/Jt;

	//B行列の係数を決定する


	const double gg=am*(Jp+mp*Lcp*Lcp)/Jt;
	const double hh=am*mp*Lcp*Lr/Jt;


	  dMatrix<4,4> A = (double[4][4]){{0,0,1,0},
									  {0,0,0,1},
									  {0,a,b,c},
									  {0,d,e,f}};

cout<<A[4][4]<<endl;
	 dVector<4> B =(double[4]){0, 0, gg, hh};




	//こっからゲインを計算する
	 	const int seed = 1;
	 	const double Hz = 100;
	 	const double Ts = 1.0 / Hz;


	 	dMatrix<4, 4> At=A;
	 	dVector<4> Bt=B;

	 	dMatrix<4, 4> Q = (double[4][4] ) { { 10, 0, 0, 0 },
											{ 0, 2, 0, 0 },
											{ 0, 0,0.5, 0 },
											{ 0, 0, 0, 0 } };

	 	double R = 0.05;
	 	cout<<"A行列:"<<endl;
		printf("%lf",At);
		cout<<"B行列："<<endl;
		printf("%lf",Bt);
		cout<<endl;
	 	dscr(At, Bt, Ts);

		cout<<"discrete:A"<<endl;
		printf("%lf",At);
		cout<<"discrete:B"<<endl;
		printf("%lf",Bt);

	 	const dVector<4> K = dlqr(At,Bt,Q,R);
		cout<<"K:"<<endl;
		printf("%lf",K);


}



