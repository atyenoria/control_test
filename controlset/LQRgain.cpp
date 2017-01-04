
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
#define CONTROLE_MODE 2	 // 0 for omission angel, 1 for omission velocity,2 for normal control
#define ERROR_MODE 0//0 for packet error rate ,1 for bit error rate
#define RMSE_MODE 1//0 for difference of the optimal value,		1 for difference of the reference

using namespace std;

// template function for calculating static variety
template<typename TYPE, std::size_t SIZE> std::size_t array_length(
		const TYPE (&)[SIZE]) {
	return SIZE;
}

int main(int argc, char **argv) {

	/***************defualt setting***************/
	int i, k, l, f;
	const int seed = 1;
	const int Loop = 1000;
	const double Hz = 100;
	const double Ts = 1.0 / Hz;
	const int Time_Range = 60;
	stringstream ss;

	const double Arm_Angle_Period = 20;
	const double Pendulum_Angle_Range = M_PI / 2;

	dMatrix<4, 4> Ac;
	dVector<4> Bc;
	pendulum(Ac, Bc);
	cout<<"\n continuous:A"<<endl;
	printf("%lf",Ac);
	cout<<"\n continuous:B"<<endl;
	printf("%lf",Bc);
	dMatrix<4, 4> Q = (double[4][4] ) {
		{ 1, 0, 0, 0 },
		{ 0, 10, 0, 0 },
		{ 0, 0,0.5, 0 },
		{ 0, 0, 0, 0 } };
	double R = 5;
	const dIdentityMatrix<4, 4> I;
	dMatrix<4, 4> A = Ac;
	dVector<4> B = Bc;
	dscr(A, B, Ts);
	cout<<"\n discrete:A"<<endl;
	printf("%lf",A);
	cout<<"\n discrete:B"<<endl;
	printf("%lf",B);
	cout<<"\n Q:"<<endl;
	printf("%lf",Q);
	cout<<"\n R:"<<endl;
	printf("%lf",R);
	const dVector<4> K = dlqr(A, B, Q, R);
	cout<<"\n K:振り子角度　アーム角度　振り子角速度　アーム角速度"<<endl;
	printf("%lf,K",K);
	cout<<"thetaで表記"<<endl;
	printf("%lf",K/180*M_PI);
	cout<<"asdfa<<"<<endl;
	cout<<endl;

}
