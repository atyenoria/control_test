#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "hrandom++.h"
#include "DTCS.h"
#include "LQ.h"
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
#include "LQ.h"
using namespace std;
#define MODE 0 //0 for no error,  1 for error

#define PEN 0
#define PEN_VELOCITY 2
#define ARM 1
#define ARM_VELOCITY 3


//振り子の初期値の設定
#define PEN_INI 27 ;//振り子角度
#define PEN_VELOCITY_INI 1 ;//振り子角速度

#define ARM_INI 3996  ;//アーム角度
#define ARM_VELOCITY_INI 1 ;//アーム角速度

#define PACKET 1500
//0,	-51,	-1,	4195,	-1
//1,	1,	1,	3999,	1
//10,	15,	0,	4157,	0
//50,	-15,	1,	4051,	1
//100,	47,	-1,	4421,	-1
//500,	7,	1,	3622,	1
//1000,	18,	-2,	3925,	-2
//1500,	27,	1,	3996,	1
//2000,	-57,	5,	3627,	5




list<string> split(string str, string delim)
{

    list<string> result;
    int cutAt;
    while( (cutAt = str.find_first_of(delim)) != str.npos )
    {
        if(cutAt > 0)
        {
            result.push_back(str.substr(0, cutAt));
        }
        str = str.substr(cutAt + 1);
    }
    if(str.length() > 0)
    {
        result.push_back(str);
    }
return result;
}




int main(){

	std::stringstream ss;

    //vector<vector<int> > values;
    vector<double> value_arm;
    vector<double> value_pendulum;
    vector<double> value_time_real;
    vector<double> value_packet;
    vector<double> interval_index;
    vector<double> mvalue;//最大値と最小値を格納
    vector<double> mvalue_time;//最大値と最小値を取る時刻を保存
    vector<double> JRR;
    vector<double> DRR;


    dVector<4> x[10000];

    string str;
    const double PI=3.14159265;
    int p;		// デリミタの指定を行う
    int q=0;		//書き込みのためのオフセットを行う
    int f=0; //振り子、アームの角度配列の指標（時刻を表していることになる）
    int b=0; //収束値を跨いでいる点を探すための指標, その他
    int c=0; //最大値最小値を求める際の窓の移動に関しての指標
    int tmp_value;//最大値、最小値の一時格納場所
    int tmp_time; //最大値、最小値を取る時刻の一時格納場所
    int n;//大ループ指標
    int k;//小ループ指標
    const double enc_pendulum=2*M_PI/10000;
    const double enc_arm=2*M_PI/15000;
   string S;
   bool Z=false;// 書き出しのはじめを決める true=書き込む
   bool N=false;//振幅運動の最大値最小値を求めるために利用する。
   double pam[3]; //読み取り値の一時格納場所]




for(n=0;n<1;n++){



	ss.str("");
	ss<<n<<"_"<<PACKET<<".csv";
	cout<<ss.str()<<endl;
	string filename;
	filename=ss.str();
	ifstream file(filename);

								if(file.fail()){
									cerr << "failed." << endl;
									exit(0);
								}

int t=0;


    while(getline(file, str)){
    	p=0;

    	if(str.find(",")!=0){

        // splitを実行します。
        list<string> strList = split(str, ",");

        // イテレータを取得します。
        list<string>::iterator iter = strList.begin();

        // 侮ｦしてみます。
        while( iter != strList.end() ) // 最後まで
        {
        	S=*iter;


            pam[p]=atof(S.c_str());
            // １つ後方に進む
            p++;
            ++iter;
        }


        if(pam[0]==1){Z=true;
        	cout<<"start time :"<<t<<endl;
        }//振り子、アームパラメータ同定で設定しなければならない

        if(Z){



											value_packet.push_back(pam[1]);
											value_arm.push_back(pam[2]);
											value_time_real.push_back(pam[0]);
											f++;





               }

    	}
        t++;
    }

    cout<<"size_t :"<<value_arm.size()<<endl;


//以前のシステム同定法

const double m=0.016; //振り上げ棒の質量[kg]
const double l=0.085; //振り上げ棒の長さの1/2[m]
const double J=0.0050; //アームの回転軸周りの慣性モーメント[kgm^2]
const double r=0.206; //アームの長さ[m]
const double g=9.797; //重力加速度[m/s^2]


const double R=8.3; //モータトルク直流抵抗[Ω]
const double Km=0.023; //モータトルク定数[Nm/A,Vs/rad]
const double Kg=7.5; //ギア比 120:16



 //操作量がアームモータの入力電圧の場合

double a = -m*r*g/J;
double b = -(Kg*Km)*(Kg*Km)/(R*J);
double c = g*(J+m*r*r)/(l*J);
double d = r*(Kg*Km)*(Kg*Km)/(R*l*J);
double e = Kg*Km/(R*J);
double f = -Kg*Km*r/(R*l*J);




dMatrix<4,4> A = (double[4][4]){{0,0,1,0},
								  {0,0,0,1},
								  {c,0,0,d},
								  {a,0,0,b}};

cout<<A[4][4]<<endl;
dVector<4> B =(double[4]){0, 0, f, e};




//こっからゲインを計算する
 	const int seed = 1;
 	const double Hz = 100;
 	const double Ts = 1.0 / Hz;


 	dMatrix<4, 4> At=A;
 	dVector<4> Bt=B;

 	dMatrix<4, 4> Q = (double[4][4] ) { { 10, 0, 0, 0 },
										{ 0, 2, 0, 0 },
										{ 0, 0,0, 0 },
										{ 0, 0, 0, 0.5 } };

 	double Rw = 0.05;

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

 	const dVector<4> K = dlqr(At,Bt,Q,Rw);
 	cout<<endl<<endl;;
 	cout<<"　振り子角度 アーム角度　振り子速度 アーム速度"<<endl;
 	cout<<"radでのゲイン"<<endl;
	printf("%lf",K);
	cout<<"thetaでのゲイン"<<endl;
	printf("%lf",K/180*PI);







	//ここから制御開始、初期値を設定してシミュレーションを行う

	const int Sample_Range=1000000;//サンプル数×サンプリング時間

	dVector<4> xref[int(Sample_Range) + 2];

		for (k = 0; k <= Sample_Range; k++) {
			xref[k] = 0.0;
			xref[k][ARM] = 0.0;
		}


	int i;

	double u[Sample_Range + 2];
    dVector<4> xest[Sample_Range+2];
	//初期値の設定（状態情報の順番　振り子角度　アーム角度　振り子角速度　アーム角速度　の順番）


	//アーム角度のオフセット値を計算する

	double Amax=-20000000;
	double Amin=200000000;
	for(k=600;k<value_packet.size();k++){
				if(value_arm[k]>Amax){
										Amax=value_arm[k];

										}
				if(value_arm[k]<Amin){
												Amin=value_arm[k];

												}
		}
	cout<<"Amax :"<<Amax<<endl;
	cout<<"Amin :"<<Amin<<endl;
	double offset=(Amax+Amin)/2;








	//実験結果の初期値　と　データから読み取った初期値が合っているかどうか
	const double pen_ini=PEN_INI;
	const double pen_velocity_ini=PEN_VELOCITY_INI;

	const double arm_ini=ARM_INI;
	const double arm_velocity_ini=ARM_VELOCITY_INI;

	//初期値の確認
	cout<<"arm_ini :"<<arm_ini<<endl;
	cout<<"value_arm[0] :"<<value_arm[0]<<endl;

	if(arm_ini !=value_arm[0]){
										cerr << "初期値が違う" << endl;
										exit(0);
									}




	//実験データにオフセットを加える
	for(k=0;k<value_packet.size();k++){
		value_arm[k]=value_arm[k]-offset;
	}

	//シミュレーションの初期値を設定
	x[0][PEN]=pen_ini*enc_pendulum;//振り子角度
	x[0][PEN_VELOCITY]=pen_velocity_ini*enc_pendulum/0.004;//振り子角速度
	x[0][ARM]=arm_ini*enc_arm-offset*enc_arm;//アーム角度
	x[0][ARM_VELOCITY]=arm_velocity_ini*enc_arm/0.004;//アーム角速度



	xest[0][PEN]=pen_ini*enc_pendulum;//振り子角度
	xest[0][PEN_VELOCITY]=pen_velocity_ini*enc_pendulum/0.004;//振り子角速度
	xest[0][ARM]=arm_ini*enc_arm-offset*enc_arm;//アーム角度
	xest[0][ARM_VELOCITY]=arm_velocity_ini*enc_arm/0.004;//アーム角速度



	double RMSE=0;
	int packet_count=0;

								for(k=0;k<value_packet.size();k++){





									u[k] = K * (xref[k] - xest[k]);
									//電圧の分解能を含める
									u[k] = double(float(u[k]));
									if(value_packet[k]==1){
										    u[k]=0;
										    packet_count++;
										   // cout<<"error"<<endl;
										    }

									x[k + 1] = At * x[k] + Bt * u[k] ;

									//エンコーダの分解能を含める
									xest[k+1]=x[k+1];
									xest[k+1][ARM]=double (int(x[k+1][ARM]/(2*M_PI)*15000)*2*M_PI/15000);
									xest[k+1][ARM_VELOCITY]= double (int(x[k+1][ARM_VELOCITY]*0.004/(2*M_PI)*15000)*2*M_PI/15000/0.004);
									if(k<150){
									cout<<"arm angle experiment:"<<double(value_arm[k])*enc_arm<<"rad"<<endl;
									cout<<"arm angle simulation:"<<xest[k+1][ARM]<<"rad"<<endl<<endl;
									}
									xest[k+1][PEN]=double (int(x[k+1][PEN]/(2*M_PI)*10000)*2*M_PI/10000);
									xest[k+1][PEN_VELOCITY]= double (int(x[k+1][PEN_VELOCITY]*0.004/(2*M_PI)*15000)*2*M_PI/15000/0.004);


									RMSE=RMSE+(xest[k+1][ARM]-value_arm[k+1]*enc_arm)*(xest[k+1][ARM]-value_arm[k+1]*enc_arm);

									}
									cout<<"RMSE_raw :"<<RMSE<<endl;
									RMSE=sqrt(RMSE/value_packet.size());
									cout<<"RMSE:"<<RMSE<<"rad"<<endl;
									cout<<"RMSE:"<<RMSE*180/M_PI<<"deg"<<endl;
									cout<<"packe error count :"<<packet_count<<endl;


										u[k] = K * (xref[k] - xest[k]);



										ss.clear();
										ss.str("");


										ss<<n<<"_"<<PACKET<<"_output.csv";
										cout<<ss.str()<<endl;
										string Filename;
										Filename=ss.str();

										std::ofstream ofs(Filename);

										//ここから出力ファイルのフォーマットを決定する
										ofs << "#time (0.01s)\t" << "#simulation arm angle(deg)\t"<<"#experiment arm angle(deg)"<<endl;

										for(k=0;k<value_packet.size();k++){
												double tmp1=xest[k][ARM];
												double tmp2=value_arm[k];

												ofs <<k*0.01<< ","<<tmp1*180/M_PI<<","<<tmp2*enc_arm*180/M_PI<<endl;
											}

}

//アーム角度の実験結果のデータ確認用
//for(k=0;k<200;k++){
//	cout<<"arm angle :"<<value_arm[k]*enc_arm*180/M_PI<<endl;
//						}


										cout<<"------初期値-------"<<endl;

										cout<<"pen angle :"<<x[0][PEN]<<endl;
										cout<<"pen velocity :"<<x[0][PEN_VELOCITY]<<endl;
										cout<<"arm angle :"<<x[0][ARM]<<endl;


										cout<<"arm velocity :"<<x[0][ARM_VELOCITY]<<endl;


										return 0;
       	   	   	   }












