/*

 * edit.cpp
 *
 *  Created on: May 21, 2014
 *      Author: atyenoria
 */
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
using namespace std;


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
		ss.str("");
		ss<<1<<"_"<<90<<"_"<<1<<"kHz.csv";
		const double Kr=-200;
		cout<<ss.str()<<endl;
		string filename;
		filename=ss.str();
    ifstream file(filename);
    //vector<vector<int> > values;
    vector<double> value_arm;
    vector<double> value_pendulum;
    vector<double> value_time_real;
    vector<double> interval_index;
    vector<double> mvalue;//最大値と最小値を格納
    vector<double> mvalue_time;//最大値と最小値を取る時刻を保存

    string str;
    const double PI=3.14159265;
    int p;		// デリミタの指定を行う
    int q=0;		//書き込みのためのオフセットを行う
    int f=0; //振り子、アームの角度配列の指標（時刻を表していることになる）
    int b=0; //収束値を跨いでいる点を探すための指標, その他
    int c=0; //最大値最小値を求める際の窓の移動に関しての指標
    int tmp_value;//最大値、最小値の一時格納場所
    int tmp_time; //最大値、最小値を取る時刻の一時格納場所
   string S;
   bool Z=false;// 書き出しのはじめを決める true=書き込む
   bool N=false;//振幅運動の最大値最小値を求めるために利用する。
   double pam[3]; //読み取り値の一時格納場所

    if(file.fail()){
        cerr << "failed." << endl;
        exit(0);
    }


    while(getline(file, str)){
p=0;
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


        if(pam[0]==1){Z=true;}//振り子、アームパラメータ同定で設定しなければならない

        if(Z){

        //ofs << pam[0]<<","<< pam[1]<<","<<pam[2]<<endl;
        value_arm.push_back(pam[2]);//ベクトルにすべてアームの角度をおさめる。
        value_pendulum.push_back(pam[1]);
        value_time_real.push_back(pam[0]);
        f++;
       // cout<<f<<endl;


}
q++;
    }

    cout<<"size_t :"<<value_arm.size()<<endl;




b=0;

for(b=0;b<value_arm.size();b=b+100){//この値で窓の長さを決定する

							if( ((3<value_arm[b]) && (3>value_arm[b+100]))||
									((3>value_arm[b])&&(3<value_arm[b+100]))){


										interval_index.push_back(b);
															}
					}

for(b=0;b<interval_index.size();b++){
	cout<<interval_index[b]<<endl;
}


  for(b=0;b<interval_index.size();b++){
	  tmp_value=0;

	  for(c=int(interval_index[b]);
			  c<int(interval_index[b+1]);
			  c++){
								  if(N){
									  //最大値を求める
										  if(tmp_value<value_arm[c]){//振り子とアームで違う

											 tmp_value=value_arm[c];

											  tmp_time=value_time_real[c];
									  }

								  }else{//最小値を求める
												  if(tmp_value>value_arm[c]){//振り子とアームで違う

													  tmp_value=value_arm[c];

													  tmp_time=value_time_real[c];
												  }



								  	  }

	  	  }
	  mvalue.push_back(tmp_value);
		  mvalue_time.push_back(tmp_time);
		  N=!N;//最大値最小値を交互に求める

  }

  double sum;
  int sum_index;
  ofstream ofs("arm_100kHz_0.001_step_response.csv");
  for(b=0;b<mvalue.size();b=b+1	){
	  cout<<"time :"<<mvalue_time[b]<<"   arm angle:"<<mvalue[b]<<endl;
      ofs <<mvalue_time[b]/10000<<","<< double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;

  }

cout<<"収束値 アームの角度:"<<value_arm[value_arm.size()-200]<<"   at:"<<value_arm.size()
		<<endl;//バグルからオフセットを加えた
cout<<"収束値　振り子角度"<<value_pendulum[value_pendulum.size()-200]<<endl;
cout<<"end"<<endl;




//JrとDrの計算
const double am=0.0207831;
const double bm=0.00358509;
const double Tp=mvalue_time[0]*0.001;
const double Amax=double(abs(int(mvalue[0])-int(value_arm[value_arm.size()-200])))/15000*2*PI;
cout<<"arm[0]"<<value_arm[0]<<endl;

const double theta=double(abs(int(value_arm[2])-int(value_arm[value_arm.size()-200])))/15000*2*PI;;
const double Jr=0.00343267; //double(double(am*Kr*(Tp*Tp))/double(PI*PI+log(Amax/theta))*(log(Amax/theta)));
const double Dr=0.00424722;//(-2*Jr*log(Amax/theta))/Tp-bm;
cout<<"file name:"<<ss.str()<<endl;
cout<<"am:	"<<am<<endl;
cout<<"bm:	"<<bm<<endl;
cout<<"Tp:	"<<Tp<<endl;
cout<<"Amax:	"<<Amax<<endl;
cout<<"Kr:	"<<Kr<<endl;
cout<<"theta:	"<<theta<<endl;
cout<<"Jr:	"<<Jr<<endl;
cout<<"Dr:	"<<Dr<<endl;

cout<<"--------------------"<<endl;
//JpとDpの計算


const double mp=0.016; //振り子の質量[kg]
const double lp=0.213;  //振り子の長さ[m]
const double Lcp=0.085; //振り子の重心までの距離[m]
const double g=9.797;   //重力加速度[m/s^2]
const double lamda=0.892; //減衰率λ＊＊＊＊変数＊＊＊＊＊＊＊＊
const double T=0.75; //周期T[s]＊＊＊＊＊変数＊＊＊＊
const double Lr=0.206;   //アームの長さ[m]
const double Jp=((mp*g*Lcp*T*T)/(4*PI*PI+log(lamda)*log(lamda))-mp*Lcp*Lcp);	//振り子の重心周りの慣性モーメント[kg*m^2]
const double Dp=((-2*(Jp+mp*Lcp*Lcp)*log(lamda))/T);   //振り子の回転軸の粘性摩擦[kg*m^2/s]

cout<<"mp:	"<<mp<<endl;
cout<<"lp:	"<<lp<<endl;
cout<<"Lcp:	"<<Lcp<<endl;
cout<<"lamda:	"<<lamda<<endl;
cout<<"T:	"<<T<<endl;
cout<<"Jp:	"<<Jp<<endl;
cout<<"Dp:	"<<Dp<<endl;




//ゲインの計算＊＊＊＊＊＊＊＊＊＊


const double Jt=Jp*mp*Lr*Lr+Jr*Jp+Jr*mp*Lcp*Lcp;



//A行列の係数を決定する
const double aa=mp*mp*Lcp*Lcp*Lr*g/Jt;
const double bb=-(Jp+mp*Lcp*Lcp)*(Dr+bm)/Jt;
const double cc=-mp*Lcp*Lr*Dp/Jt;
const double dd=mp*Lcp*g*(Jr+mp*Lr*Lr)/Jt;
const double ee=mp*Lcp*Lr*(Dr+bm)/Jt;
const double ff=-(Jr+mp*Lr*Lr)*Dp/Jt;


//B行列の係数を決定する


const double gg=am*(Jp+mp*Lcp*Lcp)/Jt;
const double hh=am*mp*Lcp*Lr/Jt;


  dMatrix<4,4> A = (double[4][4]){{0,0,1,0},
								  {0,0,0,1},
								  {0,aa,bb,cc},
								  {0,dd,ee,ff}};

cout<<A[4][4]<<endl;
 dVector<4> B =(double[4]){0, 0, gg, hh};




//こっからゲインを計算する
 	const int seed = 1;
 	const double Hz = 100;
 	const double Ts = 1.0 / Hz;


 	dMatrix<4, 4> At=A;
 	dVector<4> Bt=B;

 	dMatrix<4, 4> Q = (double[4][4] ) { { 8, 0, 0, 0 },
										{ 0, 10, 0, 0 },
										{ 0, 0,0.5, 0 },
										{ 0, 0, 0, 0 } };

 	double R = 0.006;
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
 	cout<<endl<<endl;;
 	cout<<"アーム角度　振り子角度　アーム速度　振り子速度"<<endl;
 	cout<<"radでのゲイン"<<endl;
	printf("%lf",K);
	cout<<"thetaでのゲイン"<<endl;
	printf("%lf",K/180*PI);






return 0;
}

