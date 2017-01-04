/*
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

    //vector<vector<int> > values;
    vector<double> value_arm;
    vector<double> value_pendulum;
    vector<double> value_time_real;
    vector<double> interval_index;
    vector<double> mvalue;//最大値と最小値を格納
    vector<double> mvalue_time;//最大値と最小値を取る時刻を保存
    vector<double> JRR;
    vector<double> DRR;

    string str;
    const double PI=3.14159265;
    int i;
    int p;		// デリミタの指定を行う
    int q=0;		//書き込みのためのオフセットを行う
    int f=0; //振り子、アームの角度配列の指標（時刻を表していることになる）
    int b=0; //収束値を跨いでいる点を探すための指標, その他
    int c=0; //最大値最小値を求める際の窓の移動に関しての指標
    int tmp_value;//最大値、最小値の一時格納場所
    int tmp_time; //最大値、最小値を取る時刻の一時格納場所
    int n;//ループ指標
   string S;
   bool Z=false;// 書き出しのはじめを決める true=書き込む
   bool N=false;//振幅運動の最大値最小値を求めるために利用する。
   double pam[3]; //読み取り値の一時格納場所



	const double Kr=-40;//パラメーター同定に必要なゲイン定数
	double KR=abs(Kr);


for(n=0;n<10;n++){

	ss.str("");
	ss<<n<<"_"<<KR<<"_"<<1<<"kHz.csv";
	cout<<ss.str()<<endl;
	string filename;
	filename=ss.str();
	ifstream file(filename);

    if(file.fail()){
        cerr << "failed." << endl;
        exit(0);
    }

    Z=false;
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

 N=false;

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

//  for(b=0;b<mvalue.size();b=b+1){
//	  cout<<"time :"<<mvalue_time[b]<<"   arm angle:"<<mvalue[b]<<endl;
//      ofs <<mvalue_time[b]/10000<<","<< double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;
//
//  }

cout<<"収束値 アームの角度:"<<value_arm[value_arm.size()-200]<<"   at:"<<value_arm.size()
		<<endl;//バグルからオフセットを加えた
cout<<"収束値　振り子角度"<<value_pendulum[value_pendulum.size()-200]<<endl;
cout<<"end"<<endl;




//JrとDrの計算
const double am=0.0207831;
const double bm=0.00358509;
const double Tp=mvalue_time[0]*0.001;
const double Amax=double(abs(int(mvalue[0])))/15000*2*PI;
cout<<"arm[0]"<<value_arm[0]<<endl;
const double theta=double(abs(int(value_arm[2])))/15000*2*PI;
const double Jr=double(double(am*Kr*(Tp*Tp))/double(PI*PI+log(Amax/theta))*(log(Amax/theta)));
const double Dr=(-2*Jr*log(Amax/theta))/Tp-bm;
cout<<"file name:"<<ss.str()<<endl;
cout<<value_arm[2]<<endl;
cout<<"am:	"<<am<<endl;
cout<<"bm:	"<<bm<<endl;
cout<<"Tp:	"<<Tp<<endl;
cout<<"Amax:	"<<Amax<<endl;
cout<<"Kr:	"<<Kr<<endl;
cout<<"theta:	"<<theta<<endl;
cout<<"Jr:	"<<Jr<<endl;
cout<<"Dr:	"<<Dr<<endl;
cout<<"--------------------"<<endl;
JRR.push_back(Jr);
DRR.push_back(Dr);


for(i=0;i<100;i++){
	cout<<"arm_value :"<<value_arm[i]<<endl;
}
}
cout<<JRR[0]<<endl;
cout<<JRR[1]<<endl;


cout<<"ここからJr"<<endl;
for(i=0;i<20;i++){
	cout<<JRR[i]<<endl;
}

cout<<"ここからDr"<<endl;
for(i=0;i<20;i++){
	cout<<DRR[i]<<endl;
}
cout<<"Kr :"<<KR<<endl;
return 0;

}





