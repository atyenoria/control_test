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
    ifstream file("pendulum_1_0_100kHz.csv");
    //vector<vector<int> > values;
    vector<double> value_arm;
    vector<double> value_pendulum;
    vector<double> value_time_real;
    vector<double> interval_index;
    vector<double> mvalue;//最大値と最小値を格納
    vector<double> mvalue_time;//最大値と最小値を取る時刻を保存

    string str;
    int p;		// デリミタの指定を行う
    int q=0;		//書き込みのためのオフセットを行う
    int f=0; //振り子、アームの角度配列の指標（時刻を表していることになる）
    int b=0; //収束値を跨いでいる点を探すための指標, その他
    int c=0; //最大値最小値を求める際の窓の移動に関しての指標
    int tmp_value;//最大値、最小値の一時格納場所
    int tmp_time; //最大値、最小値を取る時刻の一時格納場所
   string S;
   bool Z=false;// 書き出しのはじめを決める true=書き込む
   bool N=true;//振幅運動の最大値最小値を求めるために利用する。
   double pam[3]; //読み取り値の一時格納場所

    if(file.fail()){
        cerr << "failed." << endl;
        exit(0);}

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
if(q>3000){Z=true;
        if(pam[0]==5000){Z=true;}//振り子、アームパラメータ同定で設定しなければならない

        if(Z){

        //ofs << pam[0]<<","<< pam[1]<<","<<pam[2]<<endl;
        value_arm.push_back(pam[2]);//ベクトルにすべてアームの角度をおさめる。
        value_pendulum.push_back(pam[1]);
        value_time_real.push_back(pam[0]);

    }
}
q++;
    }

b=0;

for(b=0;b<value_arm.size();b=b+100){//この値で窓の長さを決定する


    if( ((3<value_pendulum[b]) && (3>value_pendulum[b+100]))||
		    ((3>value_pendulum[b])&&(3<value_pendulum[b+100]))){


	interval_index.push_back(b);
}
}

for(b=0;b<interval_index.size();b++){
	cout<<interval_index[b]<<endl;
}

   //最大値の計算とその刻み幅を求める
   cout<<"\n"<<value_pendulum.size()<<endl;

  for(b=0;b<interval_index.size();b++){
	  tmp_value=0;
	  for(c=int(interval_index[b]);
			  c<int(interval_index[b+1]);
			  c++){
								  if(N){

									  if(tmp_value<value_pendulum[c]){
										        tmp_value=value_pendulum[c];
											  tmp_time=value_time_real[c];
									  }

								  }else{//最小値を求める

												  if(tmp_value>value_pendulum[c]){
													   tmp_value=value_pendulum[c];
													  tmp_time=tmp_time=value_time_real[c];
												  }


								  	  }

	  	  }
	  mvalue.push_back(tmp_value);
		  mvalue_time.push_back(tmp_time);
		  N=!N;//最大値最小値を交互に求める

  }

//  double sum;
//  int sum_index;
  ofstream ofs("apendulum_100kHz_attenuation_factor.csv");
  for(b=0;b<mvalue.size();b=b+1	){
	  cout<<"time :"<<mvalue_time[b]<<"   pendulum angle:"<<mvalue[b]<<endl;
      cout<<"間隔　："<<mvalue_time[b]-mvalue_time[b+2]<<"  減衰比　:"<<double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;
      ofs <<mvalue_time[b]/10000<<","<< double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;
  }

cout<<"収束値 アームの角度:"<<value_arm[value_arm.size()-200]<<"   "<<value_arm.size()
		<<endl;//バグルからオフセットを加えた
cout<<"収束値　振り子角度"<<value_pendulum[value_pendulum.size()-200]<<endl;
cout<<"end"<<endl;



const double PI=3.14159265;
const double mp=0.016; //振り子の質量[kg]
const double lp=0.213;  //振り子の長さ[m]
const double Lcp=0.085; //振り子の重心までの距離[m]
const double g=9.797;   //重力加速度[m/s^2]
const double lamda=0.888; //減衰率λ
const double T=(mvalue_time[7]-mvalue_time[5])*0.0001; //周期T[s]

const double Jp=(mp*g*Lcp*T*T)/(4*PI*PI+log(lamda)*log(lamda))-mp*Lcp*Lcp;	//振り子の重心周りの慣性モーメント[kg*m^2]
const double Dp=(-2*(Jp+mp*Lcp*Lcp)*log(lamda))/T;   //振り子の回転軸の粘性摩擦[kg*m^2/s]




cout<<"mp:	"<<mp<<endl;
cout<<"lp:	"<<lp<<endl;
cout<<"Lcp:	"<<Lcp<<endl;
cout<<"lamda:	"<<lamda<<endl;
cout<<"T:	"<<T<<endl;
cout<<"Jp:	"<<Jp<<endl;
cout<<"Dp:	"<<Dp<<endl;


    return 0;
}

