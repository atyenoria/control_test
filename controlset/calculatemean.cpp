/*
 * calculate mean.cpp
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
	    ifstream file("pendulum_100kHz.csv");
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
	    double tmp_mean=0;
	   string S;
	   bool Z=false;// 書き出しのはじめを決める true=書き込む
	   bool N=true;//振幅運動の最大値最小値を求めるために利用する。
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
	        tmp_mean=tmp_mean+pam[0];

	      f++;
}

	   cout<<tmp_mean/f<<endl;
	}

