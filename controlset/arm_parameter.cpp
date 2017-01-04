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
    vector<double> mvalue;//�����ͤȺǾ��ͤ��Ǽ
    vector<double> mvalue_time;//�����ͤȺǾ��ͤ���������¸
    vector<double> JRR;
    vector<double> DRR;

    string str;
    const double PI=3.14159265;
    int i;
    int p;		// �ǥ�ߥ��λ����Ԥ�
    int q=0;		//�񤭹��ߤΤ���Υ��ե��åȤ�Ԥ�
    int f=0; //����ҡ�������γ�������λ�ɸ�ʻ����ɽ���Ƥ��뤳�Ȥˤʤ��
    int b=0; //��«�ͤ�٤��Ǥ�������õ������λ�ɸ, ����¾
    int c=0; //�����ͺǾ��ͤ����ݤ���ΰ�ư�˴ؤ��Ƥλ�ɸ
    int tmp_value;//�����͡��Ǿ��ͤΰ����Ǽ���
    int tmp_time; //�����͡��Ǿ��ͤ������ΰ����Ǽ���
    int n;//�롼�׻�ɸ
   string S;
   bool Z=false;// �񤭽Ф��ΤϤ������� true=�񤭹���
   bool N=false;//������ư�κ����ͺǾ��ͤ���뤿������Ѥ��롣
   double pam[3]; //�ɤ߼���ͤΰ����Ǽ���



	const double Kr=-40;//�ѥ�᡼����Ʊ���ɬ�פʥ��������
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
        // split��¹Ԥ��ޤ���
        list<string> strList = split(str, ",");

        // ���ƥ졼����������ޤ���
        list<string>::iterator iter = strList.begin();

        // ����Ƥߤޤ���
        while( iter != strList.end() ) // �Ǹ�ޤ�
        {
        	S=*iter;


            pam[p]=atof(S.c_str());
            // ���ĸ����˿ʤ�
            p++;
            ++iter;
        }


        if(pam[0]==1){Z=true;}//����ҡ�������ѥ�᡼��Ʊ������ꤷ�ʤ���Фʤ�ʤ�

        if(Z){

        //ofs << pam[0]<<","<< pam[1]<<","<<pam[2]<<endl;
        value_arm.push_back(pam[2]);//�٥��ȥ�ˤ��٤ƥ�����γ��٤򤪤���롣
        value_pendulum.push_back(pam[1]);
        value_time_real.push_back(pam[0]);
        f++;
       // cout<<f<<endl;


}
q++;
    }

    cout<<"size_t :"<<value_arm.size()<<endl;




b=0;

for(b=0;b<value_arm.size();b=b+100){//�����ͤ����Ĺ������ꤹ��

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
									  //�����ͤ����
										  if(tmp_value<value_arm[c]){//����Ҥȥ�����ǰ㤦

											 tmp_value=value_arm[c];

											  tmp_time=value_time_real[c];
									  }

								  }else{//�Ǿ��ͤ����
												  if(tmp_value>value_arm[c]){//����Ҥȥ�����ǰ㤦

													  tmp_value=value_arm[c];

													  tmp_time=value_time_real[c];
												  }



								  	  }

	  	  }
	  mvalue.push_back(tmp_value);
		  mvalue_time.push_back(tmp_time);
		  N=!N;//�����ͺǾ��ͤ��ߤ˵���

  }

  double sum;
  int sum_index;
  ofstream ofs("arm_100kHz_0.001_step_response.csv");

//  for(b=0;b<mvalue.size();b=b+1){
//	  cout<<"time :"<<mvalue_time[b]<<"   arm angle:"<<mvalue[b]<<endl;
//      ofs <<mvalue_time[b]/10000<<","<< double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;
//
//  }

cout<<"��«�� ������γ���:"<<value_arm[value_arm.size()-200]<<"   at:"<<value_arm.size()
		<<endl;//�Х��뤫�饪�ե��åȤ�ä���
cout<<"��«�͡�����ҳ���"<<value_pendulum[value_pendulum.size()-200]<<endl;
cout<<"end"<<endl;




//Jr��Dr�η׻�
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


cout<<"��������Jr"<<endl;
for(i=0;i<20;i++){
	cout<<JRR[i]<<endl;
}

cout<<"��������Dr"<<endl;
for(i=0;i<20;i++){
	cout<<DRR[i]<<endl;
}
cout<<"Kr :"<<KR<<endl;
return 0;

}





