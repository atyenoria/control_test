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
    vector<double> mvalue;//�����ͤȺǾ��ͤ��Ǽ
    vector<double> mvalue_time;//�����ͤȺǾ��ͤ���������¸

    string str;
    const double PI=3.14159265;
    int p;		// �ǥ�ߥ��λ����Ԥ�
    int q=0;		//�񤭹��ߤΤ���Υ��ե��åȤ�Ԥ�
    int f=0; //����ҡ�������γ�������λ�ɸ�ʻ����ɽ���Ƥ��뤳�Ȥˤʤ��
    int b=0; //��«�ͤ�٤��Ǥ�������õ������λ�ɸ, ����¾
    int c=0; //�����ͺǾ��ͤ����ݤ���ΰ�ư�˴ؤ��Ƥλ�ɸ
    int tmp_value;//�����͡��Ǿ��ͤΰ����Ǽ���
    int tmp_time; //�����͡��Ǿ��ͤ������ΰ����Ǽ���
   string S;
   bool Z=false;// �񤭽Ф��ΤϤ������� true=�񤭹���
   bool N=false;//������ư�κ����ͺǾ��ͤ���뤿������Ѥ��롣
   double pam[3]; //�ɤ߼���ͤΰ����Ǽ���

    if(file.fail()){
        cerr << "failed." << endl;
        exit(0);
    }


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
  for(b=0;b<mvalue.size();b=b+1	){
	  cout<<"time :"<<mvalue_time[b]<<"   arm angle:"<<mvalue[b]<<endl;
      ofs <<mvalue_time[b]/10000<<","<< double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;

  }

cout<<"��«�� ������γ���:"<<value_arm[value_arm.size()-200]<<"   at:"<<value_arm.size()
		<<endl;//�Х��뤫�饪�ե��åȤ�ä���
cout<<"��«�͡�����ҳ���"<<value_pendulum[value_pendulum.size()-200]<<endl;
cout<<"end"<<endl;




//Jr��Dr�η׻�
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
//Jp��Dp�η׻�


const double mp=0.016; //����Ҥμ���[kg]
const double lp=0.213;  //����Ҥ�Ĺ��[m]
const double Lcp=0.085; //����Ҥνſ��ޤǤε�Υ[m]
const double g=9.797;   //���ϲ�®��[m/s^2]
const double lamda=0.892; //����Ψ�ˡ��������ѿ�����������������
const double T=0.75; //����T[s]�����������ѿ���������
const double Lr=0.206;   //�������Ĺ��[m]
const double Jp=((mp*g*Lcp*T*T)/(4*PI*PI+log(lamda)*log(lamda))-mp*Lcp*Lcp);	//����Ҥνſ�����δ����⡼����[kg*m^2]
const double Dp=((-2*(Jp+mp*Lcp*Lcp)*log(lamda))/T);   //����Ҥβ�ž����Ǵ���໤[kg*m^2/s]

cout<<"mp:	"<<mp<<endl;
cout<<"lp:	"<<lp<<endl;
cout<<"Lcp:	"<<Lcp<<endl;
cout<<"lamda:	"<<lamda<<endl;
cout<<"T:	"<<T<<endl;
cout<<"Jp:	"<<Jp<<endl;
cout<<"Dp:	"<<Dp<<endl;




//������η׻���������������������


const double Jt=Jp*mp*Lr*Lr+Jr*Jp+Jr*mp*Lcp*Lcp;



//A����η�������ꤹ��
const double aa=mp*mp*Lcp*Lcp*Lr*g/Jt;
const double bb=-(Jp+mp*Lcp*Lcp)*(Dr+bm)/Jt;
const double cc=-mp*Lcp*Lr*Dp/Jt;
const double dd=mp*Lcp*g*(Jr+mp*Lr*Lr)/Jt;
const double ee=mp*Lcp*Lr*(Dr+bm)/Jt;
const double ff=-(Jr+mp*Lr*Lr)*Dp/Jt;


//B����η�������ꤹ��


const double gg=am*(Jp+mp*Lcp*Lcp)/Jt;
const double hh=am*mp*Lcp*Lr/Jt;


  dMatrix<4,4> A = (double[4][4]){{0,0,1,0},
								  {0,0,0,1},
								  {0,aa,bb,cc},
								  {0,dd,ee,ff}};

cout<<A[4][4]<<endl;
 dVector<4> B =(double[4]){0, 0, gg, hh};




//���ä��饲�����׻�����
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
 	cout<<"A����:"<<endl;
	printf("%lf",At);
	cout<<"B����"<<endl;
	printf("%lf",Bt);
	cout<<endl;
 	dscr(At, Bt, Ts);

	cout<<"discrete:A"<<endl;
	printf("%lf",At);
	cout<<"discrete:B"<<endl;
	printf("%lf",Bt);

 	const dVector<4> K = dlqr(At,Bt,Q,R);
 	cout<<endl<<endl;;
 	cout<<"��������١�����ҳ��١�������®�١������®��"<<endl;
 	cout<<"rad�ǤΥ�����"<<endl;
	printf("%lf",K);
	cout<<"theta�ǤΥ�����"<<endl;
	printf("%lf",K/180*PI);






return 0;
}

