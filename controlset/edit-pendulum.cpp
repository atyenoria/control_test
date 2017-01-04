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
    vector<double> mvalue;//�����ͤȺǾ��ͤ��Ǽ
    vector<double> mvalue_time;//�����ͤȺǾ��ͤ���������¸

    string str;
    int p;		// �ǥ�ߥ��λ����Ԥ�
    int q=0;		//�񤭹��ߤΤ���Υ��ե��åȤ�Ԥ�
    int f=0; //����ҡ�������γ�������λ�ɸ�ʻ����ɽ���Ƥ��뤳�Ȥˤʤ��
    int b=0; //��«�ͤ�٤��Ǥ�������õ������λ�ɸ, ����¾
    int c=0; //�����ͺǾ��ͤ����ݤ���ΰ�ư�˴ؤ��Ƥλ�ɸ
    int tmp_value;//�����͡��Ǿ��ͤΰ����Ǽ���
    int tmp_time; //�����͡��Ǿ��ͤ������ΰ����Ǽ���
   string S;
   bool Z=false;// �񤭽Ф��ΤϤ������� true=�񤭹���
   bool N=true;//������ư�κ����ͺǾ��ͤ���뤿������Ѥ��롣
   double pam[3]; //�ɤ߼���ͤΰ����Ǽ���

    if(file.fail()){
        cerr << "failed." << endl;
        exit(0);}

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
if(q>3000){Z=true;
        if(pam[0]==5000){Z=true;}//����ҡ�������ѥ�᡼��Ʊ������ꤷ�ʤ���Фʤ�ʤ�

        if(Z){

        //ofs << pam[0]<<","<< pam[1]<<","<<pam[2]<<endl;
        value_arm.push_back(pam[2]);//�٥��ȥ�ˤ��٤ƥ�����γ��٤򤪤���롣
        value_pendulum.push_back(pam[1]);
        value_time_real.push_back(pam[0]);

    }
}
q++;
    }

b=0;

for(b=0;b<value_arm.size();b=b+100){//�����ͤ����Ĺ������ꤹ��


    if( ((3<value_pendulum[b]) && (3>value_pendulum[b+100]))||
		    ((3>value_pendulum[b])&&(3<value_pendulum[b+100]))){


	interval_index.push_back(b);
}
}

for(b=0;b<interval_index.size();b++){
	cout<<interval_index[b]<<endl;
}

   //�����ͤη׻��Ȥ��ι���������
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

								  }else{//�Ǿ��ͤ����

												  if(tmp_value>value_pendulum[c]){
													   tmp_value=value_pendulum[c];
													  tmp_time=tmp_time=value_time_real[c];
												  }


								  	  }

	  	  }
	  mvalue.push_back(tmp_value);
		  mvalue_time.push_back(tmp_time);
		  N=!N;//�����ͺǾ��ͤ��ߤ˵���

  }

//  double sum;
//  int sum_index;
  ofstream ofs("apendulum_100kHz_attenuation_factor.csv");
  for(b=0;b<mvalue.size();b=b+1	){
	  cout<<"time :"<<mvalue_time[b]<<"   pendulum angle:"<<mvalue[b]<<endl;
      cout<<"�ֳ֡���"<<mvalue_time[b]-mvalue_time[b+2]<<"  �����桡:"<<double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;
      ofs <<mvalue_time[b]/10000<<","<< double(double(abs(int(mvalue[b+2])))/double(abs(int(mvalue[b]))))<<endl;
  }

cout<<"��«�� ������γ���:"<<value_arm[value_arm.size()-200]<<"   "<<value_arm.size()
		<<endl;//�Х��뤫�饪�ե��åȤ�ä���
cout<<"��«�͡�����ҳ���"<<value_pendulum[value_pendulum.size()-200]<<endl;
cout<<"end"<<endl;



const double PI=3.14159265;
const double mp=0.016; //����Ҥμ���[kg]
const double lp=0.213;  //����Ҥ�Ĺ��[m]
const double Lcp=0.085; //����Ҥνſ��ޤǤε�Υ[m]
const double g=9.797;   //���ϲ�®��[m/s^2]
const double lamda=0.888; //����Ψ��
const double T=(mvalue_time[7]-mvalue_time[5])*0.0001; //����T[s]

const double Jp=(mp*g*Lcp*T*T)/(4*PI*PI+log(lamda)*log(lamda))-mp*Lcp*Lcp;	//����Ҥνſ�����δ����⡼����[kg*m^2]
const double Dp=(-2*(Jp+mp*Lcp*Lcp)*log(lamda))/T;   //����Ҥβ�ž����Ǵ���໤[kg*m^2/s]




cout<<"mp:	"<<mp<<endl;
cout<<"lp:	"<<lp<<endl;
cout<<"Lcp:	"<<Lcp<<endl;
cout<<"lamda:	"<<lamda<<endl;
cout<<"T:	"<<T<<endl;
cout<<"Jp:	"<<Jp<<endl;
cout<<"Dp:	"<<Dp<<endl;


    return 0;
}

