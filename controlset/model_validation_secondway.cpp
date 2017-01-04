/*
 * model_validation.cpp
 *
 *  Created on: Jun 28, 2014
 *      Author: atyenoria
 */




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
#define ARM 1
#define PEN 0
#define ARM_VELOCITY 3
#define PEN_VELOCITY 2




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
    vector<double> mvalue;//�����ͤȺǾ��ͤ��Ǽ
    vector<double> mvalue_time;//�����ͤȺǾ��ͤ���������¸
    vector<double> JRR;
    vector<double> DRR;


    dVector<4> x[10000];

    string str;
    const double PI=3.14159265;
    int p;		// �ǥ�ߥ��λ����Ԥ�
    int q=0;		//�񤭹��ߤΤ���Υ��ե��åȤ�Ԥ�
    int f=0; //����ҡ�������γ�������λ�ɸ�ʻ����ɽ���Ƥ��뤳�Ȥˤʤ��
    int b=0; //��«�ͤ�٤��Ǥ�������õ������λ�ɸ, ����¾
    int c=0; //�����ͺǾ��ͤ����ݤ���ΰ�ư�˴ؤ��Ƥλ�ɸ
    int tmp_value;//�����͡��Ǿ��ͤΰ����Ǽ���
    int tmp_time; //�����͡��Ǿ��ͤ������ΰ����Ǽ���
    int n;//�롼�׻�ɸ
    int k;//�롼���Ѥ�k
    string S;
    bool Z=false;// �񤭽Ф��ΤϤ������� true=�񤭹���
    bool N=false;//������ư�κ����ͺǾ��ͤ���뤿������Ѥ��롣
    double pam[3]; //�ɤ߼���ͤΰ����Ǽ���



	const double Kr=-90;//�ѥ�᡼����Ʊ���ɬ�פʥ��������
	double KR=abs(Kr);


for(n=0;n<1;n++){

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


					switch (MODE) {

									case 0:// no error


											value_arm.push_back(pam[2]);//�٥��ȥ�ˤ��٤ƥ�����γ��٤򤪤���롣
											value_pendulum.push_back(pam[1]);
											value_time_real.push_back(pam[0]);
											f++;
											break;

									case 1:// error

											value_packet.push_back(pam[2]);
											value_arm.push_back(pam[1]);
											value_time_real.push_back(pam[0]);
											f++;
											break;


								   }

               }
    }

    cout<<"size_t :"<<value_arm.size()<<endl;


//���ĤΥѥ�᡼��Ʊ���ɬ�פȤ��륷���ƥ�Ʊ��ˡ
const double am=0.0207831;
const double bm=0.00358509;
const double Jr=0.005023319;//������δ����⡼����
const double Dr=0.0133366;//�������Ǵ���໤


cout<<"am:	"<<am<<endl;
cout<<"bm:	"<<bm<<endl;
cout<<"Jr:	"<<Jr<<endl;
cout<<"Dr:	"<<Dr<<endl;
cout<<"--------------------"<<endl;


const double mp=0.016; //����Ҥμ���[kg]
const double lp=0.213;  //����Ҥ�Ĺ��[m]
const double Lcp=0.085; //����Ҥνſ��ޤǤε�Υ[m]
const double g=9.797;   //���ϲ�®��[m/s^2]
const double lamda=0.92; //����Ψ�ˡ��������ѿ�����������������
const double T=0.75; //����T[s]�����������ѿ���������
const double Lr=0.206;   //�������Ĺ��[m]
const double Jp=7.41803e-05;	//����Ҥνſ�����δ����⡼����[kg*m^2]
const double Dp=5.78395e-05;   //����Ҥβ�ž����Ǵ���໤[kg*m^2/s]



cout<<"mp:	"<<mp<<endl;
cout<<"lp:	"<<lp<<endl;
cout<<"Lcp:	"<<Lcp<<endl;
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

//���ֹ���ν��֤��ѹ����뢪������ҡ������ࡡ�����®�١�������®�١���

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

 	dMatrix<4, 4> Q = (double[4][4] ) { { 2, 0, 0, 0 },
										{ 0, 10, 0, 0 },
										{ 0, 0,0, 0 },
										{ 0, 0, 0, 0.5 } };

 	double R = 0.005;

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
 	cout<<"������ҳ��� ��������١������®�� ������®��"<<endl;
 	cout<<"rad�ǤΥ�����"<<endl;
	printf("%lf",K);
	cout<<"theta�ǤΥ�����"<<endl;
	printf("%lf",K/180*PI);







	//�����������泫�ϡ�����ͤ����ꤷ�ƥ��ߥ�졼������Ԥ�

	const int Sample_Range=100;//����ץ���ߥ���ץ�󥰻���

	dVector<4> xref[Sample_Range + 2];

		for (k = 0; k <= Sample_Range; k++) {
			xref[k] = 0.0;
			xref[k][ARM] = M_PI/2;
		}


	int i;
	dVector<4> xest[Sample_Range + 2];
	double u[Sample_Range + 2];

	//����ͤ�����ʾ��־���ν��֡�����ҳ��١���������١�����ҳ�®�١��������®�١��ν��֡�
	x[0][ARM]=0;//���������
	x[0][PEN]=0;//����ҳ���
	x[0][ARM_VELOCITY]=0;//�������®��
	x[0][PEN_VELOCITY]=0;//����ҳ�®��


	xest[0][ARM]=0;//���������
	xest[0][PEN]=0;//����ҳ���
	xest[0][ARM_VELOCITY]=0;//�������®��
	xest[0][PEN_VELOCITY]=0;//����ҳ�®��

	double RMSE=0;

	switch (MODE) {

		         case 0:

								for(k=0;k<Sample_Range;k++){



									u[k] = K * (xref[k] - xest[k]);
									//�Ű���ʬ��ǽ��ޤ��
									u[k] = double(float(u[k]));
									x[k + 1] = At * x[k] + Bt * u[k] ;

									//���󥳡�����ʬ��ǽ��ޤ��
									xest[k+1]=x[k+1];
								//	cout<<"before: "<<xest[k+1][0]<<endl;
									xest[k+1][ARM]=double (int(x[k+1][0]/(2*M_PI)*15000)*2*M_PI/15000);
									//cout<<"after arm: "<<xest[k+1][ARM]<<endl;
									//cout<<"after pendulum: "<<xest[k+1][PEN]<<endl;
									xest[k+1][ARM_VELOCITY]= double (int(x[k+1][2]*0.004/(2*M_PI)*15000)*2*M_PI/15000/0.004);


									xest[k+1][PEN]=double (int(x[k+1][1]/(2*M_PI)*10000)*2*M_PI/10000);
									xest[k+1][PEN_VELOCITY]= double (int(x[k+1][3]*0.004/(2*M_PI)*15000)*2*M_PI/15000/0.004);


									RMSE=(xest[k+1][ARM]-value_arm[k+1])*(xest[k+1][ARM]-value_arm[k+1]);

									}

									RMSE=sqrt(RMSE/Sample_Range);
									cout<<"RMSE:"<<RMSE<<endl;


									break;


									case 1:


										for(k=0;k<Sample_Range;k++){

										u[k] = K * (xref[k] - xest[k]);
										if(value_packet[k]==1){
											u[k]=0;
										}

										//�Ű���ʬ��ǽ��ޤ��
										u[k] = double(float(u[k]));
										x[k + 1] = At * x[k] + Bt * u[k] ;

										//���󥳡�����ʬ��ǽ��ޤ��
										xest[k+1]=x[k+1];
									//	cout<<"before: "<<xest[k+1][0]<<endl;
										xest[k+1][ARM]=double (int(x[k+1][0]/(2*M_PI)*15000)*2*M_PI/15000);
										cout<<"after penaaaaadulum: "<<xest[k+1][0]<<endl;
										cout<<"after pendulum: "<<xest[k+1][1]<<endl;
										xest[k+1][ARM_VELOCITY]= double (int(x[k+1][2]*0.004/(2*M_PI)*15000)*2*M_PI/15000/0.004);


										xest[k+1][PEN]=double (int(x[k+1][1]/(2*M_PI)*10000)*2*M_PI/10000);
										xest[k+1][PEN_VELOCITY]= double (int(x[k+1][3]*0.004/(2*M_PI)*15000)*2*M_PI/15000/0.004);


										RMSE=(xest[k+1][ARM]-value_arm[k+1])*(xest[k+1][ARM]-value_arm[k+1]);

										}

										RMSE=sqrt(RMSE/Sample_Range);
										cout<<"RMSE:"<<RMSE<<endl;

										break;
										}

}



ss.clear();
ss.str("");

ss << "reference_normal_control_biterorr_"
		<< sigma_w[number_disturbance];
std::ofstream ofs(ss.str().c_str());
ofs << "#packet error rate\t" << "#fall down rate\t"
		<< "#arm angel EMSE\t" << "#pendulum angel RMSE\t"
		<< "#all degree" << endl;
for (i = 0; i < count_error_index; i++) {
	for (k = 0; k < 4; k++) {
		ofs << error[i][k] << ",";
	}
	ofs << endl;
}
ofs.close();


return 0;

}





