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
	    double tmp_mean=0;
	   string S;
	   bool Z=false;// �񤭽Ф��ΤϤ������� true=�񤭹���
	   bool N=true;//������ư�κ����ͺǾ��ͤ���뤿������Ѥ��롣
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
	        tmp_mean=tmp_mean+pam[0];

	      f++;
}

	   cout<<tmp_mean/f<<endl;
	}

