/*
 * gnuplot.cpp
 *
 *  Created on: 2014/04/21
 *      Author: root
 */


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
#include "hrandom++.h"
#include "DTCS.h"
#include "PendulumR.h"
#include "LQ.h"
#define MODE 3
//1 for displaying objects and saving normal setting
//0 for writing eps
//2 for loading defualt setting
//3 for writing using loading

#define FILE_OUTPUT1 "reference_normal_control_biterorr_0.01"
#define FILE_OUTPUT2 "reference_omission_velocity_biterorr_0.01"
#define FILE_OUTPUT3 "reference_omission_angle_biterorr_0.01"
#define FILE_OUTPUT4 "reference_omission_velocity_biterorr_0.001"
#define FILE_OUTPUT5 "reference_omission_angle_biterorr_0.0001"
#define FILE_OUTPUT6 "reference_omission_angle_biterorr_0.001"
#define FILE_OUTPUT7 "reference_omission_angle_biterorr_0.0001"
#define FILE_OUTPUT8 "reference_omission_angle_biterorr_0.0001"
#define FILE_OUTPUT9 "reference_normal_control_biterorr_1e-05"
#define FILE_OUTPUT10 "reference_normal_control_biterorr_1e-05"
#define FILE_OUTPUT11 "reference_normal_control_biterorr_1e-05"
#define FILE_OUTPUT12 "reference_normal_control_biterorr_1e-05"
#define MULTI_PLOT true// 1 above 0 below

#define FILE_OUTPUTX "estimation0.01"





using namespace std;

int main(int argc, char **argv){

FILE* gnuplot=popen("gnuplot","w");

fprintf(gnuplot,

		" set grid\n"

		" set datafile separator\",\"\n"


		);

char str1[90]=FILE_OUTPUT1;
char str2[90]=FILE_OUTPUT2;
char str3[90]=FILE_OUTPUT3;
char str4[90]=FILE_OUTPUT4;
char str5[90]=FILE_OUTPUT5;
char str6[90]=FILE_OUTPUT6;
char str7[90]=FILE_OUTPUT7;
char str8[90]=FILE_OUTPUT8;
char str9[90]=FILE_OUTPUT9;
char str10[90]=FILE_OUTPUT10;
char str11[90]=FILE_OUTPUT11;
char str12[90]=FILE_OUTPUT12;
char strx[90]=FILE_OUTPUTX;
switch(MODE){

case 0:
if(MULTI_PLOT){
fprintf(gnuplot,"set terminal postscript eps enhanced\n");
fprintf(gnuplot,"set output \"step2a.eps\"\n");
fprintf(gnuplot,"set xlabel \"Packet Error Rate\"\n");
fprintf(gnuplot,"set ylabel \"ARM ANGLE RMSE(deg)\"\n");
fprintf(gnuplot,"set logscale x\n");
fprintf(gnuplot,"set multiplot\n");
fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\"\n",str1,str1,str2,str2,str3,str3);
fprintf(gnuplot,"quit\n");
}else{
fprintf(gnuplot,"set multiplot layout 3,1\n");
fprintf(gnuplot,"set logscale x\n");
fprintf(gnuplot,"plot \"%s\" using 1:2 with linespoints\n",strx);
fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints\n",strx);
fprintf(gnuplot,"plot \"%s\" using 1:4 with linespoints\n",strx);
}
break;




case 1:

	if(MULTI_PLOT){
	fprintf(gnuplot,
			" set grid\n"
			" set datafile separator\",\"\n"

			"set logscale x\n"
			"set multiplot\n"
			"set format e\n"
			"set key below\n"
			//"set yrange [8:13]n"
			"set xrange [0.001:1]\n"
			"set ytics 0.2\n"

			);



	//fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints linecolor rgbcolor \"red\" title \"%s\", \"%s\" using 1:3 with linespoints linecolor rgbcolor 'red'title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\",\"%s\" using 1:3 with linespoints linestyle 2 title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\"\n",str1,str1,str2,str2,str3,str3,str4,str4,str5,str5,str6,str6);
	fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints linecolor rgbcolor \"red\" title \"%s\",\"%s\" using 1:3 with linespoints linecolor rgbcolor \"blue\" title \"%s\",\"%s\" using 1:3 with linespoints linecolor rgbcolor \"green\" title \"%s\"\n",str2,str1,str1,str2,str3,str3);
	fprintf(gnuplot,"save \"save.plt\"\n");

	}else{
	fprintf(gnuplot,"set multiplot layout 3,1\n");
	fprintf(gnuplot,"set logscale x\n");
	fprintf(gnuplot,"plot \"%s\" using 1:2 with linespoints\n",strx);
	fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints\n",strx);
	fprintf(gnuplot,"plot \"%s\" using 1:4 with linespoints\n",strx);
	}

break;





case 2:



	if(MULTI_PLOT){
	fprintf(gnuplot,"load \"load.plt\"\n"
			);


	}else{
	fprintf(gnuplot,"set multiplot layout 3,1\n");
	fprintf(gnuplot,"set logscale x\n");
	fprintf(gnuplot,"plot \"%s\" using 1:2 with linespoints\n",strx);
	fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints\n",strx);
	fprintf(gnuplot,"plot \"%s\" using 1:4 with linespoints\n",strx);
	}

break;


case 3:
	fprintf(gnuplot,"set output \"graph2.png\"\n");
	//fprintf(gnuplot,"set terminal postscript eps enhanced color\n");
	fprintf(gnuplot,"set terminal png\n");


	if(MULTI_PLOT){
	fprintf(gnuplot,"load \"load.plt\"\n"
			);


	fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints linecolor rgbcolor \"red\" title \"%s\", \"%s\" using 1:3 with linespoints linetype 2 linecolor rgbcolor 'red'title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\",\"%s\" using 1:3 with linespoints linestyle 2 title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\", \"%s\" using 1:3 with linespoints title \"%s\"\n",str1,str1,str2,str2,str3,str3,str4,str4,str5,str5,str6,str6);



	}else{
	fprintf(gnuplot,"set multiplot layout 3,1\n");
	fprintf(gnuplot,"set logscale x\n");
	fprintf(gnuplot,"plot \"%s\" using 1:2 with linespoints\n",strx);
	fprintf(gnuplot,"plot \"%s\" using 1:3 with linespoints\n",strx);
	fprintf(gnuplot,"plot \"%s\" using 1:4 with linespoints\n",strx);
	}



break;}














fflush(gnuplot);



fprintf(gnuplot,"quit\n");
pclose(gnuplot);



time_t rawtime;
  struct tm * timeinfo;
cout<<"\n\n\n"<<endl;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  printf ("Current local time and date: %s", asctime(timeinfo));




}
