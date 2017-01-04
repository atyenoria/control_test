// Time-stamp: <2012-12-26 17:00:16 kobayasi>

#ifndef __PENDULUM_H__
#define __PENDULUM_H__

#define PENDULUM_MODEL    "REALTECH RTC05"

#define PENDULUM_ANGLE    0 //����Ҥγ��٤�ɽ������ǥå���
#define ARM_ANGLE         1 //������γ��١�����
#define PENDULUM_VELOCITY 2 //����Ҥγ�®�١���
#define ARM_VELOCITY      3 //������γ�®�١���


inline
void pendulum(dMatrix<4,4> &A, dVector<4> &B) {
  /*
   * ��ž����Ω���ҤΥ�ǥ�
   * REALTECH RTC05
   */
  
  // ����(40cm)
  //const double m=0.054; //����夲���μ���[kg]
  //const double l=0.20;  //����夲����Ĺ����1/2[m]
  
  // ����(20cm)
  const double m=0.016; //����夲���μ���[kg]
  //const double l=0.10;  //����夲����Ĺ����1/2[m]
  const double l=0.085;
  //const double J=0.0048; //������β�ž������δ����⡼����[kgm^2]
  const double J=0;
  //const double r=0.20;   //�������Ĺ��[m]
  const double r=0.206;
  //const double g=9.81;   //���ϲ�®��[m/s^2]
  const double g=9.797;

  /*���
  // ����̤�������β�ž�ȥ륯�ξ��
  const double a = -m*r*g/J;
  const double b = 0;
  const double c = g*(J+m*r*r)/(l*J);
  const double d = 0;
  const double e = 1/J;
  const double f = -r/(l*J);
  */
  
  const double R=8.3;    //�⡼���ȥ륯ľή��[��]
  const double Km=0.023; //�⡼���ȥ륯���[Nm/A,Vs/rad]
  const double Kg=7.5;   //������ 120:16
  
  // ����̤�������⡼���������Ű��ξ��
  const double a = -m*r*g/J;
  const double b = -(Kg*Km)*(Kg*Km)/(R*J);
  const double c = g*(J+m*r*r)/(l*J);
  const double d = r*(Kg*Km)*(Kg*Km)/(R*l*J);
  const double e = Kg*Km/(R*J);
  const double f = -Kg*Km*r/(R*l*J);
  
  
  /*
   * ���ֶ��֥�ǥ�
   * x'(t) = A*x(t) + B*u(t)
   * y(t)  = C*x(t) + D*u(t)
   *
   * ���ϡ�
   * u(t) =  ������⡼���������Ű���or ������β�ž�ȥ륯��
   * ���֡�
   * x(t) = [����Ҥγ��١���(t),
   *         ������γ��١���(t),  % ���ͽ�ȤϽ�����ۤʤ�Τ���ա�
   *         ����Ҥγ�®�٦�'(t), %   �ץ����Ū�Թ��y[k][0,1]=x[k][0,1]��
   *         ������γ�®�٦�'(t)];
   */
  
  A = (double[4][4]){{ 0, 0, 1, 0},
		     { 0, 0, 0, 1},
		     { c, 0, 0, d},
		     { a, 0, 0, b}};
  B = (double[4]){0, 0, f, e};
  
  
  /* �ʤ������ͽ�Ǥϡ�
   * x(t) = [������γ��١���(t),
   *         ���Ҥγ��١�����(t),
   *         ������γ�®�٦�'(t),
   *         ���Ҥγ�®�١���'(t)];
   
  A = (double[4][4]){{ 0, 0, 1, 0},
		     { 0, 0, 0, 1},
		     { 0, a, b, 0},
		     { 0, c, d, 0}};
  B = (double[4]){0, 0, e, f};
  
  */
}


inline
void pendulum(dMatrix<4,4> &A, dMatrix<4,1> &B) {
  dVector<4> tmp;
  
  pendulum(A, tmp);
  for (int i=0; i<4; i++) B(i,0) = tmp(i);
}

#endif
