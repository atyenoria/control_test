// Time-stamp: <2012-12-26 16:59:57 kobayasi>

#ifndef __PENDULUM_H__
#define __PENDULUM_H__

#define PENDULUM_MODEL    "N.J.Ploplys"

#define PENDULUM_ANGLE    0 //����Ҥγ��٤�ɽ������ǥå���
#define ARM_ANGLE         1 //������γ��١�����
#define PENDULUM_VELOCITY 2 //����Ҥγ�®�١���
#define ARM_VELOCITY      3 //������γ�®�١���


inline
void pendulum(dMatrix<4,4> &A, dVector<4> &B) {
  /*
   * ��ž����Ω���ҤΥ�ǥ�
   * N. J. Ploplys, "Wireless Feedback Control of Mechanical Systems"
   * Master Thesis, Department of Mechanical Engineering,
   *   University of Illinois, Champaign, IL, 2003
   */
  
  const double M=0.004;  //����ҤΤ����μ���[kg]
  const double m=0.025;  //����夲���μ���[kg]
  const double l=0.241;  //����夲����Ĺ��[m]
  const double r=0.152;  //�������Ĺ��[m]
  const double J=1.21e-3;//������β�ž������δ����⡼����[kgm^2]
  const double g=9.81;   //���ϲ�®��[m/s^2]
  
  const double a=(M+m/3)*l*l;
  const double b=J+(M+m)*r*r;
  const double c=(M+m/2)*r*l;
  const double d=(M+m/2)*g*l;
  
  /*
   * Ϣ³���֥����ƥ����ʬ�����������������
   * a*��''(t) + c*��''(t) -d*��'(t) = 0
   * c*��''(t) + b*��''(t) = u(t)
   *
   * ���ֶ��֥�ǥ�
   * x'(t) = A*x(t) + B*u(t)
   * y(t)  = C*x(t) + D*u(t)
   *
   * ���ϡ�
   * u(t) =  ������β�ž�ȥ륯
   * ���֡�
   * x(t) = [����Ҥγ��١���(t),
   *         ������γ��١���(t),  % Ploplys�ȤϵդʤΤ���ա�
   *         ����Ҥγ�®�٦�'(t), %   �ץ����Ū�Թ��y[k][0,1]=x[k][0,1]��
   *         ������γ�®�٦�'(t)];
   */
  
  const double tmp=1.0/(a*b-c*c);
  
  A = (double[4][4]){{       0, 0, 1, 0},
		     {       0, 0, 0, 1},
		     { b*d*tmp, 0, 0, 0},
		     {-c*d*tmp, 0, 0, 0}};
  B = (double[4]){0, 0, -c*tmp, a*tmp};
}


inline
void pendulum(dMatrix<4,4> &A, dMatrix<4,1> &B) {
  dVector<4> tmp;
  
  pendulum(A, tmp);
  for (int i=0; i<4; i++) B(i,0) = tmp(i);
}

#endif
