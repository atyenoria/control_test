// Time-stamp: <2012-12-26 16:59:57 kobayasi>

#ifndef __PENDULUM_H__
#define __PENDULUM_H__

#define PENDULUM_MODEL    "N.J.Ploplys"

#define PENDULUM_ANGLE    0 //振り子の角度を表すインデックス
#define ARM_ANGLE         1 //アームの角度　　〃
#define PENDULUM_VELOCITY 2 //振り子の角速度　〃
#define ARM_VELOCITY      3 //アームの角速度　〃


inline
void pendulum(dMatrix<4,4> &A, dVector<4> &B) {
  /*
   * 回転型倒立振子のモデル
   * N. J. Ploplys, "Wireless Feedback Control of Mechanical Systems"
   * Master Thesis, Department of Mechanical Engineering,
   *   University of Illinois, Champaign, IL, 2003
   */
  
  const double M=0.004;  //振り子のおもりの質量[kg]
  const double m=0.025;  //振り上げ棒の質量[kg]
  const double l=0.241;  //振り上げ棒の長さ[m]
  const double r=0.152;  //アームの長さ[m]
  const double J=1.21e-3;//アームの回転軸周りの慣性モーメント[kgm^2]
  const double g=9.81;   //重力加速度[m/s^2]
  
  const double a=(M+m/3)*l*l;
  const double b=J+(M+m)*r*r;
  const double c=(M+m/2)*r*l;
  const double d=(M+m/2)*g*l;
  
  /*
   * 連続時間システムの微分方程式（線形近似）
   * a*θ''(t) + c*φ''(t) -d*θ'(t) = 0
   * c*θ''(t) + b*φ''(t) = u(t)
   *
   * 状態空間モデル
   * x'(t) = A*x(t) + B*u(t)
   * y(t)  = C*x(t) + D*u(t)
   *
   * 入力：
   * u(t) =  アームの回転トルク
   * 状態：
   * x(t) = [振り子の角度　θ(t),
   *         アームの角度　φ(t),  % Ploplysとは逆なので注意：
   *         振り子の角速度θ'(t), %   プログラム的都合（y[k][0,1]=x[k][0,1]）
   *         アームの角速度φ'(t)];
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
