// Time-stamp: <2012-12-26 17:00:16 kobayasi>

#ifndef __PENDULUM_H__
#define __PENDULUM_H__

#define PENDULUM_MODEL    "REALTECH RTC05"

#define PENDULUM_ANGLE    0 //振り子の角度を表すインデックス
#define ARM_ANGLE         1 //アームの角度　　〃
#define PENDULUM_VELOCITY 2 //振り子の角速度　〃
#define ARM_VELOCITY      3 //アームの角速度　〃


inline
void pendulum(dMatrix<4,4> &A, dVector<4> &B) {
  /*
   * 回転型倒立振子のモデル
   * REALTECH RTC05
   */
  
  // 振子(40cm)
  //const double m=0.054; //振り上げ棒の質量[kg]
  //const double l=0.20;  //振り上げ棒の長さの1/2[m]
  
  // 振子(20cm)
  const double m=0.016; //振り上げ棒の質量[kg]
  //const double l=0.10;  //振り上げ棒の長さの1/2[m]
  const double l=0.085;
  //const double J=0.0048; //アームの回転軸周りの慣性モーメント[kgm^2]
  const double J=0;
  //const double r=0.20;   //アームの長さ[m]
  const double r=0.206;
  //const double g=9.81;   //重力加速度[m/s^2]
  const double g=9.797;

  /*ß
  // 操作量がアームの回転トルクの場合
  const double a = -m*r*g/J;
  const double b = 0;
  const double c = g*(J+m*r*r)/(l*J);
  const double d = 0;
  const double e = 1/J;
  const double f = -r/(l*J);
  */
  
  const double R=8.3;    //モータトルク直流抵抗[Ω]
  const double Km=0.023; //モータトルク定数[Nm/A,Vs/rad]
  const double Kg=7.5;   //ギア比 120:16
  
  // 操作量がアームモータの入力電圧の場合
  const double a = -m*r*g/J;
  const double b = -(Kg*Km)*(Kg*Km)/(R*J);
  const double c = g*(J+m*r*r)/(l*J);
  const double d = r*(Kg*Km)*(Kg*Km)/(R*l*J);
  const double e = Kg*Km/(R*J);
  const double f = -Kg*Km*r/(R*l*J);
  
  
  /*
   * 状態空間モデル
   * x'(t) = A*x(t) + B*u(t)
   * y(t)  = C*x(t) + D*u(t)
   *
   * 入力：
   * u(t) =  アームモータの入力電圧（or アームの回転トルク）
   * 状態：
   * x(t) = [振り子の角度　θ(t),
   *         アームの角度　φ(t),  % 仕様書とは順序が異なるので注意：
   *         振り子の角速度θ'(t), %   プログラム的都合（y[k][0,1]=x[k][0,1]）
   *         アームの角速度φ'(t)];
   */
  
  A = (double[4][4]){{ 0, 0, 1, 0},
		     { 0, 0, 0, 1},
		     { c, 0, 0, d},
		     { a, 0, 0, b}};
  B = (double[4]){0, 0, f, e};
  
  
  /* なお，仕様書では，
   * x(t) = [アームの角度　θ(t),
   *         振子の角度　　α(t),
   *         アームの角速度θ'(t),
   *         振子の角速度　α'(t)];
   
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
