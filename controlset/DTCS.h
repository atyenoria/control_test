// Time-stamp: <2013-04-26 13:31:53 kobayasi>

/*
 * [離散化]
 * dscr(A,B,C,D,t)
 * dscr(A,B,t)        ←C,Dは変わらないので簡略Ver
 * dscr(A,B,C,D,W,V,t)
 * dscr(A,B,W,t)      ←C,D,Vは変わらないので簡略版Ver
 *     (in/out) A,B,C,D,W,V：連続時間状態空間モデル(in)
 *                           離散時間状態空間モデル(out)
 *     (in) t：サンプリング間隔
 *     返り値：なし
 *
 * [LQRレギュレータ]
 * dlqr(Ad,Bd,Q,R)
 *     (in) Ad,Bd：離散時間状態空間モデル
 *     (in) Q,R：重み
 *     返り値：離散時間コントローラゲイン
 */
 
#ifndef __DTCS_H__
#define __DTCS_H__

#include "MatrixClass.h"

/*
 * 連続時間状態空間モデル
 * x'(t) = A*x(t) + B*u(t) + w(t)
 * y(t)  = C*x(t) + D*u(t) + v(t)
 *
 * w(t)：プロセス雑音
 *       E{w(t)}=0,  E{w(t)*w(t)^T}=W（共分散行列）
 * v(t)：観測雑音
 *       E{v(t)}=0,  E{v(t)*v(t)^T}=V（共分散行列）
 *
 *       E{w(t)*v(t)^T} = 0
 */
/*
 * 離散時間状態空間モデル（Ts：サンプリング間隔）
 * x[k+1] = Ad*x[k] + Bd*u[k] + w[k]
 * y[k]   = Cd*x[k] + Dd*u[k] + v[k]
 *
 * w[k]：離散時間プロセス雑音
 *       E{w[k]}=0,  E{w[k]*w[k]^T}=Wd（共分散行列）
 * v[k]：離散時間観測雑音
 *       E{v[k]}=0,  E{v[k]*v[k]^T}=Vd（共分散行列）
 */
/*
 * 連続時間システムの離散化
 *                    Ts
 * Ad=exp(A*Ts),  Bd=∫exp(A*t)dt *B,  Cd=C,  Dd=D
 *                    0
 *     Ts
 * Wd=∫exp(A*t)*W*exp(A^T*t)dt,  Vd=V
 *     0
 */

/*
 * 離散時間代数リカッチ方程式
 * DARE(Discrete-time Algebraic Riccati Equation)
 *
 * A^T X A - (A^T X B)(R + B^T X B)^-1(B^T X A) + Q = X
 *
 * X：リカッチ方程式の解
 */

/*
 * 線形二次状態フィードバックレギュレーター
 * LQR(Linear Quadratic Regulator)
 *
 * 最適レギュレータ：
 * 　コスト関数Jを最小化する状態フィードバック入力のゲインKd
 * u[k] = -Kd*x[k]
 *
 * 離散時間モデルのコスト関数
 *     ∞
 * J = Σ (x[k]^T*Q*x[k] + u[k]^T*R*u[k])
 *     k=0
 * Q, R：重み（それぞれ半正定値，正定値）
 *       Q大→早い応答  
 *       R大→小さい入力
 */

/*
 * 同一次元オブザーバ
 * xe[k+1] = Ad*xe[k] + Bd*ue[k] + Ld*(y[k]-ye[k])
 * ye[k]   = Cd*xe[k] + Dd*ue[k]
 *
 * Ld：オブザーバゲイン
 *
 * 最適オブザーバ：コスト関数Jを最小化するオブザーバゲインLd
 *
 * 離散時間カルマンフィルタのコスト関数
 * J = lim E{ (x[k]-xe[k])*(x[k]-xe[k])^T }
 *    k→∞
 *
 * LQRとの双対性：
 * Kd → Ld^T,  Ad → Ad^T,  Bd → Cd^T,  Q → Wd,  R → Vd
 */

template <int N, int L>
inline
void dscr(dMatrix<N,N> &A, dVector<N> &B,
	  dMatrix<L,N> &C, dVector<L> &D,
	  const double dt) {
  // 離散化
  // C.F.Van Loan, "Computing Integrals Involving the Matrix Expopential"
  
  // C,Dは離散化しても同じ
  
  dMatrix<N+1,N+1> H, eH;
  
  // H = [ A  B ]
  //     [ 0  0 ]
  H.clear();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      H(i,j) = A(i,j);
    }
    H(i,N) = B(i);
  }
  
  expm(H*dt, eH);
  // expm(H*dt) = [ F3 G3 ]
  //              [  0 F4 ]
  // Ad = F3
  // Bd = G3
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      A(i,j) = eH(i,j);
    }
    B(i) = eH(i,N);
  }
}

template <int N, int L>
inline
void dscr(dMatrix<N,N> &A, dVector<N> &B,
	  dMatrix<L,N> &C, dVector<L> &D,
	  dMatrix<N,N> &W, dMatrix<L,L> &V,
	  const double dt) {
  // 離散化
  // C.F.Van Loan, "Computing Integrals Involving the Matrix Expopential"
  
  // C,D,Vは離散化しても同じ
  
  dMatrix<N+N+1,N+N+1> H, eH;
  dMatrix<N,N> tmp;
  
  //     [-A^t Q  0]
  // H = [ 0   A  B]
  //     [ 0   0  0]
  H.clear();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      H(i,j) = -A(j,i);
      H(i,j+N) = W(i,j);
      H(i+N,j+N) = A(i,j);
    }
    H(i+N,N+N) = B(i);
  }
  
  expm(H*dt, eH);
  //              [ F2 G2 H2 ]
  // expm(H*dt) = [ 0  F3 G3 ]
  //              [ 0   0 F4 ]
  // Ad = F3
  // Bd = G3
  // Wd = F3^t * G2
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      tmp(i,j) = eH(i,j+N);
      A(i,j) = eH(i+N,j+N);
    }
    B(i) = eH(i+N,N+N);
  }
  W = trans(A) * tmp;
}

template <int N>
inline
void dscr(dMatrix<N,N> &A, dVector<N> &B,
	  const double dt) {
  // 離散化（簡略化Ver）
  dscr(A,B,A,B,dt); //2つ目のA,Bはダミー
}

template <int N>
inline
void dscr(dMatrix<N,N> &A, dVector<N> &B,
	  dMatrix<N,N> &W,
	  const double dt) {
  // 離散化（簡略化Ver）
  dscr(A,B,A,B,W,W,dt); //2つ目のA,B,Wはダミー
}


template <class E, class Eb, class Eq>
inline
dMatrix<E::Rows,E::Rows> dare(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
			      const dMatrixExpression<Eq> &Q, const double R) {
  // 離散時間代数リカッチ方程式
  // DARE(Discrete-time Algebraic Riccati Equation)
  
  // Structure-preserving Doubling Algorithm (SDA)
  // E. K.-W. Chu, "Structure-preserving algorithms for periodic discrete-time algebraic Riccati equations"
  
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  dMatrix<E::Rows,E::Rows> Ak, Ak1=A;
  dMatrix<E::Rows,E::Rows> Gk, Gk1=matmul(B,B)/R; //B R^-1 B^T
  dMatrix<E::Rows,E::Rows> Hk, Hk1=Q;
  dMatrix<E::Rows,E::Rows> W, V;
  do {
    Ak = Ak1;
    Gk = Gk1;
    Hk = Hk1;
    inv(I + Gk*Hk, W);
    
    V = W * Ak;
    Ak1 = Ak * V;
    Hk1 = Hk + trans(V) * Hk * Ak;
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //転置をうまく使うことでV2の計算を簡略化
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  return Hk1;
}

template <class E, class Eb, class Eq, class Er>
inline
dMatrix<E::Rows,E::Rows> dare(const dMatrixExpression<E> &A, const dMatrixExpression<Eb> &B,
			      const dMatrixExpression<Eq> &Q, const dMatrixExpression<Er> &R) {
  // 離散時間代数リカッチ方程式
  // DARE(Discrete-time Algebraic Riccati Equation)
  
  // Structure-preserving Doubling Algorithm (SDA)
  // E. K.-W. Chu, "Structure-preserving algorithms for periodic discrete-time algebraic Riccati equations"
  
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  dMatrix<Eb::Cols,Eb::Rows> tmp;
  solve(R,tmp = trans(B)); //B R^-1 B^T の計算用(tmp = R^-1 B^T)
  
  dMatrix<E::Rows,E::Rows> Ak, Ak1=A;
  dMatrix<E::Rows,E::Rows> Gk, Gk1=B*tmp; 
  dMatrix<E::Rows,E::Rows> Hk, Hk1=Q;
  dMatrix<E::Rows,E::Rows> W, V;
  do {
    Ak = Ak1;
    Gk = Gk1;
    Hk = Hk1;
    inv(I + Gk*Hk, W);
    
    V = W * Ak;
    Ak1 = Ak * V;
    Hk1 = Hk + trans(V) * Hk * Ak;
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //転置をうまく使うことでV2の計算を簡略化
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  return Hk1;
}

template <class E, class Eb, class Eq>
inline
dVector<Eb::Size> dlqr(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
		       const dMatrixExpression<Eq> &Q, const double R) {
  // 線形二次状態フィードバックレギュレーター
  // LQR(Linear Quadratic Regulator)
  
  // 最後のゲイン計算以外はdare()と同じ（Hk1を無駄にしないため書く）
  int i=0;
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  dMatrix<E::Rows,E::Rows> Ak, Ak1=A;
  dMatrix<E::Rows,E::Rows> Gk, Gk1=matmul(B,B)/R; //B R^-1 B^T
  dMatrix<E::Rows,E::Rows> Hk, Hk1=Q;
  dMatrix<E::Rows,E::Rows> W, V;
  do {
    Ak = Ak1;
    Gk = Gk1;
    Hk = Hk1;
    inv(I + Gk*Hk, W);
    
    V = W * Ak;
    Ak1 = Ak * V;
    Hk1 = Hk + trans(V) * Hk * Ak;
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //転置をうまく使うことでV2の計算を簡略化
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  return B*Hk1*A/(B*Hk1*B+R);
}

template <class E, class Eb, class Eq, class Er>
inline
dMatrix<Eb::Cols,Eb::Rows> dlqr(const dMatrixExpression<E> &A, const dMatrixExpression<Eb> &B,
				const dMatrixExpression<Eq> &Q, const dMatrixExpression<Er> &R) {
  // 線形二次状態フィードバックレギュレーター
  // LQR(Linear Quadratic Regulator)
  
  // 最後のゲイン計算以外はdare()と同じ（Hk1を無駄にしないため書く）
  
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  dMatrix<Eb::Cols,Eb::Rows> tmp;
  solve(R,tmp = trans(B)); //B R^-1 B^Tの初期値計算用(tmp = R^-1 B^T)
  
  dMatrix<E::Rows,E::Rows> Ak, Ak1=A;
  dMatrix<E::Rows,E::Rows> Gk, Gk1=B*tmp; 
  dMatrix<E::Rows,E::Rows> Hk, Hk1=Q;
  dMatrix<E::Rows,E::Rows> W, V;
  do {
    Ak = Ak1;
    Gk = Gk1;
    Hk = Hk1;
    inv(I + Gk*Hk, W);
    
    V = W * Ak;
    Ak1 = Ak * V;
    Hk1 = Hk + trans(V) * Hk * Ak;
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //転置をうまく使うことでV2の計算を簡略化
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  solve(trans(B)*Hk1*B+R,tmp = trans(B)*Hk1*A);
  return tmp;
}


#endif
