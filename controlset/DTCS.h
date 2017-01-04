// Time-stamp: <2013-04-26 13:31:53 kobayasi>

/*
 * [Υ����]
 * dscr(A,B,C,D,t)
 * dscr(A,B,t)        ��C,D���Ѥ��ʤ��ΤǴ�άVer
 * dscr(A,B,C,D,W,V,t)
 * dscr(A,B,W,t)      ��C,D,V���Ѥ��ʤ��ΤǴ�ά��Ver
 *     (in/out) A,B,C,D,W,V��Ϣ³���־��ֶ��֥�ǥ�(in)
 *                           Υ�����־��ֶ��֥�ǥ�(out)
 *     (in) t������ץ�󥰴ֳ�
 *     �֤��͡��ʤ�
 *
 * [LQR�쥮��졼��]
 * dlqr(Ad,Bd,Q,R)
 *     (in) Ad,Bd��Υ�����־��ֶ��֥�ǥ�
 *     (in) Q,R���Ť�
 *     �֤��͡�Υ�����֥���ȥ��饲����
 */
 
#ifndef __DTCS_H__
#define __DTCS_H__

#include "MatrixClass.h"

/*
 * Ϣ³���־��ֶ��֥�ǥ�
 * x'(t) = A*x(t) + B*u(t) + w(t)
 * y(t)  = C*x(t) + D*u(t) + v(t)
 *
 * w(t)���ץ�������
 *       E{w(t)}=0,  E{w(t)*w(t)^T}=W�ʶ�ʬ�������
 * v(t)����¬����
 *       E{v(t)}=0,  E{v(t)*v(t)^T}=V�ʶ�ʬ�������
 *
 *       E{w(t)*v(t)^T} = 0
 */
/*
 * Υ�����־��ֶ��֥�ǥ��Ts������ץ�󥰴ֳ֡�
 * x[k+1] = Ad*x[k] + Bd*u[k] + w[k]
 * y[k]   = Cd*x[k] + Dd*u[k] + v[k]
 *
 * w[k]��Υ�����֥ץ�������
 *       E{w[k]}=0,  E{w[k]*w[k]^T}=Wd�ʶ�ʬ�������
 * v[k]��Υ�����ִ�¬����
 *       E{v[k]}=0,  E{v[k]*v[k]^T}=Vd�ʶ�ʬ�������
 */
/*
 * Ϣ³���֥����ƥ��Υ����
 *                    Ts
 * Ad=exp(A*Ts),  Bd=��exp(A*t)dt *B,  Cd=C,  Dd=D
 *                    0
 *     Ts
 * Wd=��exp(A*t)*W*exp(A^T*t)dt,  Vd=V
 *     0
 */

/*
 * Υ����������ꥫ�å�������
 * DARE(Discrete-time Algebraic Riccati Equation)
 *
 * A^T X A - (A^T X B)(R + B^T X B)^-1(B^T X A) + Q = X
 *
 * X���ꥫ�å��������β�
 */

/*
 * �����󼡾��֥ե����ɥХå��쥮��졼����
 * LQR(Linear Quadratic Regulator)
 *
 * ��Ŭ�쥮��졼����
 * �������ȴؿ�J��Ǿ���������֥ե����ɥХå����ϤΥ�����Kd
 * u[k] = -Kd*x[k]
 *
 * Υ�����֥�ǥ�Υ����ȴؿ�
 *     ��
 * J = �� (x[k]^T*Q*x[k] + u[k]^T*R*u[k])
 *     k=0
 * Q, R���Ťߡʤ��줾��Ⱦ�����͡������͡�
 *       Q�碪�ᤤ����  
 *       R�碪����������
 */

/*
 * Ʊ�켡�����֥�����
 * xe[k+1] = Ad*xe[k] + Bd*ue[k] + Ld*(y[k]-ye[k])
 * ye[k]   = Cd*xe[k] + Dd*ue[k]
 *
 * Ld�����֥����Х�����
 *
 * ��Ŭ���֥����С������ȴؿ�J��Ǿ������륪�֥����Х�����Ld
 *
 * Υ�����֥���ޥ�ե��륿�Υ����ȴؿ�
 * J = lim E{ (x[k]-xe[k])*(x[k]-xe[k])^T }
 *    k����
 *
 * LQR�Ȥ���������
 * Kd �� Ld^T,  Ad �� Ad^T,  Bd �� Cd^T,  Q �� Wd,  R �� Vd
 */

template <int N, int L>
inline
void dscr(dMatrix<N,N> &A, dVector<N> &B,
	  dMatrix<L,N> &C, dVector<L> &D,
	  const double dt) {
  // Υ����
  // C.F.Van Loan, "Computing Integrals Involving the Matrix Expopential"
  
  // C,D��Υ�������Ƥ�Ʊ��
  
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
  // Υ����
  // C.F.Van Loan, "Computing Integrals Involving the Matrix Expopential"
  
  // C,D,V��Υ�������Ƥ�Ʊ��
  
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
  // Υ�����ʴ�ά��Ver��
  dscr(A,B,A,B,dt); //2���ܤ�A,B�ϥ��ߡ�
}

template <int N>
inline
void dscr(dMatrix<N,N> &A, dVector<N> &B,
	  dMatrix<N,N> &W,
	  const double dt) {
  // Υ�����ʴ�ά��Ver��
  dscr(A,B,A,B,W,W,dt); //2���ܤ�A,B,W�ϥ��ߡ�
}


template <class E, class Eb, class Eq>
inline
dMatrix<E::Rows,E::Rows> dare(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
			      const dMatrixExpression<Eq> &Q, const double R) {
  // Υ����������ꥫ�å�������
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
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //ž�֤򤦤ޤ��Ȥ����Ȥ�V2�η׻����ά��
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  return Hk1;
}

template <class E, class Eb, class Eq, class Er>
inline
dMatrix<E::Rows,E::Rows> dare(const dMatrixExpression<E> &A, const dMatrixExpression<Eb> &B,
			      const dMatrixExpression<Eq> &Q, const dMatrixExpression<Er> &R) {
  // Υ����������ꥫ�å�������
  // DARE(Discrete-time Algebraic Riccati Equation)
  
  // Structure-preserving Doubling Algorithm (SDA)
  // E. K.-W. Chu, "Structure-preserving algorithms for periodic discrete-time algebraic Riccati equations"
  
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  dMatrix<Eb::Cols,Eb::Rows> tmp;
  solve(R,tmp = trans(B)); //B R^-1 B^T �η׻���(tmp = R^-1 B^T)
  
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
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //ž�֤򤦤ޤ��Ȥ����Ȥ�V2�η׻����ά��
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  return Hk1;
}

template <class E, class Eb, class Eq>
inline
dVector<Eb::Size> dlqr(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
		       const dMatrixExpression<Eq> &Q, const double R) {
  // �����󼡾��֥ե����ɥХå��쥮��졼����
  // LQR(Linear Quadratic Regulator)
  
  // �Ǹ�Υ�����׻��ʳ���dare()��Ʊ����Hk1��̵�̤ˤ��ʤ�����񤯡�
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
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //ž�֤򤦤ޤ��Ȥ����Ȥ�V2�η׻����ά��
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  return B*Hk1*A/(B*Hk1*B+R);
}

template <class E, class Eb, class Eq, class Er>
inline
dMatrix<Eb::Cols,Eb::Rows> dlqr(const dMatrixExpression<E> &A, const dMatrixExpression<Eb> &B,
				const dMatrixExpression<Eq> &Q, const dMatrixExpression<Er> &R) {
  // �����󼡾��֥ե����ɥХå��쥮��졼����
  // LQR(Linear Quadratic Regulator)
  
  // �Ǹ�Υ�����׻��ʳ���dare()��Ʊ����Hk1��̵�̤ˤ��ʤ�����񤯡�
  
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  dMatrix<Eb::Cols,Eb::Rows> tmp;
  solve(R,tmp = trans(B)); //B R^-1 B^T�ν���ͷ׻���(tmp = R^-1 B^T)
  
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
    Gk1 = Gk + Ak * Gk * trans(Ak * W); //ž�֤򤦤ޤ��Ȥ����Ȥ�V2�η׻����ά��
  } while (fnorm(Hk1-Hk) > DBL_EPSILON*fnorm(Hk));
  
  // Hk1 is a solution of the discrite-time Riccati eq.
  solve(trans(B)*Hk1*B+R,tmp = trans(B)*Hk1*A);
  return tmp;
}


#endif
