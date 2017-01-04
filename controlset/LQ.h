// Time-stamp: <2013-04-26 13:31:44 kobayasi>

// 最適フィルタ（カルマンフィルタ）とセットで使う
 
#ifndef __LQ_H__
#define __LQ_H__


// no packetloss (=LQR)
/*
template <class E, class Eb, class Eq>
inline
dVector<Eb::Size> dlqr(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
		       const dMatrixExpression<Eq> &Q, const double R) {
  dVector<Eb::Size> G;
  dMatrix<E::Rows,E::Cols> K, Kn=Q;
  
  G = B*Kn*A/(B*Kn*B+R);
  do {
    K = Kn;
    Kn = trans(A)*K*(A - matmul(B,G)) + Q;
    
    G = B*Kn*A/(B*Kn*B+R);
  } while (fnorm(Kn-K) > DBL_EPSILON*fnorm(K));
  
  return G;
}
*/

// TCP (packetloss rate of u is a)
template <class E, class Eb, class Eq>
inline
dVector<Eb::Size> dlqr_tcp(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
			   const dMatrixExpression<Eq> &Q, const double R,
			   const double a) {
  // L. Schenato, et.al.,
  // "Foundation of Control and Estimation over Lossy Networks"
  
  dVector<Eb::Size> G;
  dMatrix<E::Rows,E::Cols> K, Kn=Q;
  
  G = B*Kn*A/(B*Kn*B+R);
  do {
    K = Kn;
    Kn = trans(A)*K*(A - (1-a)*matmul(B,G)) + Q;
    
    G = B*Kn*A/(B*Kn*B+R);
  } while (fnorm(Kn-K) > DBL_EPSILON*fnorm(K));
  
  return G;
}


// UDP (packetloss of u is a, packetloss rate of x is b)
// C=I(すべて観測可能)かつV=0(観測誤差なし)のときのみ最適
template <class E, class Eb, class Eq>
inline
dVector<Eb::Size> dlqr_udp(const dMatrixExpression<E> &A, const dVectorExpression<Eb> &B,
			   const dMatrixExpression<Eq> &Q, const double R,
			   const double a, const double b) {
  // L. Schenato, et.al.,
  // "Foundation of Control and Estimation over Lossy Networks"
  
  dVector<Eb::Size> G;
  dMatrix<E::Rows,E::Cols> S, Sn=Q;
  dMatrix<E::Rows,E::Cols> T, Tn=0.0;
  
  const double c = a*b;
  
  G = B*Sn*A/(B*((1-c)*Sn + c*Tn)*B+R);
  do {
    S = Sn;
    T = Tn;
    
    Sn = trans(A)*S*(A - (1-a)*matmul(B,G)) + Q;
    Tn = trans(A)*(b*T + (1-b)*S)*A + Q;
    
    G = B*Sn*A/(B*((1-c)*Sn + c*Tn)*B+R);
  } while (fnorm(Sn-S) > DBL_EPSILON*fnorm(S));
  
  return G;
}


#endif
