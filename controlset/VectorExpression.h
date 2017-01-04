/* Time-stamp: <2013-02-26 18:20:20 kobayasi>
   
   ベクトルの高速演算のためのExpressionTemplate（VectorClass.hから利用）
   ver.2013-02-20  by K.Kobayashi
   
   fabs(B) ：各要素の絶対値をとったベクトルを返す
   round(B)：　〃　　四捨五入　〃
   ceil(B) ：　〃　　天井関数　〃
   floor(B)：　〃　　床関数　　〃   
   square(B)： 〃　　二乗　　　〃   
   
   subvector(B,i,N)：ベクトルの(i)から部分ベクトル<N>を返す
                     （Bを直接書き換え可能：subvector(B,...)=X）
		     
   norm(B)   ：ユークリッドノルム(L2ノルム)を返す（= l2norm(B)）
   l1norm(B) ：L1ノルムを返す
   l0norm(B) ：L0ノルムを返す
   infnorm(B)：無限大ノルムを返す
   max(B)    ：要素の最大値を返す
   min(B)    ：要素の最小値を返す
   sum(B)    ：要素の和を返す
   
   print(B)：　　 printfを使ってベクトルを表示する
   printf(fmt,B)：↑の値の表示フォーマット指定版
*/


// ベクトルを表すクラス
template<class E>
class dVectorExpression {
public:
  const E &operator() () const {
    return *static_cast<const E *>(this);
  }
  E &operator() () {
    return *static_cast<E *>(this);
  }
};


// 演算処理のための構造体（行列の演算でも利用）
struct OperationAdd {
  template<class T>
  static T apply(T ll, T rr) {
    return ll + rr;
  }
};
struct OperationSub {
  template<class T>
  static T apply(T ll, T rr) {
    return ll - rr;
  }
};
struct OperationMul {
  template<class T>
  static T apply(T ll, T rr) {
    return ll * rr;
  }
};
struct OperationDiv {
  template<class T>
  static T apply(T ll, T rr) {
    return ll / rr;
  }
};
struct OperationNeg {
  template<class T>
  static T apply(T rr) {
    return -rr;
  }
};
struct OperationFabs {
  template<class T>
  static T apply(T rr) {
    return fabs(rr);
  }
};
struct OperationRound {
  template<class T>
  static T apply(T rr) {
    return round(rr);
  }
};
struct OperationCeil {
  template<class T>
  static T apply(T rr) {
    return ceil(rr);
  }
};
struct OperationFloor {
  template<class T>
  static T apply(T rr) {
    return floor(rr);
  }
};
struct OperationSquare {
  template<class T>
  static T apply(T rr) {
    return rr*rr;
  }
};


// 演算コストに応じて演算結果をキャッシュする構造体
// （ベクトル*行列のように複数要素に渡る演算があるときにキャッシュ：Costs=1）
template<class R, bool Costs_>
class dVectorCache:
  public dVectorExpression< dVectorCache<R, Costs_> > {
public:
  static const bool Costs = Costs_;
  static const int Size = R::Size;
private:
  const R &rv_; //no cache
public:
  dVectorCache(const R &rv): rv_(rv) {}
  
  double operator[] (int i) const { return rv_[i]; }
  double operator() (int i) const { return rv_(i); }
};
template<class R>
class dVectorCache<R, 1>:
  public dVectorExpression< dVectorCache<R, 1> > {
public:
  static const bool Costs = 0;
  static const int Size = R::Size;
private:
  const dVector<R::Size> rv_; //cache
public:
  dVectorCache(const R &rv): rv_(rv) {}
  
  double operator[] (int i) const { return rv_[i]; }
  double operator() (int i) const { return rv_(i); }
};


// ベクトルの演算式を表すクラス（ベクトルとベクトルの加減算）
template<class L, class Ope, class R>
class dVectorOperation:
  public dVectorExpression< dVectorOperation<L, Ope, R> > {
public:
  static const bool Costs = L::Costs | R::Costs;
  static const int Size = L::Size;
private:
  const L &lv_;
  const R &rv_;
public:
  dVectorOperation(const L &lv, const R &rv): lv_(lv), rv_(rv) {}
  
  double operator[] (int i) const { return Ope::apply(lv_[i], rv_[i]); }
  double operator() (int i) const { return Ope::apply(lv_(i), rv_(i)); }
};

// ベクトルの演算式を表すクラス（ベクトルと数値の乗除算）
template<class L, class Ope>
class dVectorOperation<L, Ope, double>:
  public dVectorExpression< dVectorOperation<L, Ope, double> > {
public:
  static const bool Costs = L::Costs;
  static const int Size = L::Size;
private:
  const L &lv_;
  const double rd_;//  __attribute__((aligned(8))); //速くなる？
public:
  dVectorOperation(const L &lv, const double rd): lv_(lv), rd_(rd) {}
  
  double operator[] (int i) const { return Ope::apply(lv_[i], rd_); }
  double operator() (int i) const { return Ope::apply(lv_(i), rd_); }
};

// ベクトルの演算式を表すクラス（ベクトルの負演算など）
template<class Ope, class R>
class dVectorOperationS:
  public dVectorExpression< dVectorOperationS<Ope, R> > {
public:
  static const bool Costs = R::Costs;
  static const int Size = R::Size;
private:
  const R &rv_;
public:
  dVectorOperationS(const R &rv): rv_(rv) {}
  
  double operator[] (int i) const { return Ope::apply(rv_[i]); }
  double operator() (int i) const { return Ope::apply(rv_(i)); }
};

// 加算（ベクトル + ベクトル）
template<class L, class R>
inline
dVectorOperation<L, OperationAdd, R> operator + (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  return dVectorOperation<L, OperationAdd, R>( lv(), rv() );
}
// 減算（ベクトル - ベクトル）
template<class L, class R>
inline
dVectorOperation<L, OperationSub, R> operator - (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  return dVectorOperation<L, OperationSub, R>( lv(), rv() );
}
// 乗算（ベクトル * 数値）
template<class L>
inline
dVectorOperation<L, OperationMul, double> operator * (const dVectorExpression<L> &lv, const double rd) {
  return dVectorOperation<L, OperationMul, double>( lv(), rd );
}
// 乗算（数値 * ベクトル）＝（ベクトル * 数値）
template<class R>
inline
dVectorOperation<R, OperationMul, double> operator * (const double ld, const dVectorExpression<R> &rv) {
  return dVectorOperation<R, OperationMul, double>( rv(), ld );
}
// 除算（ベクトル / 数値）
template<class L>
inline
dVectorOperation<L, OperationDiv, double> operator / (const dVectorExpression<L> &lv, const double rd) {
  return dVectorOperation<L, OperationDiv, double>( lv(), rd );
}
// 負（- ベクトル）
template<class R>
inline
dVectorOperationS<OperationNeg, R> operator - (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationNeg, R>( rv() );
}
// 等号（ベクトル == ベクトル）// 要素はdouble型なので演算結果を比較するのはおすすめしない
template<class L, class R>
inline
bool operator == (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  for (int i = 0; i < L::Size; i++) {
    if ( lv()[i] != rv()[i] ) return 0;
  }
  return 1;
}
// 不等号（ベクトル != ベクトル）
template<class L, class R>
inline
bool operator != (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  for (int i = 0; i < L::Size; i++) {
    if ( lv()[i] != rv()[i] ) return 1;
  }
  return 0;
}

// 絶対 fabs(ベクトル)
template<class R>
inline
dVectorOperationS<OperationFabs, R> fabs (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationFabs, R>( rv() );
}
// 四捨五入 round(ベクトル)
template<class R>
inline
dVectorOperationS<OperationRound, R> round (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationRound, R>( rv() );
}
// 天井関数 ceil(ベクトル)
template<class R>
inline
dVectorOperationS<OperationCeil, R> ceil (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationCeil, R>( rv() );
}
// 床関数 floor(ベクトル)
template<class R>
inline
dVectorOperationS<OperationFloor, R> floor (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationFloor, R>( rv() );
}
// 各要素の二乗 square(ベクトル)
template<class R>
inline
dVectorOperationS<OperationSquare, R> square (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationSquare, R>( rv() );
}

// 内積（ベクトル * ベクトル）
template<class L, class R>
inline
double operator * (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  double result = 0;
  for (int i = 0; i < L::Size; i++) {
    result += lv()[i] * rv()[i];
  }
  return result;
}


/**************************************************************************/


// ユークリッドノルム（L2ノルム） norm(ベクトル)
template<class R>
inline
double norm (const dVectorExpression<R> &rv) {
  return sqrt( rv * rv );
}

#define l2norm(rv)  norm(rv)

// L1ノルム l1norm(ベクトル)
template<class R>
inline
double l1norm (const dVectorExpression<R> &rv) {
  return sum( fabs(rv) );
}

// L0ノルム l0norm(ベクトル)
template<class R>
inline
double l0norm (const dVectorExpression<R> &rv) {
  double result = 0;
  double tmp;
  for (int i=0; i<R::Size; i++) {
    tmp = rv()[i];
    if (tmp != 0) result++;
  }
  return result;
}

// 無限大ノルム infnorm(ベクトル)
template<class R>
inline
double infnorm (const dVectorExpression<R> &rv) {
  return max( fabs(rv) );
}

// 最大値 max(ベクトル)
template<class R>
inline
double max (const dVectorExpression<R> &rv) {
  double result = -DBL_MAX;
  double tmp;
  for (int i=0; i<R::Size; i++) {
    tmp = rv()[i];
    if (result < tmp) result = tmp;
  }
  return result;
}

// 最小値 min(ベクトル)
template<class R>
inline
double min (const dVectorExpression<R> &rv) {
  double result = DBL_MAX;
  double tmp;
  for (int i=0; i<R::Size; i++) {
    tmp = rv()[i];
    if (result > tmp) result = tmp;
  }
  return result;
}

// 要素の和 sum(ベクトル)
template<class R>
inline
double sum (const dVectorExpression<R> &rv) {
  double result = 0;
  for (int i=0; i<R::Size; i++) {
    result += rv()[i];
  }
  return result;
}


/**************************************************************************/

// 表示用
template<class E>
void printf(const char *fmt, const dVectorExpression<E> &dv) {
  printf("{");
  for (int i = 0; i < E::Size-1; i++) {
    printf(fmt, dv()[i]);
    printf(", ");
  }
  printf(fmt, dv()[E::Size-1]);
  printf("}\n");
}

// 表示用
template<class E>
inline
void print(const dVectorExpression<E> &dv) {
  printf("% e", dv);
}
