/* Time-stamp: <2013-04-26 13:38:54 kobayasi>
 
   行列の高速演算のためのExpressionTemplate（MatrixClass.hから利用）
   ver.2013-04-26  by K.Kobayashi
   
   fabs(A) ：各要素の絶対値をとった行列を返す
   round(A)：　〃　　四捨五入　〃
   ceil(A) ：　〃　　天井関数　〃
   floor(A)：　〃　　床関数　　〃
   square(A)： 〃　　二乗　　　〃   （←行列の二乗ではないので注意）
   
   submatrix(A,i,j,N,M)：行列Aの(i,j)から部分行列<N,M>を返す
                        （Aを直接書き換え可能：submatrix(A,...)=X）
   
   rowvector(A,i)：i行ベクトルを返す
   colvector(A,j)：j列ベクトルを返す
   
   diag(A)  ：対角要素ベクトルを返す
   diag(B)  ：ベクトルを対角要素に持つ行列を返す
   trans(A) ：転置行列を返す
   trans(B) ：ベクトルを列行列<1,Size>として返す
   matrix(B)：ベクトルを行行列<Size,1>として返す
   matmul(B1,B2)：列ベクトル*行ベクトル=行列を返す（=matrix(B1)*trans(B2)）
   
   norm1(A)　：1-ノルムを返す
   infnorm(A)：無限大ノルムを返す
   fnorm(A)　：フロベニウスノルムを返す
   max(A)    ：最大の要素を返す
   min(A)    ：最小の要素を返す
   trace(A)  ：対角成分の和を返す
   
   chol(A,L)   ：Cholesky分解する（行列Lに書き込み）
   lu(A,p,LU)  ：LU分解する（配列p,行列LUに書き込み）
   det(A)      ：行列式を返す
   solve(A,x=y)：連立一次方程式（A*x=y）を解く（行列xに書き込み）
   inv(A,iA)   ：逆行列A^-1を求める（行列iAに書き込み）
   inv(A)      ：    〃    を返す
   expm(A,eA)  ：行列の指数関数e^Aを求める（行列eAに書き込み）
   expm(A)     ：    〃           を返す
   signm(A,S)  ：行列の符号関数を求める（行列Sに書き込み）
   signm(A)    ：    〃        を返す
   sqrtm(A,sA) ：行列の平方根√A求める（行列eAに書き込み）
   sqrtm(A)    ：    〃        を返す
   qr(A,Q,R)   ：QR分解する（行列Q,Rに書き込み）
   eig(A,er,ei)：固有値を求める（固有値のベクトルを実部er,虚部eiに書き込み）
   
   print(A)：　　 printfを使って行列を表示する
   printf(fmt,A)：↑の値の表示フォーマット指定版
*/


// 行列を表すクラス
template<class E>
class dMatrixExpression {
public:
  const E &operator() () const {
    return *static_cast<const E *>(this);
  }
  E &operator() () {
    return *static_cast<E *>(this);
  }
};


// 演算コストに応じて演算結果をキャッシュする構造体
// （行列*行列のように複数要素に渡る演算があるときにキャッシュ：Costs=1）
template<class R, bool Costs_>
class dMatrixCache:
  public dMatrixExpression< dMatrixCache<R, Costs_> > {
public:
  static const bool Costs = Costs_;
  static const int Rows = R::Rows;
  static const int Cols = R::Cols;
  typedef typename R::rowvector_type rowvector_type;
private:
  const R &rm_; //no cache
public:
  dMatrixCache(const R &rm): rm_(rm) {}
  
  rowvector_type operator[] (int i) const { return rm_[i]; }
  rowvector_type operator() (int i) const { return rm_(i); }
  double operator() (int i, int j)  const { return rm_(i,j); }
};
template<class R>
class dMatrixCache<R, 1>:
  public dMatrixExpression< dMatrixCache<R, 1> > {
public:
  static const bool Costs = 0;
  static const int Rows = R::Rows;
  static const int Cols = R::Cols;
  typedef dVector<R::Cols> rowvector_type;
private:
  const dMatrix<R::Rows,R::Cols> rm_; //cache
public:
  dMatrixCache(const R &rm): rm_(rm) {}
  
  const rowvector_type &operator[] (int i) const { return rm_[i]; }
  const rowvector_type &operator() (int i) const { return rm_(i); }
  double operator() (int i, int j)  const { return rm_(i,j); }
};


// 行列の演算式を表すクラス（行列と行列の加減算）
template<class L, class Ope, class R>
class dMatrixOperation:
  public dMatrixExpression< dMatrixOperation<L, Ope, R> > {
public:
  static const bool Costs = L::Costs | R::Costs;
  static const int Rows = L::Rows;
  static const int Cols = L::Cols;
  typedef dVectorOperation<typename L::rowvector_type, Ope, typename R::rowvector_type> rowvector_type;
private:
  const L &lm_;
  const R &rm_;
public:
  dMatrixOperation(const L &lm, const R &rm): lm_(lm), rm_(rm) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(lm_[i], rm_[i]); }
  rowvector_type operator() (int i) const { return rowvector_type(lm_(i), rm_(i)); }
  double operator() (int i, int j)  const { return Ope::apply(lm_(i,j), rm_(i,j)); }
};

// 行列の演算式を表すクラス（行列と数値の乗除算）
template<class L, class Ope>
class dMatrixOperation<L, Ope, double>:
  public dMatrixExpression< dMatrixOperation<L, Ope, double> > {
public:
  static const bool Costs = L::Costs;
  static const int Rows = L::Rows;
  static const int Cols = L::Cols;
  typedef dVectorOperation<typename L::rowvector_type, Ope, double> rowvector_type;
private:
  const L &lm_;
  const double rd_;//  __attribute__((aligned(8))); //速くなる？
public:
  dMatrixOperation(const L &lm, const double rd): lm_(lm), rd_(rd) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(lm_[i], rd_); }
  rowvector_type operator() (int i) const { return rowvector_type(lm_(i), rd_); }
  double operator() (int i, int j)  const { return Ope::apply(lm_(i,j), rd_); }
};

// 行列の演算式を表すクラス（行列の負演算など）
template<class Ope, class R>
class dMatrixOperationS:
  public dMatrixExpression< dMatrixOperationS<Ope, R> > {
public:
  static const bool Costs = R::Costs;
  static const int Rows = R::Rows;
  static const int Cols = R::Cols;
  typedef dVectorOperationS<Ope, typename R::rowvector_type> rowvector_type;
private:
  const R &rm_;
public:
  dMatrixOperationS(const R &rm): rm_(rm) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(rm_[i]); }
  rowvector_type operator() (int i) const { return rowvector_type(rm_(i)); }
  double operator() (int i, int j)  const { return Ope::apply(rm_(i,j)); }
};

// 加算（行列 + 行列）
template<class L, class R>
inline
dMatrixOperation<L, OperationAdd, R> operator + (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  return dMatrixOperation<L, OperationAdd, R>( lm(), rm() );
}
// 減算（行列 - 行列）
template<class L, class R>
inline
dMatrixOperation<L, OperationSub, R> operator - (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  return dMatrixOperation<L, OperationSub, R>( lm(), rm() );
}
// 乗算（行列 * 数値）
template<class L>
inline
dMatrixOperation<L, OperationMul, double> operator * (const dMatrixExpression<L> &lm, const double rd) {
  return dMatrixOperation<L, OperationMul, double>( lm(), rd );
}
// 乗算（数値 * 行列）＝（行列 * 数値）
template<class R>
inline
dMatrixOperation<R, OperationMul, double> operator * (const double ld, const dMatrixExpression<R> &rm) {
  return dMatrixOperation<R, OperationMul, double>( rm(), ld );
}
// 除算（行列 / 数値）
template<class L>
inline
dMatrixOperation<L, OperationDiv, double> operator / (const dMatrixExpression<L> &lm, const double rd) {
  return dMatrixOperation<L, OperationDiv, double>( lm(), rd );
}
// 負（- 行列）
template<class R>
inline
dMatrixOperationS<OperationNeg, R> operator - (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationNeg, R>( rm() );
}
// 等号（行列 == 行列）// 要素はdouble型なので演算結果を比較するのはおすすめしない
template<class L, class R>
inline
bool operator == (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  for (int i = 0; i < L::Rows; i++) {
    for (int j = 0; j < L::Cols; j++) {
      if ( lm()(i,j) != rm()(i,j) ) return 0;
    }
  }
  return 1;
}
// 不等号（行列 != 行列）
template<class L, class R>
inline
bool operator != (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  for (int i = 0; i < L::Rows; i++) {
    for (int j = 0; j < L::Cols; j++) {
      if ( lm()(i,j) != rm()(i,j) ) return 1;
    }
  }
  return 0;
}

// 絶対値 fabs(行列)
template<class R>
inline
dMatrixOperationS<OperationFabs, R> fabs (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationFabs, R>( rm() );
}
// 四捨五入 round(行列)
template<class R>
inline
dMatrixOperationS<OperationRound, R> round (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationRound, R>( rm() );
}
// 天井関数 ceil(行列)
template<class R>
inline
dMatrixOperationS<OperationCeil, R> ceil (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationCeil, R>( rm() );
}
// 床関数 floor(行列)
template<class R>
inline
dMatrixOperationS<OperationFloor, R> floor (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationFloor, R>( rm() );
}
// 各要素の二乗 square(行列)
template<class R>
inline
dMatrixOperationS<OperationSquare, R> square (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationSquare, R>( rm() );
}


// 行列の演算式を表すクラス（行列の行ベクトル）
template<class L>
class dMatrixOpeRow:
  public dVectorExpression< dMatrixOpeRow<L> > {
public:
  static const bool Costs = L::Costs;
  static const int Size = L::Cols; //長さは列数
private:
  const L &lm_;
  const int ri_;
public:
  dMatrixOpeRow(const L &lm, const int ri): lm_(lm), ri_(ri) {}
  
  double operator[] (int j) const { return lm_(ri_,j); }
  double operator() (int j) const { return lm_(ri_,j); }
};

// 行列の演算式を表すクラス（行列の列ベクトル）
template<class L>
class dMatrixOpeCol:
  public dVectorExpression< dMatrixOpeCol<L> > {
public:
  static const bool Costs = L::Costs;
  static const int Size = L::Rows; //長さは行数
private:
  const L &lm_;
  const int rj_;
public:
  dMatrixOpeCol(const L &lm, const int rj): lm_(lm), rj_(rj) {}
  
  double operator[] (int i) const { return lm_(i,rj_); }
  double operator() (int i) const { return lm_(i,rj_); }
};

// i行ベクトル row(行列,行番号)
template<class L>
inline
dMatrixOpeRow<L> rowvector (const dMatrixExpression<L> &lm, int i) {
  return dMatrixOpeRow<L>( lm(), i );
}

// j列ベクトル col(行列,列番号)
template<class L>
inline
dMatrixOpeCol<L> colvector (const dMatrixExpression<L> &lm, int j) {
  return dMatrixOpeCol<L>( lm(), j );
}


// 行列*列ベクトルの演算式を表すクラス
template<class L, class R>
class dMatrixOpeMulL:
  public dVectorExpression< dMatrixOpeMulL<L, R> > {
public:
  static const bool Costs = 1;
  static const int Size = L::Rows; //長さは行列Lの行数になる
private:
  const dMatrixCache<L,L::Costs> lm_; //キャッシュで高速化
  const dVectorCache<R,R::Costs> rv_; //キャッシュで高速化
  //const L &lm_;
  //const R &rv_;
public:
  dMatrixOpeMulL(const L &lm, const R &rv): lm_(dMatrixCache<L,L::Costs>(lm)), rv_(dVectorCache<R,R::Costs>(rv)) {}
  //dMatrixOpeMulL(const L &lm, const R &rv): lm_(lm), rv_(rv) {}
  
  double operator[] (int i) const {
    /*
    double result = 0;
    for (int j = 0; j < L::Cols; j++) {
      result += lm_(i,j) * rv_[j]; //i番の要素は行列Lのi行*列ベクトル
    }
    return result;
    */
    return rowvector(lm_,i) * rv_;
  }
  double operator() (int i) const {
    return rowvector(lm_,i) * rv_;
  }
};

// 行ベクトル*行列の演算式を表すクラス
template<class L, class R>
class dMatrixOpeMulR:
  public dVectorExpression< dMatrixOpeMulR<L, R> > {
public:
  static const bool Costs = 1;
  static const int Size = R::Cols; //長さは行列Rの列数になる
private:
  const dVectorCache<L,L::Costs> lv_; //キャッシュで高速化
  const dMatrixCache<R,R::Costs> rm_; //キャッシュで高速化
  //const L &lv_;
  //const R &rm_;
public:
  dMatrixOpeMulR(const L &lv, const R &rm): lv_(dVectorCache<L,L::Costs>(lv)), rm_(dMatrixCache<R,R::Costs>(rm)) {}
  //dMatrixOpeMulR(const L &lv, const R &rm): lv_(lv), rm_(rm) {}
  
  double operator[] (int j) const {
    /*
    double result = 0;
    for (int i = 0; i < R::Rows; i++) {
      result += lv_[i] * rm_(i,j); //j番の要素は行ベクトル*行列Rのj列
    }
    return result;
    */
    return lv_ * colvector(rm_,j);
  }
  double operator() (int j) const {
    return lv_ * colvector(rm_,j);
  }
};

// 行列*行列の演算式を表すクラス
template<class L, class R>
class dMatrixOpeMul:
  public dMatrixExpression< dMatrixOpeMul<L, R> > {
public:
  static const bool Costs = 1;
  static const int Rows = L::Rows; //行数は行列Lの行数になる
  static const int Cols = R::Cols; //列数は行列Rの列数になる
  typedef dMatrixOpeMulR<typename L::rowvector_type, R> rowvector_type;
private:
  const dMatrixCache<L,L::Costs> lm_; //キャッシュで高速化
  const dMatrixCache<R,R::Costs> rm_; //キャッシュで高速化
  //const L &lm_;
  //const R &rm_;
public:
  dMatrixOpeMul(const L &lm, const R &rm): lm_(dMatrixCache<L,L::Costs>(lm)), rm_(dMatrixCache<R,R::Costs>(rm)) {}
  //dMatrixOpeMul(const L &lm, const R &rm): lm_(lm), rm_(rm) {}
  
  //行ベクトル*行列
  rowvector_type operator[] (int i) const { return rowvector_type(lm_[i], rm_); }
  rowvector_type operator() (int i) const { return rowvector_type(lm_(i), rm_); }
  
  double operator() (int i, int j) const {
    /*
    double result = 0;
    for (int k = 0; k < L::Cols; k++) {
      result += lm_(i,k) * rm_(k,j); //(i,j)番の要素は行列Lのi行*行列Rのj列
    }
    return result;
    */
    return rowvector(lm_,i) * colvector(rm_,j);
  }
};

// 乗算（行列 * 列ベクトル）（勝手に列ベクトルと解釈）
template<class L, class R>
inline
dMatrixOpeMulL<L, R> operator * (const dMatrixExpression<L> &lm, const dVectorExpression<R> &rv) {
  return dMatrixOpeMulL<L, R>( lm(), rv() );
}
// 乗算（行ベクトル * 行列）（勝手に行ベクトルと解釈）
template<class L, class R>
inline
dMatrixOpeMulR<L, R> operator * (const dVectorExpression<L> &lv, const dMatrixExpression<R> &rm) {
  return dMatrixOpeMulR<L, R>( lv(), rm() );
}
// 乗算（行列 * 行列）
template<class L, class R>
inline
dMatrixOpeMul<L, R> operator * (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  return dMatrixOpeMul<L, R>( lm(), rm() );
}


// 行列の演算式を表すクラス（行列の対角要素ベクトル）
template<class R>
class dMatrixOpeDiag:
  public dVectorExpression< dMatrixOpeDiag<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Size = R::Rows; //長さは行数
private:
  const R &rm_;
public:
  dMatrixOpeDiag(const R &rm): rm_(rm) {}
  
  double operator[] (int i) const { return rm_(i,i); }
  double operator() (int i) const { return rm_(i,i); }
};

// ベクトルの演算式を表すクラス（ベクトルを対角要素に持つ行列）
template<class R>
class dVectorOpeDiag:
  public dMatrixExpression< dVectorOpeDiag<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Rows = R::Size;
  static const int Cols = R::Size;
  typedef dVectorOperation<dStdBasisVector<R::Size>, OperationMul, double> rowvector_type;
private:
  const R &rv_;
public:
  dVectorOpeDiag(const R &rv): rv_(rv) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(dStdBasisVector<R::Size>(i), rv_[i]); }
  rowvector_type operator() (int i) const { return rowvector_type(dStdBasisVector<R::Size>(i), rv_(i)); }
  double operator() (int i, int j) const { return (i==j)? rv_[i] : 0; }
};

// 対角要素ベクトル diag(行列)
template<class R>
inline
dMatrixOpeDiag<R> diag (const dMatrixExpression<R> &rm) {
  return dMatrixOpeDiag<R>( rm() );
}

// ベクトルを対角要素に持つ行列 diag(ベクトル)
template<class R>
inline
dVectorOpeDiag<R> diag (const dVectorExpression<R> &rv) {
  return dVectorOpeDiag<R>( rv() );
}


// 行列の演算式を表すクラス（行列の転置）
template<class R>
class dMatrixOpeTrans:
  public dMatrixExpression< dMatrixOpeTrans<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Rows = R::Cols; //転置なので逆になる
  static const int Cols = R::Rows; //
  typedef dMatrixOpeCol<R> rowvector_type;
private:
  const R &rm_;
public:
  dMatrixOpeTrans(const R &rm): rm_(rm) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(rm_, i); }
  rowvector_type operator() (int i) const { return rowvector_type(rm_, i); }
  double operator() (int i, int j) const { return rm_(j,i); }
};

// 転置 trans(行列)
template<class R>
inline
dMatrixOpeTrans<R> trans (const dMatrixExpression<R> &rm) {
  return dMatrixOpeTrans<R>( rm() );
}

// ベクトルの演算式を表すクラス（ベクトルを行列<1,Size>とみなす）
template<class R>
class dVectorOpeTrans:
  public dMatrixExpression< dVectorOpeTrans<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Rows = 1;
  static const int Cols = R::Size;
  typedef R rowvector_type;
private:
  const R &rv_;
public:
  dVectorOpeTrans(const R &rv): rv_(rv) {}
  
  rowvector_type operator[] (int i) const { return rv_; }//i=0以外も同じになるので注意
  rowvector_type operator() (int i) const { return rv_; }
  double operator() (int i, int j) const { return rv_(j); }
};

// 転置（ベクトルを行列<1,Size>とみなす）trans(ベクトル)
template<class R>
inline
dVectorOpeTrans<R> trans (const dVectorExpression<R> &rv) {
  return dVectorOpeTrans<R>( rv() );
}

// ベクトルの演算式を表すクラス（ベクトルを行列<Size,1>とみなす）
template<class R>
class dVectorOpeMatrix:
  public dMatrixExpression< dVectorOpeMatrix<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Rows = R::Size;
  static const int Cols = 1;
  typedef dSubVector<R, 1> rowvector_type;
private:
  const R &rv_;
public:
  dVectorOpeMatrix(const R &rv): rv_(rv) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(rv_, i); }
  rowvector_type operator() (int i) const { return rowvector_type(rv_, i); }
  double operator() (int i, int j) const { return rv_[i]; }
};

// ベクトルを行列<Size,1>と見なす matrix(ベクトル)
template<class R>
inline
dVectorOpeMatrix<R> matrix (const dVectorExpression<R> &rv) {
  return dVectorOpeMatrix<R>( rv() );
}


// ベクトルの演算式を表すクラス（列ベクトル*行ベクトル→行列）
template<class L, class R>
class dVectorOpeMul:
  public dMatrixExpression< dVectorOpeMul<L, R> > {
public:
  static const bool Costs = L::Costs | R::Costs;
  static const int Rows = L::Size;
  static const int Cols = R::Size;
  typedef dVectorOperation<R, OperationMul, double> rowvector_type;
private:
  const L &lv_;
  const R &rv_;
public:
  dVectorOpeMul(const L &lv, const R &rv): lv_(lv), rv_(rv) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(rv_, lv_[i]); }
  rowvector_type operator() (int i) const { return rowvector_type(rv_, lv_(i)); }
  double operator() (int i, int j) const { return lv_[i] * rv_[j]; }
};

// 乗算（列ベクトル*行ベクトル→行列）（勝手に列ベクトル,行ベクトルと解釈）
template<class L, class R>
inline
dVectorOpeMul<L, R> matmul (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  return dVectorOpeMul<L, R>( lv(), rv() );
}


/**************************************************************************/


// 1-ノルム norm1(行列)
template<class R>
inline
double norm1 (const dMatrixExpression<R> &rm) {
  // 最大列和： max_j( sum_i( |rm(i,j)| ))
  int i, j;
  double norm=0, sum;
  
  for(j=0; j<R::Cols; j++) {
    sum = 0;
    for(i=0; i<R::Rows; i++) {
      sum += fabs(rm()(i,j));
    }
    if (norm < sum) norm = sum;
  }
  return norm;
}

// 無限大ノルム infnorm(行列)
template<class R>
inline
double infnorm (const dMatrixExpression<R> &rm) {
  // 最大行和： max_i( sum_j( |rm(i,j)| ))
  int i, j;
  double norm=0, sum;
  
  for(i=0; i<R::Rows; i++) {
    sum = 0;
    for(j=0; j<R::Cols; j++) {
      sum += fabs(rm()(i,j));
    }
    if (norm < sum) norm = sum;
  }
  return norm;
}

// フロベニウスノルム fnorm(行列)
template<class R>
inline
double fnorm (const dMatrixExpression<R> &rm) {
  // 二乗和平方根： √sum_i( sum_j( rm(i,j)^2 ))
  int i, j;
  double sum=0;
  
  for(i=0; i<R::Rows; i++) {
    for(j=0; j<R::Cols; j++) {
      sum += OperationSquare::apply(rm()(i,j));
    }
  }
  return sqrt(sum);
}

// 最大値 max(行列)
template<class R>
inline
double max (const dMatrixExpression<R> &rm) {
  double result = -DBL_MAX;
  double tmp;
  for (int i=0; i<R::Rows; i++) {
    for (int j=0; j<R::Cols; j++) {
      tmp = rm()(i,j);
      if (result < tmp) result = tmp;
    }
  }
  return result;
}

// 最小値 min(行列)
template<class R>
inline
double min (const dMatrixExpression<R> &rm) {
  double result = DBL_MAX;
  double tmp;
  for (int i=0; i<R::Rows; i++) {
    for (int j=0; j<R::Cols; j++) {
      tmp = rm()(i,j);
      if (result > tmp) result = tmp;
    }
  }
  return result;
}

// トレース（対角成分の和） trance(行列)
template<class R>
inline
double trace (const dMatrixExpression<R> &rm) {
  double sum=0;
  for (int i=0; i<R::Rows; i++) {
    sum += rm()(i,i);
  }
  return sum;
}


/**************************************************************************/

// Cholesky分解 chol(行列,下三角行列) A = L*trans(L)
template <class E>
int chol(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &L) {
  int i, j, k;
  double sum;
  
  for (j=0; j<E::Rows; j++) {
    // L(j,j) の計算
    sum = A()(j,j);
    for (k=0; k<j; k++) {
      sum -= L(j,k) * L(j,k);
    }
    if (sum == 0) return 0;
    
    L(j,j) = sqrt(sum);
    
    // L(i,j) (i>j) の計算
    for (i=j+1; i<E::Rows; i++) {
      sum = A()(i,j);
      for (k=0; k<j; k++) {
        sum -= L(i,k) * L(j,k);
      }
      L(i,j) = sum / L(j,j);
      L(j,i) = 0; //上三角成分は0
   }
  }
  return 1;
}

// LU分解 lu(行列,行交換情報の配列,LU行列) A = L*U
template <class E>
inline
double lu(const dMatrixExpression<E> &A, int pivot[], dMatrix<E::Rows,E::Rows> &LU) {
  int i, j, k, ii, kk;
  double tmp;
  double det;
  double w;
  double weight[E::Rows];
  
  LU = A;
  for (i=0; i<E::Rows; i++) {
    w = 0;
    for (j=0; j<E::Rows; j++) { //絶対値最大の要素を求める
      tmp = fabs(LU(i,j));
      if (w<tmp) w = tmp;
    }
    if (w==0) return 0; //singular
    
    weight[i] = 1/w; //行の重みは最大絶対値の逆数
    pivot[i] = i;    //行交換情報の初期値
  }
  
  det = 1;
  for (i=0; i<E::Rows; i++) {
    w = -1;
    for (k=i; k<E::Rows; k++) { //より下の各行について
      kk = pivot[k];            //重み×絶対値が最大の行を見つける
      tmp = fabs(LU(kk,i)) * weight[kk];
      if (w<tmp) { w=tmp; j=k; } //jは重み×絶対値が最大の行
    }
    ii = pivot[j];
    if (j != i) {
      pivot[j] = pivot[i]; pivot[i] = ii; //行番号を交換
      det = -det; //行を交換すれば行列式の符号が変わる
    }
    w = LU(ii,i); det *= w; //行列式はUの対角成分の積
    if (w==0) return 0;
    
    for (k=i+1; k<E::Rows; k++) {  //Gauss消去法
      kk = pivot[k];
      tmp = (LU(kk,i) /= w);      //下三角行列Lの要素を計算
      #pragma ivdep
      for (j=i+1; j<E::Rows; j++) {//上三角行列Uの要素を計算
	LU(kk,j) -= tmp * LU(ii,j);
      }
    }
  }
  return det; //戻り値は行列式
}

// LU分解で行列式を求める det(行列)
template <class E>
inline
double det(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Cols> LU;
  int pivot[E::Rows]; 
  
  return lu(A, pivot, LU);
}

// LU分解を用いた連立一次方程式計算 solve(LU行列,行交換情報の配列,解)
template <int Size_, int Num_>
inline
void solve(const dMatrix<Size_, Size_> &LU, const int pivot[], dMatrix<Size_,Num_> &solA) {
  //・A*x = y という連立一次方程式は，L*U*x = y を解くことと等価
  //  L*c = y → U*x = c と計算することで求められる（L,Uは下三角,上三角なのでc,xの計算簡単）
  //・係数Aが同じ複数(Num_個)の連立方程式を同時に解ける
  //・x = A^-1 * y を計算するのはこちらの方が逆行列を求めるより速い＆高精度
  //  逆行列を求めるのはyを単位行列(ただし行交換必要)にした場合と等価
  
  int i, j, k, ii;
  double c[Size_];
  double sum;
  
  for (k=0; k<Num_; k++) {
    //前進代入
    for (i=0; i<Size_; i++) {
      ii = pivot[i];
      sum = solA(ii,k);
      for (j=i-1; j>=0; j--) {
	sum -= LU(ii,j) * c[j];
      }
      c[i] = sum;
    }
    //後退代入
    for (i=Size_-1; i>=0; i--) {
      ii = pivot[i];
      sum = c[i];
      for (j=i+1; j<Size_; j++) {
	sum -= LU(ii,j) * solA(j,k);
      }
      solA(i,k) = sum / LU(ii,i);
    }
  }
}

// 連立一次方程式（A*x=y）を解く solve(A,x = y) → 解x=A^-1*y
template <class E, int Num_>
inline
int solve(const dMatrixExpression<E> &A, dMatrix<E::Rows,Num_> &solA) {
  dMatrix<E::Rows,E::Cols> LU;
  int pivot[E::Rows]; 
  
  if (lu(A, pivot, LU) == 0) return 0;
  
  solve(LU, pivot, solA);
  
  return 1;
}

// 連立一次方程式（A*x=y）を解く（２×２用）
template <class E, int Num_>
inline
int solve(const dMatrixExpression<E> &A, dMatrix<2,Num_> &solA) {
  dMatrix<2,2> invA;
  dMatrix<2,Num_> tmp; //同じ変数に代入できないので代入経由用
  
  //逆行列を求めた方が速い
  if (inv(A, invA) == 0) return 0;
  
  tmp = invA*solA;
  solA = tmp;
  return 1;
}

// 逆行列 inv(行列,逆行列)
template <class E>
inline
double inv(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &invA) {
  int i, j, k, ii;
  double sum;
  
  dMatrix<E::Rows,E::Rows> LU;
  int pivot[E::Rows]; 
  
  double det = lu(A, pivot, LU);
  if (det == 0) return 0;
  
  //単位行列の行交換をしつつsolve()
  for (k=0; k<E::Rows; k++) {
    //前進代入
    for (i=0; i<E::Rows; i++) {
      ii = pivot[i];
      sum = (ii == k);//行交換された単位行列
      for (j=i-1; j>=0; j--) {
	sum -= LU(ii,j) * invA(j,k);
      }
      invA(i,k) = sum;
    }
    //後退代入
    for (i=E::Rows-1; i>=0; i--) {
      ii = pivot[i];
      sum = invA(i,k);
      for (j=i+1; j<E::Rows; j++) {
	sum -= LU(ii,j) * invA(j,k);
      }
      invA(i,k) = sum / LU(ii,i);
    }
  }
  return det; //戻り値は行列式
}

// 逆行列（２×２用）
template <class E>
inline
double inv(const dMatrixExpression<E> &A, dMatrix<2,2> &invA) {
  invA(0,0) = A()(1,1);
  invA(0,1) =-A()(1,0);
  invA(1,0) =-A()(0,1);
  invA(1,1) = A()(0,0);
  
  double det = invA(0,0)*invA(1,1)-invA(0,1)*invA(1,0);
  if (det == 0) return 0;
  invA /= det;  return det; //戻り値は行列式
}

// 逆行列：オブジェクトのコピーが生じるので変数に代入する場合は inv(A,invA) を使うべし
template <class E>
inline
dMatrix<E::Rows,E::Rows> inv(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> invA;
  
  inv(A, invA);
  
  return invA;
}

// 行列の指数関数 expm(行列,exp行列)
template <class E>
inline
int expm(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &U){
  // exp(A)の計算
  // テイラー展開で計算すると誤差大→パデ近似を用いた手法
  // exp(A) = { Pade(A/2^s) }^2s
  // 参考:
  // Computing Matrix Exponentials using C++ with boost library
  // https://www.dbtsai.com/blog/2008-11-25-matrix-exponential/
  
  const int p = 6; // recommended and gererally satisfactory
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  int i, j, k;
  dMatrix<E::Rows,E::Rows> P, Q, A2;
  dMatrix<E::Rows,E::Rows> tmp;//同じ変数に代入できないので代入経由用
  
  // Calcuate Pade coefficient
  double c[p+1];
  c[0] = 1;
  for(i=0; i<p; i++) { 
    c[i+1] = c[i] * ((p-i)/((i+1)*(2.0*p-i)));
  }
  
  // Calcuate the infinty norm of A,
  // which is defined as the largest row sum of a matrix
  double norm = infnorm(A);
  
  // If norm = 0, and all A elements are not NaN or infinity but zero, 
  // then U should be identity.
  if (norm == 0.0) {
    if (A == 0.0*I) U = I;
    else            U = 0.0;
    return 0;
  }
  
  // Scaling, seek s such that ||A/2^s|| < 1/2, and set scale = 1/2^s
  int s = 0;
  if(norm > 0.5) {
    s = (int)((log2(norm) + 2.0));
    
    U = A/(1<<s);
  } else {
    U = A;
  }
  
  // Horner evaluation of the irreducible fraction.
  // Initialise P (numerator) and Q (denominator) 
  A2 = U * U;
  Q = c[p];//*I;
  P = c[p-1];//*I;
  int odd = 1;
  for(k=p-2; k>=0; k--) {
    if (odd == 1) {
      tmp = Q * A2 + c[k] * I;
      Q = tmp;
    } else {
      tmp = P * A2 + c[k] * I;
      P = tmp;
    }
    odd = 1-odd;
  }
  if (odd == 1) {
    tmp = Q * U;
    Q = tmp;
  } else {
    tmp = P * U;
    P = tmp;
  }
  
  solve(Q-P,tmp = 2.0*P); //tmp=2.0*(Q-P)^-1*P
  
  if (odd == 1) {
    U = -(I + tmp); 
  } else {
    U =  (I + tmp);
  }
  
  // Squaring 
  for(i=0; i<s; i++) {
    tmp = U * U;
    U = tmp;
  }
  return 1;
}

// 行列の指数関数：オブジェクトのコピーが生じるので変数に代入する場合は expm(A,eA) を使うべし
template <class E>
inline
dMatrix<E::Rows,E::Rows> expm(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> eA;
  
  expm(A, eA);
  
  return eA;
}


// 行列の符号関数 signm(行列,sign行列)
template <class E>
inline
int signm(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &S){
  // Scaled Newton iteration
  // 参考:
  // Ralph Byers, Hongguo Xu, "A New Scaling for Newton's Iteration for
  // the Polar Decomposition and its Backward Stability"
  
  dMatrix<E::Rows,E::Rows> Sk, Zk; S=A;
  double mu;
  
  do {
    Sk = S;
    
    //if (inv(Sk, Zk) == 0) return 0; //高速化のためエラー処理を省く
    inv(Sk, Zk);
    
    mu = sqrt(fnorm(Zk)/fnorm(Sk)); 
    //mu = pow((norm1(Zk)*infnorm(Zk))/(norm1(Sk)*infnorm(Sk)), 0.25); 
    //mu = 1.0/pow(fabs(det), 1.0/E::Rows); 
    //mu = 1; //scalingなし
    
    S = 0.5 * (mu*Sk + 1.0/mu*Zk);
  } while (fnorm(S-Sk) > DBL_EPSILON*fnorm(Sk));
  
  return 1;
}

// 行列の符号関数：オブジェクトのコピーが生じるので変数に代入する場合は signm(A,S) を使うべし
template <class E>
inline
dMatrix<E::Rows,E::Rows> signm(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> S;
  
  signm(A, S);
  
  return S;
}


// 行列の平方根 sqrtm(行列)
template <class E>
inline
int sqrtm(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &sA) {
  // 参考:
  // Eugene D. Denman, "The Matrix Sign Function and Computations in Systems"
    
  int i,j;
  
  // H = [ 0 A ]
  //     [ I 0 ]
  dMatrix<E::Rows*2,E::Rows*2> H; H.clear();
  
  submatrix(H, E::Rows,0, E::Rows,E::Rows) = 1.0;//*I;
  submatrix(H, 0,E::Rows, E::Rows,E::Rows) = A;
  
  // F = signm(H) + [I  0] = [  I  S12 ]
  //                [0 -I]   [ S21  -I ]
  dMatrix<E::Rows*2,E::Rows*2> F;
  signm(H, F);
  /*
  for (i=0; i<E::Rows; i++) { //√A の計算では不要
    F(i,i) += 1;
    F(i+E::Rows,i+E::Rows) -= 1;
  }
  */
  sA = submatrix(F, 0,E::Rows, E::Rows,E::Rows); // √A = S12
  
  return 1;
}

// 行列の平方根：オブジェクトのコピーが生じるので変数に代入する場合は sqrtm(A,S) を使うべし
template <class E>
inline
dMatrix<E::Rows,E::Rows> sqrtm(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> sA;
  
  sqrtm(A, sA);
  
  return sA;
}


// QR分解 qr(行列,Q行列,R行列) A = Q*R
template <class E>
inline
int qr(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &Q, dMatrix<E::Rows,E::Cols> &R) {
  // from Numerical Recipies
  
  int i,j,k;
  double tmp,sum,scale,sigma,c;
  
  Q = 1.0;//*I;
  R = A;
  
  for (k=0; k<E::Cols-1; k++) {
    for (scale=0,i=k; i<E::Rows; i++) {
      tmp = fabs(R(i,k));
      if (scale < tmp) scale = tmp; //最大値を求める
    }
    if (scale == 0.0) { //Singular case
      return 0; 
    } else {
      for (sum=0,i=k; i<E::Rows; i++) {
	sum += OperationSquare::apply( R(i,k)/=scale );
      }
      sigma = (R(k,k) >= 0.0)? sqrt(sum) : -sqrt(sum);
      R(k,k) += sigma;
      c = sigma*R(k,k);
      
      for (j=k+1; j<E::Cols; j++) {
	for (sum=0,i=k; i<E::Rows; i++) sum += R(i,k)*R(i,j);
	sum /= c;
	for (i=k; i<E::Rows; i++) R(i,j) -= sum*R(i,k);
      }
      
      for (j=0; j<E::Rows; j++) {
	for (sum=0,i=k; i<E::Rows; i++) sum += R(i,k)*Q(j,i);
	sum /= c;
	for (i=k; i<E::Rows; i++) Q(j,i) -= sum*R(i,k);
      }
      
      R(k,k) = -scale*sigma;
    }
  }
  for (j=0; j<E::Cols; j++) {
    for (i=j+1; i<E::Rows; i++) {
      R(i,j) = 0; //下三角は0
    }
  }
  return 1;
}


// Balancing (subfunction for eig)
template<int N>
inline
void balance(dMatrix<N,N> &A) {
  // from Numerical Recipies
  
  const double RADIX = 2.0;
  
  int last,j,i;
  double s,r,g,f,c,sqrdx;
  
  sqrdx=RADIX*RADIX;
  last=0;
  do {
    last=1;
    for (i=0;i<N;i++) { //Calculate row and column norms.
      r=c=0.0;
      for (j=0;j<N;j++) {
	if (j != i) {
	  c += fabs(A(j,i));
	  r += fabs(A(i,j));
	}
      }
      if (c && r) { //If both are nonzero,
	g=r/RADIX;
	f=1.0;
	s=c+r;
	while (c<g) {
	  // find the integer power of the machine radix that
	  //  comes closest to balancing the matrix.
	  f *= RADIX;
	  c *= sqrdx;
	}
	g=r*RADIX;
	while (c>g) {
	  f /= RADIX;
	  c /= sqrdx;
	}
	if ((c+r)/f < 0.95*s) {
	  last=0;
	  g=1.0/f;
	   //Apply similarity transformation
	  for (j=0;j<N;j++) A(i,j) *= g;
	  for (j=0;j<N;j++) A(j,i) *= f;
	}
      }
    }
  } while (last == 0);
}

// Hessenberg行列に変換（subfunction for eig）
template<int N>
inline
void hess(dMatrix<N,N> &A) {
  // from Numerical Recipies
  
  int m,j,i;
  double  y,x,tmp;
  
  for (m=1;m<N-1;m++) { //m is called r + 1 in the text.
    x = 0.0;
    i = m;
    for (j=m; j<N; j++) { //Find the pivot.
      if (fabs(A(j,m-1)) > fabs(x)) {
	x = A(j,m-1);
	i = j;
      }
    }
    if (i != m) {// Interchange rows and columns.
      for (j=m-1; j<N; j++) {
	tmp = A(i,j);
	A(i,j) = A(m,j);
	A(m,j) = tmp;
      }
      for (j=0; j<N; j++) {
	tmp = A(j,i);
	A(j,i) = A(j,m);
	A(j,m) = tmp;
      }
    }
    if (x != 0.0) { //Carry out the elimination.
      for (i=m+1; i<N; i++) {
	if ((y=A(i,m-1)) != 0.0) {
	  y /= x;
	  A(i,m-1) = y;
	  for (j=m; j<N; j++) A(i,j) -= y*A(m,j);
	  for (j=0; j<N; j++) A(j,m) += y*A(j,i);
	}
      }
    }
  }
  /*
  for (j=0; j<N; j++) {
    for (i=j+2; i<N; i++) {
      A(i,j) = 0; // i>j+1は0（eigの計算では不要）
    }
  }
  */
}

// 固有値 eig (行列,実部ベクトルei,虚部ベクトルei)
template<class E>
inline
int eig(const dMatrixExpression<E> &A, dVector<E::Rows> &er, dVector<E::Rows> &ei) {
  // from Numerical Recipies
  
  dMatrix<E::Rows,E::Rows> a = A;
  
  balance(a);
  hess(a);
  
  int nn,m,l,k,j,its,i,mmin;
  double z,y,x,w,v,u,t,s,r,q,p,anorm;
  
  // Compute matrix norm for possible use in
  //  locatings single small subdiagonal element.
  anorm=0.0; 
  for (i=0; i<E::Rows; i++) {
    for (j=(i>1)? i-1:0; j<E::Rows; j++) {
      anorm += fabs(a(i,j));
    }
  }
  nn = E::Rows-1;
  t = 0.0; //Gets changed only by an exceptional shift.
  while (nn >= 0) { //Begin search for next eigenvalue.
    its = 0;
    do {
      //Begin iteration: look for single small subdiagonal element.
      for (l=nn; l>0; l--) {
	s = fabs(a(l-1,l-1)) + fabs(a(l,l)); 
	if (s == 0.0) s = anorm;
	if ((fabs(a(l,l-1)) + s) == s) break;
      }
      x = a(nn,nn);
      if (l == nn) {       //One root found.
	ei[nn] = 0.0;
	er[nn] = x+t;
	nn--;
      } else {
	y = a(nn-1,nn-1);
	w = a(nn,nn-1) * a(nn-1,nn);
	if (l == nn-1) {   //Two roots found...
	  p = 0.5*(y-x);
	  q = p*p+w;
	  z = sqrt(fabs(q));
	  x += t;
	  if (q >= 0.0) {  //...a real pair.
	    z = p + ((p>0.0)? fabs(z):-fabs(z));
	    ei[nn-1] = ei[nn] = 0.0;
	    er[nn-1] = er[nn] = x+z;
	    if (z != 0.0) er[nn] = x-w/z;
	  } else {         //...a complex pair.
	    ei[nn-1]= -(ei[nn]=z);
	    er[nn-1] = er[nn] = x+p;
	  }
	  nn -= 2;
	} else {           //No roots found. Continue iteration.
	  if (its == 30) return 0; //Too many iterations in hqr
	  if (its == 10 || its == 20) { //Form exceptional shift.
	    t += x;
	    for (i=0; i<=nn; i++) a(i,i) -= x;
	    s = fabs(a(nn,nn-1)) + fabs(a(nn-1,nn-2));
	    y = x = 0.75*s;
	    w = -0.4375*s*s;
	  }
	  its++;
	  
	  //Form shift and then look for
	  // 2 consecutive small subdiagonal elements.
	  for (m=nn-2; m>=l; m--) {
	    z = a(m,m);
	    r = x-z;
	    s = y-z;
	    p = (r*s-w)/a(m+1,m) + a(m,m+1);
	    q = a(m+1,m+1)-z-r-s;
	    r = a(m+2,m+1);
	    s = fabs(p)+fabs(q)+fabs(r); //Scale to prevent overflow or underflow.
	    p /= s; 
	    q /= s;
	    r /= s;
	    if (m == l) break;
	    u = fabs(a(m,m-1))*(fabs(q)+fabs(r));
	    v = fabs(p)*(fabs(a(m-1,m-1))+fabs(z)+fabs(a(m+1,m+1)));
	    if (u+v == v) break;
	  }
	  for (i=m; i<nn-1; i++) {
	    a(i+2,i) = 0.0;
	    if (i != m) a(i+2,i-1) = 0.0;
	  }
	  for (k=m; k<nn; k++) {
	    //Double QR step on rows l to nn and columns m to nn.
	    if (k != m) {
	      p = a(k,k-1); //Begin setup of Householder vector.
	      q = a(k+1,k-1); 
	      r = 0.0;
	      if (k != nn-1) r = a(k+2,k-1);
	      if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
		p /= x; //Scale to prevent overflow or underflow.
		q /= x; 
		r /= x;
	      }
	    }
	    s = (p>0)? fabs(sqrt(p*p+q*q+r*r)):-fabs(sqrt(p*p+q*q+r*r));
	    if (s != 0.0) {
	      if (k == m) {
		if (l != m)
		  a(k,k-1) = -a(k,k-1);
	      } else
		a(k,k-1) = -s*x;
	      p += s;
	      x = p/s;
	      y = q/s;
	      z = r/s;
	      q /= p;
	      r /= p;
	      for (j=k; j<=nn; j++) { //Row modification.
		p = a(k,j) + q*a(k+1,j);
		if (k != nn-1) {
		  p += r*a(k+2,j);
		  a(k+2,j) -= p*z;
		}
		a(k+1,j) -= p*y;
		a(k,j) -= p*x;
	      }
	      mmin = (nn<k+3)? nn : k+3;
	      for (i=l; i<=mmin; i++) { //Column modification.
		p = x*a(i,k) + y*a(i,k+1);
		if (k != nn) {
		  p += z*a(i,k+2);
		  a(i,k+2) -= p*r;
		}
		a(i,k+1) -= p*q;
		a(i,k) -= p;
	      }
	    }
	  }
	}
      }
    } while (l <= nn);
  }
  return 1;
}


/**************************************************************************/

// 表示用
template<class E>
void printf(const char *fmt, const dMatrixExpression<E> &dm) {
  for (int i = 0; i < E::Rows; i++) {
    printf((i==0)? "{{":" {");
    for (int j = 0; j < E::Cols-1; j++) {
      printf(fmt, dm()(i,j));
      printf(", ");
    }
    printf(fmt, dm()(i,E::Cols-1));
    printf((i==E::Rows-1)? "}}\n":"},\n");
  }
}

// 表示用
template<class E>
inline
void print(const dMatrixExpression<E> &dm) {
  printf("% e", dm);
}
