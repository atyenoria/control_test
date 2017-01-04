/* Time-stamp: <2012-12-26 10:56:05 kobayasi>
   
   行列クラス
   ver.2012-12-23  by K.Kobayashi
   
   高速な四則演算ができる２次元配列（行列）
   
   ＜使い方＞
   dMatrix<3,3> m;  // 定義：3行x3列のdouble型行列
   
   m = 1.0;         // すべての対角成分に代入．他は0．
                    // 　単位行列が隠れているだけ（ m=1.0*I ）
   m(1,2) = 1.0;    // (1,2)の要素に代入．m[1][2] = 1.0 でも結果は同だが，
                    // 　行ベクトルから値を取ることになるため，m(1,2)の方が速い．
   d = (m * m)(1,2);// 演算して直接(1,2)を取ることもOK
                    // 　(1,2)に必要な分の計算しか行わないので高速．
		    
   ＜注意＞
   ・演算の際に行列やベクトルのサイズが正しいかどうかは判断しない（高速化のため）
   ・演算内容によって列ベクトルだったり行ベクトルだったりする
   ・行列*行列，(行)ベクトル*行列，行列*(列)ベクトルの結果を演算に使った変数に代入できない
   
   使える演算は MatrixExpression.h を参照．
*/

#ifndef __MATRIXCLASS__
#define __MATRIXCLASS__

#include "VectorClass.h"

// 行列を表すクラス(前方宣言)
template<class E>
class dMatrixExpression;

// 単位行列クラス(前方宣言)
template <int Rows_, int Cols_>
class dIdentityMatrix;


// 行列クラス
template <int Rows_, int Cols_>
class dMatrix:
  public dMatrixExpression< dMatrix<Rows_, Cols_> > {
public:
  static const bool Costs = 0;
  static const int Rows = Rows_;
  static const int Cols = Cols_;
  typedef dVector<Cols_> rowvector_type;
private:
  rowvector_type data_[Rows_];
public:
  // コンストラクタ（引数なし）
  dMatrix() {}
  
  // コピーコンストラクタ（行列をコピー）
  dMatrix(const dMatrix &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] = dm[i];
    }
  }
  dMatrix(const dIdentityMatrix<Rows_,Cols_> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
      data_[i][i] = 1;
    }
  }
  template<class E>
  dMatrix(const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
	data_[i][j] = dm()(i,j);
      }
    }
  }
  // コピーコンストラクタ（2次元配列をコピー）
  dMatrix(const double dpp[Rows_][Cols_]) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] = dpp[i];
    }
  }
  // コピーコンストラクタ（数値(*単位行列)で埋める）
  dMatrix(const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
      data_[i][i] = d;
    }
  }
  
  /** function *************************************************************/
  
  // 行列の行数
  int rows() const {
    return Rows_;
  }
  // 行列の列数
  int cols() const {
    return Cols_;
  }
  // 行列をクリア
  void clear() {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
    }
  }
  // 内部のベクトル配列を返す
  rowvector_type *data() {
    return data_;
  }
  
  /** operator *************************************************************/
  
  // i行のベクトル
  const rowvector_type &operator[] (int i) const {
    return data_[i];
  }
  rowvector_type &operator[] (int i) {
    return data_[i];
  }
  const rowvector_type &operator() (int i) const {
    return data_[i];
  }
  rowvector_type &operator() (int i) {
    return data_[i];
  }
  // (i,j)番目の要素
  double operator() (int i, int j) const {
    return data_[i][j];
  }
  double &operator() (int i, int j) {
    return data_[i][j];
  }
  // 代入（ = 行列）
  dMatrix &operator = (const dMatrix &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] = dm[i];
    }
    return *this;
  }
  dMatrix &operator = (const dIdentityMatrix<Rows_,Cols_> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
      data_[i][i] = 1;
    }
    return *this;
  }
  template<class E>
  dMatrix &operator = (const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
	data_[i][j] = dm()(i,j);
      }
    }
    return *this;
  }
  // 代入（ = ２次元配列）
  dMatrix &operator = (const double dpp[Rows_][Cols_]) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] = dpp[i];
    }
    return *this;
  }
  // 代入（ = 数値(*単位行列)）
  dMatrix &operator = (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
      data_[i][i] = d;
    }
    return *this;
  }
  // 代入（ += 行列）
  template<class E>
  dMatrix &operator += (const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
	data_[i][j] += dm()(i,j);
      }
    }
    return *this;
  }
  // 代入（ -= 行列）
  template<class E>
  dMatrix &operator -= (const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
	data_[i][j] -= dm()(i,j);
      }
    }
    return *this;
  }
  // 代入（ *= 行列）
  template<class E>
  dMatrix &operator *= (const dMatrixExpression<E> &dm) {
    dMatrix<Rows_, Cols_> result = (*this) * dm;
    return *this = result; //非効率だが上書き計算できないのでこうする
  }
  // 代入（ *= 数値）
  dMatrix &operator *= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] *= d;
    }
    return *this;
  }
  // 代入（ /= 数値）
  dMatrix &operator /= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] /= d;
    }
    return *this;
  }
};


// 単位行列のクラス（内部配列を持たないのでメモリ節約）
template <int Rows_, int Cols_>
class dIdentityMatrix:
  public dMatrixExpression< dIdentityMatrix<Rows_, Cols_> > {
public:
  static const bool Costs = 0;
  static const int Rows = Rows_;
  static const int Cols = Cols_;
  typedef dStdBasisVector<Cols_> rowvector_type;
public:
  // コンストラクタ（引数なし）
  dIdentityMatrix() {}
  
  /** function *************************************************************/
  
  // 行列の行数
  int rows() const {
    return Rows_;
  }
  // 行列の列数
  int cols() const {
    return Cols_;
  }
  
  /** operator *************************************************************/
  
  // i行のベクトル
  rowvector_type operator[] (int i) const {
    return rowvector_type(i);
  }
  rowvector_type operator() (int i) const {
    return rowvector_type(i);
  }
  // (i,j)番目の要素
  double operator() (int i, int j) const {
    return (i == j);
  }
};


// 部分行列を表すクラス
template <class L, int Rows_, int Cols_>
class dSubMatrix:
  public dMatrixExpression< dSubMatrix<L, Rows_, Cols_> > {
public:
  static const bool Costs = 0;
  static const int Rows = Rows_;
  static const int Cols = Cols_;
  typedef dSubVector<typename L::rowvector_type, Cols_> rowvector_type;
private:
  L &lm_;
  const int rowoffset_;
  const int coloffset_;
public:
  dSubMatrix(const L &lm, const int rowoffset, const int coloffset): lm_(*const_cast<L*>(&lm)), rowoffset_(rowoffset), coloffset_(coloffset) {}  //行列のconstをはずす
  
  /** function *************************************************************/
  
  // 行列の行数
  int rows() const {
    return Rows_;
  }
  // 行列の列数
  int cols() const {
    return Cols_;
  }
  // 行列をクリア
  void clear() {
    #pragma ivdep
    for (int i = rowoffset_; i < rowoffset_+Rows_; i++) {
      memset(lm_[i].data()+coloffset_, 0, sizeof(double)*Cols_);
    }
  }
  
  /** operator *************************************************************/
  
  // i行のベクトル
  rowvector_type operator[] (int i) const {
    return subvector(lm_[i+rowoffset_], coloffset_, Cols_);
  }
  rowvector_type operator() (int i) const {
    return subvector(lm_[i+rowoffset_], coloffset_, Cols_);
  }
  // (i,j)番目の要素
  double operator() (int i, int j) const {
    return lm_(i+rowoffset_, j+coloffset_);
  }
  double &operator() (int i, int j) {
    return lm_(i+rowoffset_, j+coloffset_);
  }
  // 代入（ = 行列）
  dSubMatrix &operator = (const dMatrix<Rows_,Cols_> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] = dm[i];
    }
    return *this;
  }
  dSubMatrix &operator = (const dIdentityMatrix<Rows_,Cols_> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i].clear();
      (*this)(i,i) = 1;
    }
    return *this;
  }
  template<class E>
  dSubMatrix &operator = (const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
  	(*this)(i,j) = dm()(i,j);
      }
    }
    return *this;
  }
  // 代入（ = ２次元配列）
  dSubMatrix &operator = (const double dpp[Rows_][Cols_]) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] = dpp[i];
    }
    return *this;
  }
  // 代入（ = 数値(*単位行列)）
  dSubMatrix &operator = (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i].clear();
      (*this)(i,i) = d;
    }
    return *this;
  }
  // 代入（ += 行列）
  template<class E>
  dSubMatrix &operator += (const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
  	(*this)(i,j) += dm()(i,j);
      }
    }
    return *this;
  }
  // 代入（ -= 行列）
  template<class E>
  dSubMatrix &operator -= (const dMatrixExpression<E> &dm) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      #pragma ivdep
      for (int j = 0; j < Cols_; j++) {
  	(*this)(i,j) -= dm()(i,j);
      }
    }
    return *this;
  }
  // 代入（ *= 行列）
  template<class E>
  dSubMatrix &operator *= (const dMatrixExpression<E> &dm) {
    dMatrix<Rows_, Cols_> result = (*this) * dm;
    return *this = result; //非効率だが上書き計算できないのでこうする
  }
  // 代入（ *= 数値）
  dSubMatrix &operator *= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] *= d;
    }
    return *this;
  }
  // 代入（ /= 数値）
  dSubMatrix &operator /= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] /= d;
    }
    return *this;
  }
};

// 部分行列を返す submatrix(行列，行オフセット，列オフセット，行サイズ:定数，列サイズ:定数)
#define submatrix(lm, rowoffset, coloffset, rows, cols)	\
  dSubMatrix<typeof((lm)), (rows), (cols)>((lm), (rowoffset), (coloffset))
//関数ではテンプレート引数を与えられないのでマクロで記述


#include "MatrixExpression.h"

#endif
