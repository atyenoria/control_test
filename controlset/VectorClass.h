/* Time-stamp: <2012-12-23 19:04:43 kobayasi>
   
   ベクトルクラス
   ver.2012-12-23  by K.Kobayashi
   
   高速な四則演算ができる１次元配列（ベクトル）
   
   ＜使い方＞
   dVector<2> v;  // 定義：長さ2のdouble型ベクトル
   
   v = 1.0;       // すべての要素に代入．
   v[0] = 1.0;    // 0番目の要素に代入．v(0) = 1.0 でも同じ．
   
   ＜注意＞
   ・演算の際にベクトルのサイズが正しいかどうかは判断しない（高速化のため）
   ・演算内容によって列ベクトルだったり行ベクトルだったりする
   ・(行)ベクトル*行列，行列*(列)ベクトルの結果を演算に使った変数に代入できない

   使える演算は VectorExpression.h を参照．
*/

#ifndef __VECTORECLASS__
#define __VECTORECLASS__

#include <math.h>
#include <string.h>
#include <float.h>

// ベクトルを表すクラス(前方宣言)
template<class E>
class dVectorExpression;


// ベクトルクラス
template <int Size_>
class dVector:
  public dVectorExpression< dVector<Size_> > {
public:
  static const bool Costs = 0;
  static const int Size = Size_;
private:
  double data_[Size_];//  __attribute__((aligned(8)));
public:
  // コンストラクタ（引数なし）
  dVector() {}
  
  // コピーコンストラクタ（ベクトルをコピー）
  dVector(const dVector &dv) {
    memcpy(data_, dv.data_, sizeof(double)*Size_);
  }
  template<class E>
  dVector(const dVectorExpression<E> &dv) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = dv()[i];
    }
  }
  // コピーコンストラクタ（配列をコピー）
  dVector(const double *dp) {
    memcpy(data_, dp, sizeof(double)*Size_);
  }
  // コピーコンストラクタ（数値で埋める）
  dVector(const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = d;
    }
  }
  
  /** function *************************************************************/
  
  // ベクトルの長さ
  int size() const {
    return Size_;
  }
  // ベクトルをクリア
  void clear() {
    memset(data_, 0, sizeof(double)*Size_);
  }
  // 内部の１次元配列を返す
  double *data() {
    return data_;
  }
  double *data() const {
    return data_;
  }
  
  /** operator *************************************************************/
  
  // i番目の要素
  double operator[] (int i) const {
    return data_[i];
  }
  double &operator[] (int i) {
    return data_[i];
  }
  double operator() (int i) const {
    return data_[i];
  }
  double &operator() (int i) {
    return data_[i];
  }
  // 代入（ = ベクトル）
  dVector &operator = (const dVector &dv) {
    memcpy(data_, dv.data_, sizeof(double)*Size_);
    return *this;
  }
  template<class E>
  dVector &operator = (const dVectorExpression<E> &dv) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = dv()[i];
    }
    return *this;
  }
  // 代入（ = １次元配列）
  dVector &operator = (const double *dp) {
    memcpy(data_, dp, sizeof(double)*Size_);
    return *this;
  }
  // 代入（ = 数値）
  dVector &operator = (const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = d;
    }
    return *this;
  }
  // 代入（ += ベクトル）
  template<class E>
  dVector &operator += (const dVectorExpression<E> &dv) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] += dv()[i];
    }
    return *this;
  }
  // 代入（ -= ベクトル）
  template<class E>
  dVector &operator -= (const dVectorExpression<E> &dv) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] -= dv()[i];
    }
    return *this;
  }
  // 代入（ *= 数値）
  dVector &operator *= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] *= d;
    }
    return *this;
  }
  // 代入（ /= 数値）
  dVector &operator /= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] /= d;
    }
    return *this;
  }
};


// 標準基底ベクトルのクラス（内部配列を持たないのでメモリ節約）
template <int Size_>
class dStdBasisVector:
  public dVectorExpression< dStdBasisVector<Size_> > {
public:
  static const bool Costs = 0;
  static const int Size = Size_;
private:
  const int index_;
public:
  // コンストラクタ（引数：基底のインデックス）
  dStdBasisVector(const int i): index_(i) {}
  
  /** function *************************************************************/
  
  // ベクトルの長さ
  int size() const {
    return Size_;
  }
  
  /** operator *************************************************************/
  
  // i番目の要素
  double operator[] (int i) const {
    return (i == index_)? 1 : 0;
  }
  double operator() (int i) const {
    return (i == index_)? 1 : 0;
  }
};


// 部分ベクトルを表すクラス
template <class L, int Size_>
class dSubVector:
  public dVectorExpression< dSubVector<L, Size_> > {
public:
  static const bool Costs = 0;
  static const int Size = Size_;
private:
  L &lv_;
  const int offset_;
public:
  //dSubVector(L &lv, const int offset): lv_(&lv), offset_(offset) {}
  dSubVector(const L &lv, const int offset): lv_(*const_cast<L*>(&lv)), offset_(offset) {}
  
  /** function *************************************************************/
  
  // ベクトルの長さ
  int size() const {
    return Size_;
  }
  // ベクトルをクリア
  void clear() {
    memset(data(), 0, sizeof(double)*Size_);
  }
  // 内部の１次元配列を返す
  double *data() {
    return lv_.data()+offset_;
  }
  double *data() const {
    return lv_.data()+offset_;
  }
  
  /** operator *************************************************************/
  
  // i番目の要素
  double operator[] (int i) const {
    return lv_[i+offset_];
  }
  double &operator[] (int i) {
    return lv_[i+offset_];
  }
  double operator() (int i) const {
    return lv_[i+offset_];
  }
  double &operator() (int i) {
    return lv_[i+offset_];
  }
  // 代入（ = ベクトル）
  dSubVector &operator = (dVector<Size_> &dv) {
    memcpy(data(), dv.data(), sizeof(double)*Size_);
    return *this;
  }
  template<class E>
  dSubVector &operator = (const dVectorExpression<E> &dv) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = dv()[i];
    }
    return *this;
  }
  // 代入（ = １次元配列）
  dSubVector &operator = (const double *dp) {
    memcpy(data(), dp, sizeof(double)*Size_);
    return *this;
  }
  // 代入（ = 数値）
  dSubVector &operator = (const double d) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = d;
    }
    return *this;
  }
  // 代入（ += ベクトル）
  template<class E>
  dSubVector &operator += (const dVectorExpression<E> &dv) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] += dv()[i];
    }
    return *this;
  }
  // 代入（ -= ベクトル）
  template<class E>
  dSubVector &operator -= (const dVectorExpression<E> &dv) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] -= dv()[i];
    }
    return *this;
  }
  // 代入（ *= 数値）
  dSubVector &operator *= (const double d) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] *= d;
    }
    return *this;
  }
  // 代入（ /= 数値）
  dSubVector &operator /= (const double d) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] /= d;
    }
    return *this;
  }
};

// 部分ベクトルを返す subvector(ベクトル，オフセット，サイズ:定数)
#define subvector(lv, offset, size) \
  dSubVector<typeof((lv)), (size)>((lv), (offset))
//関数ではテンプレート引数を与えられないのでマクロで記述


#include "VectorExpression.h"

#endif
