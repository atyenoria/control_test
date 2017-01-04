/* Time-stamp: <2012-12-23 19:04:43 kobayasi>
   
   �٥��ȥ륯�饹
   ver.2012-12-23  by K.Kobayashi
   
   ��®�ʻ�§�黻���Ǥ��룱��������ʥ٥��ȥ��
   
   ��Ȥ�����
   dVector<2> v;  // �����Ĺ��2��double���٥��ȥ�
   
   v = 1.0;       // ���٤Ƥ����Ǥ�������
   v[0] = 1.0;    // 0���ܤ����Ǥ�������v(0) = 1.0 �Ǥ�Ʊ����
   
   ����ա�
   ���黻�κݤ˥٥��ȥ�Υ����������������ɤ�����Ƚ�Ǥ��ʤ��ʹ�®���Τ����
   ���黻���Ƥˤ�ä���٥��ȥ���ä���ԥ٥��ȥ���ä��ꤹ��
   ��(��)�٥��ȥ�*���󡤹���*(��)�٥��ȥ�η�̤�黻�˻Ȥä��ѿ��������Ǥ��ʤ�

   �Ȥ���黻�� VectorExpression.h �򻲾ȡ�
*/

#ifndef __VECTORECLASS__
#define __VECTORECLASS__

#include <math.h>
#include <string.h>
#include <float.h>

// �٥��ȥ��ɽ�����饹(�������)
template<class E>
class dVectorExpression;


// �٥��ȥ륯�饹
template <int Size_>
class dVector:
  public dVectorExpression< dVector<Size_> > {
public:
  static const bool Costs = 0;
  static const int Size = Size_;
private:
  double data_[Size_];//  __attribute__((aligned(8)));
public:
  // ���󥹥ȥ饯���ʰ����ʤ���
  dVector() {}
  
  // ���ԡ����󥹥ȥ饯���ʥ٥��ȥ�򥳥ԡ���
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
  // ���ԡ����󥹥ȥ饯��������򥳥ԡ���
  dVector(const double *dp) {
    memcpy(data_, dp, sizeof(double)*Size_);
  }
  // ���ԡ����󥹥ȥ饯���ʿ��ͤ������
  dVector(const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = d;
    }
  }
  
  /** function *************************************************************/
  
  // �٥��ȥ��Ĺ��
  int size() const {
    return Size_;
  }
  // �٥��ȥ�򥯥ꥢ
  void clear() {
    memset(data_, 0, sizeof(double)*Size_);
  }
  // �����Σ�����������֤�
  double *data() {
    return data_;
  }
  double *data() const {
    return data_;
  }
  
  /** operator *************************************************************/
  
  // i���ܤ�����
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
  // ������ = �٥��ȥ��
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
  // ������ = �����������
  dVector &operator = (const double *dp) {
    memcpy(data_, dp, sizeof(double)*Size_);
    return *this;
  }
  // ������ = ���͡�
  dVector &operator = (const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = d;
    }
    return *this;
  }
  // ������ += �٥��ȥ��
  template<class E>
  dVector &operator += (const dVectorExpression<E> &dv) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] += dv()[i];
    }
    return *this;
  }
  // ������ -= �٥��ȥ��
  template<class E>
  dVector &operator -= (const dVectorExpression<E> &dv) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] -= dv()[i];
    }
    return *this;
  }
  // ������ *= ���͡�
  dVector &operator *= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] *= d;
    }
    return *this;
  }
  // ������ /= ���͡�
  dVector &operator /= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] /= d;
    }
    return *this;
  }
};


// ɸ�����٥��ȥ�Υ��饹���������������ʤ��Τǥ��������
template <int Size_>
class dStdBasisVector:
  public dVectorExpression< dStdBasisVector<Size_> > {
public:
  static const bool Costs = 0;
  static const int Size = Size_;
private:
  const int index_;
public:
  // ���󥹥ȥ饯���ʰ���������Υ���ǥå�����
  dStdBasisVector(const int i): index_(i) {}
  
  /** function *************************************************************/
  
  // �٥��ȥ��Ĺ��
  int size() const {
    return Size_;
  }
  
  /** operator *************************************************************/
  
  // i���ܤ�����
  double operator[] (int i) const {
    return (i == index_)? 1 : 0;
  }
  double operator() (int i) const {
    return (i == index_)? 1 : 0;
  }
};


// ��ʬ�٥��ȥ��ɽ�����饹
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
  
  // �٥��ȥ��Ĺ��
  int size() const {
    return Size_;
  }
  // �٥��ȥ�򥯥ꥢ
  void clear() {
    memset(data(), 0, sizeof(double)*Size_);
  }
  // �����Σ�����������֤�
  double *data() {
    return lv_.data()+offset_;
  }
  double *data() const {
    return lv_.data()+offset_;
  }
  
  /** operator *************************************************************/
  
  // i���ܤ�����
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
  // ������ = �٥��ȥ��
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
  // ������ = �����������
  dSubVector &operator = (const double *dp) {
    memcpy(data(), dp, sizeof(double)*Size_);
    return *this;
  }
  // ������ = ���͡�
  dSubVector &operator = (const double d) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] = d;
    }
    return *this;
  }
  // ������ += �٥��ȥ��
  template<class E>
  dSubVector &operator += (const dVectorExpression<E> &dv) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] += dv()[i];
    }
    return *this;
  }
  // ������ -= �٥��ȥ��
  template<class E>
  dSubVector &operator -= (const dVectorExpression<E> &dv) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] -= dv()[i];
    }
    return *this;
  }
  // ������ *= ���͡�
  dSubVector &operator *= (const double d) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] *= d;
    }
    return *this;
  }
  // ������ /= ���͡�
  dSubVector &operator /= (const double d) {
    double *data_ = data();
    #pragma ivdep
    for (int i = 0; i < Size_; i++) {
      data_[i] /= d;
    }
    return *this;
  }
};

// ��ʬ�٥��ȥ���֤� subvector(�٥��ȥ롤���ե��åȡ�������:���)
#define subvector(lv, offset, size) \
  dSubVector<typeof((lv)), (size)>((lv), (offset))
//�ؿ��Ǥϥƥ�ץ졼�Ȱ�����Ϳ�����ʤ��Τǥޥ���ǵ���


#include "VectorExpression.h"

#endif
