/* Time-stamp: <2012-12-26 10:56:05 kobayasi>
   
   ���󥯥饹
   ver.2012-12-23  by K.Kobayashi
   
   ��®�ʻ�§�黻���Ǥ��룲��������ʹ����
   
   ��Ȥ�����
   dMatrix<3,3> m;  // �����3��x3���double������
   
   m = 1.0;         // ���٤Ƥ��г���ʬ��������¾��0��
                    // ��ñ�̹��󤬱���Ƥ�������� m=1.0*I ��
   m(1,2) = 1.0;    // (1,2)�����Ǥ�������m[1][2] = 1.0 �Ǥ��̤�Ʊ������
                    // ���ԥ٥��ȥ뤫���ͤ��뤳�Ȥˤʤ뤿�ᡤm(1,2)������®����
   d = (m * m)(1,2);// �黻����ľ��(1,2)���뤳�Ȥ�OK
                    // ��(1,2)��ɬ�פ�ʬ�η׻������Ԥ�ʤ��Τǹ�®��
		    
   ����ա�
   ���黻�κݤ˹����٥��ȥ�Υ����������������ɤ�����Ƚ�Ǥ��ʤ��ʹ�®���Τ����
   ���黻���Ƥˤ�ä���٥��ȥ���ä���ԥ٥��ȥ���ä��ꤹ��
   ������*����(��)�٥��ȥ�*���󡤹���*(��)�٥��ȥ�η�̤�黻�˻Ȥä��ѿ��������Ǥ��ʤ�
   
   �Ȥ���黻�� MatrixExpression.h �򻲾ȡ�
*/

#ifndef __MATRIXCLASS__
#define __MATRIXCLASS__

#include "VectorClass.h"

// �����ɽ�����饹(�������)
template<class E>
class dMatrixExpression;

// ñ�̹��󥯥饹(�������)
template <int Rows_, int Cols_>
class dIdentityMatrix;


// ���󥯥饹
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
  // ���󥹥ȥ饯���ʰ����ʤ���
  dMatrix() {}
  
  // ���ԡ����󥹥ȥ饯���ʹ���򥳥ԡ���
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
  // ���ԡ����󥹥ȥ饯����2��������򥳥ԡ���
  dMatrix(const double dpp[Rows_][Cols_]) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] = dpp[i];
    }
  }
  // ���ԡ����󥹥ȥ饯���ʿ���(*ñ�̹���)�������
  dMatrix(const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
      data_[i][i] = d;
    }
  }
  
  /** function *************************************************************/
  
  // ����ιԿ�
  int rows() const {
    return Rows_;
  }
  // ��������
  int cols() const {
    return Cols_;
  }
  // ����򥯥ꥢ
  void clear() {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
    }
  }
  // �����Υ٥��ȥ�������֤�
  rowvector_type *data() {
    return data_;
  }
  
  /** operator *************************************************************/
  
  // i�ԤΥ٥��ȥ�
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
  // (i,j)���ܤ�����
  double operator() (int i, int j) const {
    return data_[i][j];
  }
  double &operator() (int i, int j) {
    return data_[i][j];
  }
  // ������ = �����
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
  // ������ = �����������
  dMatrix &operator = (const double dpp[Rows_][Cols_]) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] = dpp[i];
    }
    return *this;
  }
  // ������ = ����(*ñ�̹���)��
  dMatrix &operator = (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i].clear();
      data_[i][i] = d;
    }
    return *this;
  }
  // ������ += �����
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
  // ������ -= �����
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
  // ������ *= �����
  template<class E>
  dMatrix &operator *= (const dMatrixExpression<E> &dm) {
    dMatrix<Rows_, Cols_> result = (*this) * dm;
    return *this = result; //���Ψ������񤭷׻��Ǥ��ʤ��ΤǤ�������
  }
  // ������ *= ���͡�
  dMatrix &operator *= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] *= d;
    }
    return *this;
  }
  // ������ /= ���͡�
  dMatrix &operator /= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      data_[i] /= d;
    }
    return *this;
  }
};


// ñ�̹���Υ��饹���������������ʤ��Τǥ��������
template <int Rows_, int Cols_>
class dIdentityMatrix:
  public dMatrixExpression< dIdentityMatrix<Rows_, Cols_> > {
public:
  static const bool Costs = 0;
  static const int Rows = Rows_;
  static const int Cols = Cols_;
  typedef dStdBasisVector<Cols_> rowvector_type;
public:
  // ���󥹥ȥ饯���ʰ����ʤ���
  dIdentityMatrix() {}
  
  /** function *************************************************************/
  
  // ����ιԿ�
  int rows() const {
    return Rows_;
  }
  // ��������
  int cols() const {
    return Cols_;
  }
  
  /** operator *************************************************************/
  
  // i�ԤΥ٥��ȥ�
  rowvector_type operator[] (int i) const {
    return rowvector_type(i);
  }
  rowvector_type operator() (int i) const {
    return rowvector_type(i);
  }
  // (i,j)���ܤ�����
  double operator() (int i, int j) const {
    return (i == j);
  }
};


// ��ʬ�����ɽ�����饹
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
  dSubMatrix(const L &lm, const int rowoffset, const int coloffset): lm_(*const_cast<L*>(&lm)), rowoffset_(rowoffset), coloffset_(coloffset) {}  //�����const��Ϥ���
  
  /** function *************************************************************/
  
  // ����ιԿ�
  int rows() const {
    return Rows_;
  }
  // ��������
  int cols() const {
    return Cols_;
  }
  // ����򥯥ꥢ
  void clear() {
    #pragma ivdep
    for (int i = rowoffset_; i < rowoffset_+Rows_; i++) {
      memset(lm_[i].data()+coloffset_, 0, sizeof(double)*Cols_);
    }
  }
  
  /** operator *************************************************************/
  
  // i�ԤΥ٥��ȥ�
  rowvector_type operator[] (int i) const {
    return subvector(lm_[i+rowoffset_], coloffset_, Cols_);
  }
  rowvector_type operator() (int i) const {
    return subvector(lm_[i+rowoffset_], coloffset_, Cols_);
  }
  // (i,j)���ܤ�����
  double operator() (int i, int j) const {
    return lm_(i+rowoffset_, j+coloffset_);
  }
  double &operator() (int i, int j) {
    return lm_(i+rowoffset_, j+coloffset_);
  }
  // ������ = �����
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
  // ������ = �����������
  dSubMatrix &operator = (const double dpp[Rows_][Cols_]) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] = dpp[i];
    }
    return *this;
  }
  // ������ = ����(*ñ�̹���)��
  dSubMatrix &operator = (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i].clear();
      (*this)(i,i) = d;
    }
    return *this;
  }
  // ������ += �����
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
  // ������ -= �����
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
  // ������ *= �����
  template<class E>
  dSubMatrix &operator *= (const dMatrixExpression<E> &dm) {
    dMatrix<Rows_, Cols_> result = (*this) * dm;
    return *this = result; //���Ψ������񤭷׻��Ǥ��ʤ��ΤǤ�������
  }
  // ������ *= ���͡�
  dSubMatrix &operator *= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] *= d;
    }
    return *this;
  }
  // ������ /= ���͡�
  dSubMatrix &operator /= (const double d) {
    #pragma ivdep
    for (int i = 0; i < Rows_; i++) {
      (*this)[i] /= d;
    }
    return *this;
  }
};

// ��ʬ������֤� submatrix(���󡤹ԥ��ե��åȡ��󥪥ե��åȡ��ԥ�����:������󥵥���:���)
#define submatrix(lm, rowoffset, coloffset, rows, cols)	\
  dSubMatrix<typeof((lm)), (rows), (cols)>((lm), (rowoffset), (coloffset))
//�ؿ��Ǥϥƥ�ץ졼�Ȱ�����Ϳ�����ʤ��Τǥޥ���ǵ���


#include "MatrixExpression.h"

#endif
