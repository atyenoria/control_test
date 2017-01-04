/* Time-stamp: <2013-04-26 13:38:54 kobayasi>
 
   ����ι�®�黻�Τ����ExpressionTemplate��MatrixClass.h�������ѡ�
   ver.2013-04-26  by K.Kobayashi
   
   fabs(A) �������Ǥ������ͤ�Ȥä�������֤�
   round(A)�����������ͼθ�������
   ceil(A) ����������ŷ��ؿ�����
   floor(A)�������������ؿ�������
   square(A)�� ��������衡������   �ʢ���������ǤϤʤ��Τ���ա�
   
   submatrix(A,i,j,N,M)������A��(i,j)������ʬ����<N,M>���֤�
                        ��A��ľ�ܽ񤭴�����ǽ��submatrix(A,...)=X��
   
   rowvector(A,i)��i�ԥ٥��ȥ���֤�
   colvector(A,j)��j��٥��ȥ���֤�
   
   diag(A)  ���г����ǥ٥��ȥ���֤�
   diag(B)  ���٥��ȥ���г����Ǥ˻��Ĺ�����֤�
   trans(A) ��ž�ֹ�����֤�
   trans(B) ���٥��ȥ�������<1,Size>�Ȥ����֤�
   matrix(B)���٥��ȥ��Թ���<Size,1>�Ȥ����֤�
   matmul(B1,B2)����٥��ȥ�*�ԥ٥��ȥ�=������֤���=matrix(B1)*trans(B2)��
   
   norm1(A)����1-�Υ����֤�
   infnorm(A)��̵����Υ����֤�
   fnorm(A)�����ե�٥˥����Υ����֤�
   max(A)    ����������Ǥ��֤�
   min(A)    ���Ǿ������Ǥ��֤�
   trace(A)  ���г���ʬ���¤��֤�
   
   chol(A,L)   ��Choleskyʬ�򤹤�ʹ���L�˽񤭹��ߡ�
   lu(A,p,LU)  ��LUʬ�򤹤������p,����LU�˽񤭹��ߡ�
   det(A)      �����󼰤��֤�
   solve(A,x=y)��ϢΩ�켡��������A*x=y�ˤ�򤯡ʹ���x�˽񤭹��ߡ�
   inv(A,iA)   ���չ���A^-1�����ʹ���iA�˽񤭹��ߡ�
   inv(A)      ��    ��    ���֤�
   expm(A,eA)  ������λؿ��ؿ�e^A�����ʹ���eA�˽񤭹��ߡ�
   expm(A)     ��    ��           ���֤�
   signm(A,S)  ����������ؿ������ʹ���S�˽񤭹��ߡ�
   signm(A)    ��    ��        ���֤�
   sqrtm(A,sA) �������ʿ������A����ʹ���eA�˽񤭹��ߡ�
   sqrtm(A)    ��    ��        ���֤�
   qr(A,Q,R)   ��QRʬ�򤹤�ʹ���Q,R�˽񤭹��ߡ�
   eig(A,er,ei)����ͭ�ͤ����ʸ�ͭ�ͤΥ٥��ȥ�����er,����ei�˽񤭹��ߡ�
   
   print(A)������ printf��Ȥäƹ����ɽ������
   printf(fmt,A)�������ͤ�ɽ���ե����ޥåȻ�����
*/


// �����ɽ�����饹
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


// �黻�����Ȥ˱����Ʊ黻��̤򥭥�å��夹�빽¤��
// �ʹ���*����Τ褦��ʣ�����Ǥ��Ϥ�黻������Ȥ��˥���å��塧Costs=1��
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


// ����α黻����ɽ�����饹�ʹ���ȹ���βø�����
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

// ����α黻����ɽ�����饹�ʹ���ȿ��ͤξ������
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
  const double rd_;//  __attribute__((aligned(8))); //®���ʤ롩
public:
  dMatrixOperation(const L &lm, const double rd): lm_(lm), rd_(rd) {}
  
  rowvector_type operator[] (int i) const { return rowvector_type(lm_[i], rd_); }
  rowvector_type operator() (int i) const { return rowvector_type(lm_(i), rd_); }
  double operator() (int i, int j)  const { return Ope::apply(lm_(i,j), rd_); }
};

// ����α黻����ɽ�����饹�ʹ������黻�ʤɡ�
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

// �û��ʹ��� + �����
template<class L, class R>
inline
dMatrixOperation<L, OperationAdd, R> operator + (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  return dMatrixOperation<L, OperationAdd, R>( lm(), rm() );
}
// �����ʹ��� - �����
template<class L, class R>
inline
dMatrixOperation<L, OperationSub, R> operator - (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  return dMatrixOperation<L, OperationSub, R>( lm(), rm() );
}
// �軻�ʹ��� * ���͡�
template<class L>
inline
dMatrixOperation<L, OperationMul, double> operator * (const dMatrixExpression<L> &lm, const double rd) {
  return dMatrixOperation<L, OperationMul, double>( lm(), rd );
}
// �軻�ʿ��� * ����ˡ�ʹ��� * ���͡�
template<class R>
inline
dMatrixOperation<R, OperationMul, double> operator * (const double ld, const dMatrixExpression<R> &rm) {
  return dMatrixOperation<R, OperationMul, double>( rm(), ld );
}
// �����ʹ��� / ���͡�
template<class L>
inline
dMatrixOperation<L, OperationDiv, double> operator / (const dMatrixExpression<L> &lm, const double rd) {
  return dMatrixOperation<L, OperationDiv, double>( lm(), rd );
}
// ���- �����
template<class R>
inline
dMatrixOperationS<OperationNeg, R> operator - (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationNeg, R>( rm() );
}
// ����ʹ��� == �����// ���Ǥ�double���ʤΤǱ黻��̤���Ӥ���ΤϤ������ᤷ�ʤ�
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
// ������ʹ��� != �����
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

// ������ fabs(����)
template<class R>
inline
dMatrixOperationS<OperationFabs, R> fabs (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationFabs, R>( rm() );
}
// �ͼθ��� round(����)
template<class R>
inline
dMatrixOperationS<OperationRound, R> round (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationRound, R>( rm() );
}
// ŷ��ؿ� ceil(����)
template<class R>
inline
dMatrixOperationS<OperationCeil, R> ceil (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationCeil, R>( rm() );
}
// ���ؿ� floor(����)
template<class R>
inline
dMatrixOperationS<OperationFloor, R> floor (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationFloor, R>( rm() );
}
// �����Ǥ���� square(����)
template<class R>
inline
dMatrixOperationS<OperationSquare, R> square (const dMatrixExpression<R> &rm) {
  return dMatrixOperationS<OperationSquare, R>( rm() );
}


// ����α黻����ɽ�����饹�ʹ���ιԥ٥��ȥ��
template<class L>
class dMatrixOpeRow:
  public dVectorExpression< dMatrixOpeRow<L> > {
public:
  static const bool Costs = L::Costs;
  static const int Size = L::Cols; //Ĺ�������
private:
  const L &lm_;
  const int ri_;
public:
  dMatrixOpeRow(const L &lm, const int ri): lm_(lm), ri_(ri) {}
  
  double operator[] (int j) const { return lm_(ri_,j); }
  double operator() (int j) const { return lm_(ri_,j); }
};

// ����α黻����ɽ�����饹�ʹ������٥��ȥ��
template<class L>
class dMatrixOpeCol:
  public dVectorExpression< dMatrixOpeCol<L> > {
public:
  static const bool Costs = L::Costs;
  static const int Size = L::Rows; //Ĺ���ϹԿ�
private:
  const L &lm_;
  const int rj_;
public:
  dMatrixOpeCol(const L &lm, const int rj): lm_(lm), rj_(rj) {}
  
  double operator[] (int i) const { return lm_(i,rj_); }
  double operator() (int i) const { return lm_(i,rj_); }
};

// i�ԥ٥��ȥ� row(����,���ֹ�)
template<class L>
inline
dMatrixOpeRow<L> rowvector (const dMatrixExpression<L> &lm, int i) {
  return dMatrixOpeRow<L>( lm(), i );
}

// j��٥��ȥ� col(����,���ֹ�)
template<class L>
inline
dMatrixOpeCol<L> colvector (const dMatrixExpression<L> &lm, int j) {
  return dMatrixOpeCol<L>( lm(), j );
}


// ����*��٥��ȥ�α黻����ɽ�����饹
template<class L, class R>
class dMatrixOpeMulL:
  public dVectorExpression< dMatrixOpeMulL<L, R> > {
public:
  static const bool Costs = 1;
  static const int Size = L::Rows; //Ĺ���Ϲ���L�ιԿ��ˤʤ�
private:
  const dMatrixCache<L,L::Costs> lm_; //����å���ǹ�®��
  const dVectorCache<R,R::Costs> rv_; //����å���ǹ�®��
  //const L &lm_;
  //const R &rv_;
public:
  dMatrixOpeMulL(const L &lm, const R &rv): lm_(dMatrixCache<L,L::Costs>(lm)), rv_(dVectorCache<R,R::Costs>(rv)) {}
  //dMatrixOpeMulL(const L &lm, const R &rv): lm_(lm), rv_(rv) {}
  
  double operator[] (int i) const {
    /*
    double result = 0;
    for (int j = 0; j < L::Cols; j++) {
      result += lm_(i,j) * rv_[j]; //i�֤����ǤϹ���L��i��*��٥��ȥ�
    }
    return result;
    */
    return rowvector(lm_,i) * rv_;
  }
  double operator() (int i) const {
    return rowvector(lm_,i) * rv_;
  }
};

// �ԥ٥��ȥ�*����α黻����ɽ�����饹
template<class L, class R>
class dMatrixOpeMulR:
  public dVectorExpression< dMatrixOpeMulR<L, R> > {
public:
  static const bool Costs = 1;
  static const int Size = R::Cols; //Ĺ���Ϲ���R������ˤʤ�
private:
  const dVectorCache<L,L::Costs> lv_; //����å���ǹ�®��
  const dMatrixCache<R,R::Costs> rm_; //����å���ǹ�®��
  //const L &lv_;
  //const R &rm_;
public:
  dMatrixOpeMulR(const L &lv, const R &rm): lv_(dVectorCache<L,L::Costs>(lv)), rm_(dMatrixCache<R,R::Costs>(rm)) {}
  //dMatrixOpeMulR(const L &lv, const R &rm): lv_(lv), rm_(rm) {}
  
  double operator[] (int j) const {
    /*
    double result = 0;
    for (int i = 0; i < R::Rows; i++) {
      result += lv_[i] * rm_(i,j); //j�֤����ǤϹԥ٥��ȥ�*����R��j��
    }
    return result;
    */
    return lv_ * colvector(rm_,j);
  }
  double operator() (int j) const {
    return lv_ * colvector(rm_,j);
  }
};

// ����*����α黻����ɽ�����饹
template<class L, class R>
class dMatrixOpeMul:
  public dMatrixExpression< dMatrixOpeMul<L, R> > {
public:
  static const bool Costs = 1;
  static const int Rows = L::Rows; //�Կ��Ϲ���L�ιԿ��ˤʤ�
  static const int Cols = R::Cols; //����Ϲ���R������ˤʤ�
  typedef dMatrixOpeMulR<typename L::rowvector_type, R> rowvector_type;
private:
  const dMatrixCache<L,L::Costs> lm_; //����å���ǹ�®��
  const dMatrixCache<R,R::Costs> rm_; //����å���ǹ�®��
  //const L &lm_;
  //const R &rm_;
public:
  dMatrixOpeMul(const L &lm, const R &rm): lm_(dMatrixCache<L,L::Costs>(lm)), rm_(dMatrixCache<R,R::Costs>(rm)) {}
  //dMatrixOpeMul(const L &lm, const R &rm): lm_(lm), rm_(rm) {}
  
  //�ԥ٥��ȥ�*����
  rowvector_type operator[] (int i) const { return rowvector_type(lm_[i], rm_); }
  rowvector_type operator() (int i) const { return rowvector_type(lm_(i), rm_); }
  
  double operator() (int i, int j) const {
    /*
    double result = 0;
    for (int k = 0; k < L::Cols; k++) {
      result += lm_(i,k) * rm_(k,j); //(i,j)�֤����ǤϹ���L��i��*����R��j��
    }
    return result;
    */
    return rowvector(lm_,i) * colvector(rm_,j);
  }
};

// �軻�ʹ��� * ��٥��ȥ�ˡʾ������٥��ȥ�Ȳ���
template<class L, class R>
inline
dMatrixOpeMulL<L, R> operator * (const dMatrixExpression<L> &lm, const dVectorExpression<R> &rv) {
  return dMatrixOpeMulL<L, R>( lm(), rv() );
}
// �軻�ʹԥ٥��ȥ� * ����ˡʾ���˹ԥ٥��ȥ�Ȳ���
template<class L, class R>
inline
dMatrixOpeMulR<L, R> operator * (const dVectorExpression<L> &lv, const dMatrixExpression<R> &rm) {
  return dMatrixOpeMulR<L, R>( lv(), rm() );
}
// �軻�ʹ��� * �����
template<class L, class R>
inline
dMatrixOpeMul<L, R> operator * (const dMatrixExpression<L> &lm, const dMatrixExpression<R> &rm) {
  return dMatrixOpeMul<L, R>( lm(), rm() );
}


// ����α黻����ɽ�����饹�ʹ�����г����ǥ٥��ȥ��
template<class R>
class dMatrixOpeDiag:
  public dVectorExpression< dMatrixOpeDiag<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Size = R::Rows; //Ĺ���ϹԿ�
private:
  const R &rm_;
public:
  dMatrixOpeDiag(const R &rm): rm_(rm) {}
  
  double operator[] (int i) const { return rm_(i,i); }
  double operator() (int i) const { return rm_(i,i); }
};

// �٥��ȥ�α黻����ɽ�����饹�ʥ٥��ȥ���г����Ǥ˻��Ĺ����
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

// �г����ǥ٥��ȥ� diag(����)
template<class R>
inline
dMatrixOpeDiag<R> diag (const dMatrixExpression<R> &rm) {
  return dMatrixOpeDiag<R>( rm() );
}

// �٥��ȥ���г����Ǥ˻��Ĺ��� diag(�٥��ȥ�)
template<class R>
inline
dVectorOpeDiag<R> diag (const dVectorExpression<R> &rv) {
  return dVectorOpeDiag<R>( rv() );
}


// ����α黻����ɽ�����饹�ʹ����ž�֡�
template<class R>
class dMatrixOpeTrans:
  public dMatrixExpression< dMatrixOpeTrans<R> > {
public:
  static const bool Costs = R::Costs;
  static const int Rows = R::Cols; //ž�֤ʤΤǵդˤʤ�
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

// ž�� trans(����)
template<class R>
inline
dMatrixOpeTrans<R> trans (const dMatrixExpression<R> &rm) {
  return dMatrixOpeTrans<R>( rm() );
}

// �٥��ȥ�α黻����ɽ�����饹�ʥ٥��ȥ�����<1,Size>�Ȥߤʤ���
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
  
  rowvector_type operator[] (int i) const { return rv_; }//i=0�ʳ���Ʊ���ˤʤ�Τ����
  rowvector_type operator() (int i) const { return rv_; }
  double operator() (int i, int j) const { return rv_(j); }
};

// ž�֡ʥ٥��ȥ�����<1,Size>�Ȥߤʤ���trans(�٥��ȥ�)
template<class R>
inline
dVectorOpeTrans<R> trans (const dVectorExpression<R> &rv) {
  return dVectorOpeTrans<R>( rv() );
}

// �٥��ȥ�α黻����ɽ�����饹�ʥ٥��ȥ�����<Size,1>�Ȥߤʤ���
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

// �٥��ȥ�����<Size,1>�ȸ��ʤ� matrix(�٥��ȥ�)
template<class R>
inline
dVectorOpeMatrix<R> matrix (const dVectorExpression<R> &rv) {
  return dVectorOpeMatrix<R>( rv() );
}


// �٥��ȥ�α黻����ɽ�����饹����٥��ȥ�*�ԥ٥��ȥ뢪�����
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

// �軻����٥��ȥ�*�ԥ٥��ȥ뢪����ˡʾ������٥��ȥ�,�ԥ٥��ȥ�Ȳ���
template<class L, class R>
inline
dVectorOpeMul<L, R> matmul (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  return dVectorOpeMul<L, R>( lv(), rv() );
}


/**************************************************************************/


// 1-�Υ�� norm1(����)
template<class R>
inline
double norm1 (const dMatrixExpression<R> &rm) {
  // �������¡� max_j( sum_i( |rm(i,j)| ))
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

// ̵����Υ�� infnorm(����)
template<class R>
inline
double infnorm (const dMatrixExpression<R> &rm) {
  // ������¡� max_i( sum_j( |rm(i,j)| ))
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

// �ե�٥˥����Υ�� fnorm(����)
template<class R>
inline
double fnorm (const dMatrixExpression<R> &rm) {
  // �����ʿ������ ��sum_i( sum_j( rm(i,j)^2 ))
  int i, j;
  double sum=0;
  
  for(i=0; i<R::Rows; i++) {
    for(j=0; j<R::Cols; j++) {
      sum += OperationSquare::apply(rm()(i,j));
    }
  }
  return sqrt(sum);
}

// ������ max(����)
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

// �Ǿ��� min(����)
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

// �ȥ졼�����г���ʬ���¡� trance(����)
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

// Choleskyʬ�� chol(����,�����ѹ���) A = L*trans(L)
template <class E>
int chol(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &L) {
  int i, j, k;
  double sum;
  
  for (j=0; j<E::Rows; j++) {
    // L(j,j) �η׻�
    sum = A()(j,j);
    for (k=0; k<j; k++) {
      sum -= L(j,k) * L(j,k);
    }
    if (sum == 0) return 0;
    
    L(j,j) = sqrt(sum);
    
    // L(i,j) (i>j) �η׻�
    for (i=j+1; i<E::Rows; i++) {
      sum = A()(i,j);
      for (k=0; k<j; k++) {
        sum -= L(i,k) * L(j,k);
      }
      L(i,j) = sum / L(j,j);
      L(j,i) = 0; //�廰����ʬ��0
   }
  }
  return 1;
}

// LUʬ�� lu(����,�Ը򴹾��������,LU����) A = L*U
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
    for (j=0; j<E::Rows; j++) { //�����ͺ�������Ǥ����
      tmp = fabs(LU(i,j));
      if (w<tmp) w = tmp;
    }
    if (w==0) return 0; //singular
    
    weight[i] = 1/w; //�ԤνŤߤϺ��������ͤεտ�
    pivot[i] = i;    //�Ը򴹾���ν����
  }
  
  det = 1;
  for (i=0; i<E::Rows; i++) {
    w = -1;
    for (k=i; k<E::Rows; k++) { //��겼�γƹԤˤĤ���
      kk = pivot[k];            //�Ťߡ������ͤ�����ιԤ򸫤Ĥ���
      tmp = fabs(LU(kk,i)) * weight[kk];
      if (w<tmp) { w=tmp; j=k; } //j�ϽŤߡ������ͤ�����ι�
    }
    ii = pivot[j];
    if (j != i) {
      pivot[j] = pivot[i]; pivot[i] = ii; //���ֹ���
      det = -det; //�Ԥ�򴹤���й��󼰤���椬�Ѥ��
    }
    w = LU(ii,i); det *= w; //���󼰤�U���г���ʬ����
    if (w==0) return 0;
    
    for (k=i+1; k<E::Rows; k++) {  //Gauss�õ�ˡ
      kk = pivot[k];
      tmp = (LU(kk,i) /= w);      //�����ѹ���L�����Ǥ�׻�
      #pragma ivdep
      for (j=i+1; j<E::Rows; j++) {//�廰�ѹ���U�����Ǥ�׻�
	LU(kk,j) -= tmp * LU(ii,j);
      }
    }
  }
  return det; //����ͤϹ���
}

// LUʬ��ǹ��󼰤���� det(����)
template <class E>
inline
double det(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Cols> LU;
  int pivot[E::Rows]; 
  
  return lu(A, pivot, LU);
}

// LUʬ����Ѥ���ϢΩ�켡�������׻� solve(LU����,�Ը򴹾��������,��)
template <int Size_, int Num_>
inline
void solve(const dMatrix<Size_, Size_> &LU, const int pivot[], dMatrix<Size_,Num_> &solA) {
  //��A*x = y �Ȥ���ϢΩ�켡�������ϡ�L*U*x = y ��򤯤��Ȥ�����
  //  L*c = y �� U*x = c �ȷ׻����뤳�Ȥǵ������L,U�ϲ�����,�廰�ѤʤΤ�c,x�η׻���ñ��
  //������A��Ʊ��ʣ��(Num_��)��ϢΩ��������Ʊ���˲򤱤�
  //��x = A^-1 * y ��׻�����ΤϤ�����������չ���������®����������
  //  �չ�������Τ�y��ñ�̹���(�������Ը�ɬ��)�ˤ�����������
  
  int i, j, k, ii;
  double c[Size_];
  double sum;
  
  for (k=0; k<Num_; k++) {
    //��������
    for (i=0; i<Size_; i++) {
      ii = pivot[i];
      sum = solA(ii,k);
      for (j=i-1; j>=0; j--) {
	sum -= LU(ii,j) * c[j];
      }
      c[i] = sum;
    }
    //��������
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

// ϢΩ�켡��������A*x=y�ˤ�� solve(A,x = y) �� ��x=A^-1*y
template <class E, int Num_>
inline
int solve(const dMatrixExpression<E> &A, dMatrix<E::Rows,Num_> &solA) {
  dMatrix<E::Rows,E::Cols> LU;
  int pivot[E::Rows]; 
  
  if (lu(A, pivot, LU) == 0) return 0;
  
  solve(LU, pivot, solA);
  
  return 1;
}

// ϢΩ�켡��������A*x=y�ˤ�򤯡ʣ��ߣ��ѡ�
template <class E, int Num_>
inline
int solve(const dMatrixExpression<E> &A, dMatrix<2,Num_> &solA) {
  dMatrix<2,2> invA;
  dMatrix<2,Num_> tmp; //Ʊ���ѿ��������Ǥ��ʤ��Τ�������ͳ��
  
  //�չ�����᤿����®��
  if (inv(A, invA) == 0) return 0;
  
  tmp = invA*solA;
  solA = tmp;
  return 1;
}

// �չ��� inv(����,�չ���)
template <class E>
inline
double inv(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &invA) {
  int i, j, k, ii;
  double sum;
  
  dMatrix<E::Rows,E::Rows> LU;
  int pivot[E::Rows]; 
  
  double det = lu(A, pivot, LU);
  if (det == 0) return 0;
  
  //ñ�̹���ιԸ򴹤򤷤Ĥ�solve()
  for (k=0; k<E::Rows; k++) {
    //��������
    for (i=0; i<E::Rows; i++) {
      ii = pivot[i];
      sum = (ii == k);//�Ը򴹤��줿ñ�̹���
      for (j=i-1; j>=0; j--) {
	sum -= LU(ii,j) * invA(j,k);
      }
      invA(i,k) = sum;
    }
    //��������
    for (i=E::Rows-1; i>=0; i--) {
      ii = pivot[i];
      sum = invA(i,k);
      for (j=i+1; j<E::Rows; j++) {
	sum -= LU(ii,j) * invA(j,k);
      }
      invA(i,k) = sum / LU(ii,i);
    }
  }
  return det; //����ͤϹ���
}

// �չ���ʣ��ߣ��ѡ�
template <class E>
inline
double inv(const dMatrixExpression<E> &A, dMatrix<2,2> &invA) {
  invA(0,0) = A()(1,1);
  invA(0,1) =-A()(1,0);
  invA(1,0) =-A()(0,1);
  invA(1,1) = A()(0,0);
  
  double det = invA(0,0)*invA(1,1)-invA(0,1)*invA(1,0);
  if (det == 0) return 0;
  invA /= det;  return det; //����ͤϹ���
}

// �չ��󡧥��֥������ȤΥ��ԡ���������Τ��ѿ�������������� inv(A,invA) ��Ȥ��٤�
template <class E>
inline
dMatrix<E::Rows,E::Rows> inv(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> invA;
  
  inv(A, invA);
  
  return invA;
}

// ����λؿ��ؿ� expm(����,exp����)
template <class E>
inline
int expm(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &U){
  // exp(A)�η׻�
  // �ƥ��顼Ÿ���Ƿ׻�����ȸ��碪�ѥǶ�����Ѥ�����ˡ
  // exp(A) = { Pade(A/2^s) }^2s
  // ����:
  // Computing Matrix Exponentials using C++ with boost library
  // https://www.dbtsai.com/blog/2008-11-25-matrix-exponential/
  
  const int p = 6; // recommended and gererally satisfactory
  const dIdentityMatrix<E::Rows,E::Rows> I;
  
  int i, j, k;
  dMatrix<E::Rows,E::Rows> P, Q, A2;
  dMatrix<E::Rows,E::Rows> tmp;//Ʊ���ѿ��������Ǥ��ʤ��Τ�������ͳ��
  
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

// ����λؿ��ؿ������֥������ȤΥ��ԡ���������Τ��ѿ�������������� expm(A,eA) ��Ȥ��٤�
template <class E>
inline
dMatrix<E::Rows,E::Rows> expm(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> eA;
  
  expm(A, eA);
  
  return eA;
}


// ��������ؿ� signm(����,sign����)
template <class E>
inline
int signm(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &S){
  // Scaled Newton iteration
  // ����:
  // Ralph Byers, Hongguo Xu, "A New Scaling for Newton's Iteration for
  // the Polar Decomposition and its Backward Stability"
  
  dMatrix<E::Rows,E::Rows> Sk, Zk; S=A;
  double mu;
  
  do {
    Sk = S;
    
    //if (inv(Sk, Zk) == 0) return 0; //��®���Τ��ᥨ�顼������ʤ�
    inv(Sk, Zk);
    
    mu = sqrt(fnorm(Zk)/fnorm(Sk)); 
    //mu = pow((norm1(Zk)*infnorm(Zk))/(norm1(Sk)*infnorm(Sk)), 0.25); 
    //mu = 1.0/pow(fabs(det), 1.0/E::Rows); 
    //mu = 1; //scaling�ʤ�
    
    S = 0.5 * (mu*Sk + 1.0/mu*Zk);
  } while (fnorm(S-Sk) > DBL_EPSILON*fnorm(Sk));
  
  return 1;
}

// ��������ؿ������֥������ȤΥ��ԡ���������Τ��ѿ�������������� signm(A,S) ��Ȥ��٤�
template <class E>
inline
dMatrix<E::Rows,E::Rows> signm(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> S;
  
  signm(A, S);
  
  return S;
}


// �����ʿ���� sqrtm(����)
template <class E>
inline
int sqrtm(const dMatrixExpression<E> &A, dMatrix<E::Rows,E::Rows> &sA) {
  // ����:
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
  for (i=0; i<E::Rows; i++) { //��A �η׻��Ǥ�����
    F(i,i) += 1;
    F(i+E::Rows,i+E::Rows) -= 1;
  }
  */
  sA = submatrix(F, 0,E::Rows, E::Rows,E::Rows); // ��A = S12
  
  return 1;
}

// �����ʿ���������֥������ȤΥ��ԡ���������Τ��ѿ�������������� sqrtm(A,S) ��Ȥ��٤�
template <class E>
inline
dMatrix<E::Rows,E::Rows> sqrtm(const dMatrixExpression<E> &A) {
  dMatrix<E::Rows,E::Rows> sA;
  
  sqrtm(A, sA);
  
  return sA;
}


// QRʬ�� qr(����,Q����,R����) A = Q*R
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
      if (scale < tmp) scale = tmp; //�����ͤ����
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
      R(i,j) = 0; //�����Ѥ�0
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

// Hessenberg������Ѵ���subfunction for eig��
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
      A(i,j) = 0; // i>j+1��0��eig�η׻��Ǥ����ס�
    }
  }
  */
}

// ��ͭ�� eig (����,�����٥��ȥ�ei,�����٥��ȥ�ei)
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

// ɽ����
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

// ɽ����
template<class E>
inline
void print(const dMatrixExpression<E> &dm) {
  printf("% e", dm);
}
