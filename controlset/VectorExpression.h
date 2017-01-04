/* Time-stamp: <2013-02-26 18:20:20 kobayasi>
   
   �٥��ȥ�ι�®�黻�Τ����ExpressionTemplate��VectorClass.h�������ѡ�
   ver.2013-02-20  by K.Kobayashi
   
   fabs(B) �������Ǥ������ͤ�Ȥä��٥��ȥ���֤�
   round(B)�����������ͼθ�������
   ceil(B) ����������ŷ��ؿ�����
   floor(B)�������������ؿ�������   
   square(B)�� ��������衡������   
   
   subvector(B,i,N)���٥��ȥ��(i)������ʬ�٥��ȥ�<N>���֤�
                     ��B��ľ�ܽ񤭴�����ǽ��subvector(B,...)=X��
		     
   norm(B)   ���桼����åɥΥ��(L2�Υ��)���֤���= l2norm(B)��
   l1norm(B) ��L1�Υ����֤�
   l0norm(B) ��L0�Υ����֤�
   infnorm(B)��̵����Υ����֤�
   max(B)    �����Ǥκ����ͤ��֤�
   min(B)    �����ǤκǾ��ͤ��֤�
   sum(B)    �����Ǥ��¤��֤�
   
   print(B)������ printf��Ȥäƥ٥��ȥ��ɽ������
   printf(fmt,B)�������ͤ�ɽ���ե����ޥåȻ�����
*/


// �٥��ȥ��ɽ�����饹
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


// �黻�����Τ���ι�¤�Ρʹ���α黻�Ǥ����ѡ�
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


// �黻�����Ȥ˱����Ʊ黻��̤򥭥�å��夹�빽¤��
// �ʥ٥��ȥ�*����Τ褦��ʣ�����Ǥ��Ϥ�黻������Ȥ��˥���å��塧Costs=1��
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


// �٥��ȥ�α黻����ɽ�����饹�ʥ٥��ȥ�ȥ٥��ȥ�βø�����
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

// �٥��ȥ�α黻����ɽ�����饹�ʥ٥��ȥ�ȿ��ͤξ������
template<class L, class Ope>
class dVectorOperation<L, Ope, double>:
  public dVectorExpression< dVectorOperation<L, Ope, double> > {
public:
  static const bool Costs = L::Costs;
  static const int Size = L::Size;
private:
  const L &lv_;
  const double rd_;//  __attribute__((aligned(8))); //®���ʤ롩
public:
  dVectorOperation(const L &lv, const double rd): lv_(lv), rd_(rd) {}
  
  double operator[] (int i) const { return Ope::apply(lv_[i], rd_); }
  double operator() (int i) const { return Ope::apply(lv_(i), rd_); }
};

// �٥��ȥ�α黻����ɽ�����饹�ʥ٥��ȥ����黻�ʤɡ�
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

// �û��ʥ٥��ȥ� + �٥��ȥ��
template<class L, class R>
inline
dVectorOperation<L, OperationAdd, R> operator + (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  return dVectorOperation<L, OperationAdd, R>( lv(), rv() );
}
// �����ʥ٥��ȥ� - �٥��ȥ��
template<class L, class R>
inline
dVectorOperation<L, OperationSub, R> operator - (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  return dVectorOperation<L, OperationSub, R>( lv(), rv() );
}
// �軻�ʥ٥��ȥ� * ���͡�
template<class L>
inline
dVectorOperation<L, OperationMul, double> operator * (const dVectorExpression<L> &lv, const double rd) {
  return dVectorOperation<L, OperationMul, double>( lv(), rd );
}
// �軻�ʿ��� * �٥��ȥ�ˡ�ʥ٥��ȥ� * ���͡�
template<class R>
inline
dVectorOperation<R, OperationMul, double> operator * (const double ld, const dVectorExpression<R> &rv) {
  return dVectorOperation<R, OperationMul, double>( rv(), ld );
}
// �����ʥ٥��ȥ� / ���͡�
template<class L>
inline
dVectorOperation<L, OperationDiv, double> operator / (const dVectorExpression<L> &lv, const double rd) {
  return dVectorOperation<L, OperationDiv, double>( lv(), rd );
}
// ���- �٥��ȥ��
template<class R>
inline
dVectorOperationS<OperationNeg, R> operator - (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationNeg, R>( rv() );
}
// ����ʥ٥��ȥ� == �٥��ȥ��// ���Ǥ�double���ʤΤǱ黻��̤���Ӥ���ΤϤ������ᤷ�ʤ�
template<class L, class R>
inline
bool operator == (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  for (int i = 0; i < L::Size; i++) {
    if ( lv()[i] != rv()[i] ) return 0;
  }
  return 1;
}
// ������ʥ٥��ȥ� != �٥��ȥ��
template<class L, class R>
inline
bool operator != (const dVectorExpression<L> &lv, const dVectorExpression<R> &rv) {
  for (int i = 0; i < L::Size; i++) {
    if ( lv()[i] != rv()[i] ) return 1;
  }
  return 0;
}

// ���� fabs(�٥��ȥ�)
template<class R>
inline
dVectorOperationS<OperationFabs, R> fabs (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationFabs, R>( rv() );
}
// �ͼθ��� round(�٥��ȥ�)
template<class R>
inline
dVectorOperationS<OperationRound, R> round (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationRound, R>( rv() );
}
// ŷ��ؿ� ceil(�٥��ȥ�)
template<class R>
inline
dVectorOperationS<OperationCeil, R> ceil (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationCeil, R>( rv() );
}
// ���ؿ� floor(�٥��ȥ�)
template<class R>
inline
dVectorOperationS<OperationFloor, R> floor (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationFloor, R>( rv() );
}
// �����Ǥ���� square(�٥��ȥ�)
template<class R>
inline
dVectorOperationS<OperationSquare, R> square (const dVectorExpression<R> &rv) {
  return dVectorOperationS<OperationSquare, R>( rv() );
}

// ���ѡʥ٥��ȥ� * �٥��ȥ��
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


// �桼����åɥΥ���L2�Υ��� norm(�٥��ȥ�)
template<class R>
inline
double norm (const dVectorExpression<R> &rv) {
  return sqrt( rv * rv );
}

#define l2norm(rv)  norm(rv)

// L1�Υ�� l1norm(�٥��ȥ�)
template<class R>
inline
double l1norm (const dVectorExpression<R> &rv) {
  return sum( fabs(rv) );
}

// L0�Υ�� l0norm(�٥��ȥ�)
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

// ̵����Υ�� infnorm(�٥��ȥ�)
template<class R>
inline
double infnorm (const dVectorExpression<R> &rv) {
  return max( fabs(rv) );
}

// ������ max(�٥��ȥ�)
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

// �Ǿ��� min(�٥��ȥ�)
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

// ���Ǥ��� sum(�٥��ȥ�)
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

// ɽ����
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

// ɽ����
template<class E>
inline
void print(const dVectorExpression<E> &dv) {
  printf("% e", dv);
}
